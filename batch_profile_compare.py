"""
Batch Profile Composite Analysis
---------------------------------
Loops over all energy folders/files, runs analysis on Measurement vs MC sheets
(or any two named sheets), saves one figure per energy + scan-axis and a
summary CSV of pass-rate statistics.

Parameters are set in the CONFIG section below.
"""

import os
import io
import glob
import re
import contextlib
import platform
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # non-interactive backend — no GUI needed
from matplotlib import pyplot as plt
from scipy import interpolate as interp
from scipy.ndimage import gaussian_filter1d

from comp import dosedif, dta as dtafunc
from center import center as center_function

# ─────────────────────────────────────────────
#  CONFIG  — edit these as needed
# ─────────────────────────────────────────────

FILE_MODE = 'flat'
# 'hierarchical' : BASE_PATH contains energy subfolders (e.g. 6X/, 10FFF/),
#                  each holding one *Profile*.xlsx file.
# 'flat'         : BASE_PATH is a single folder of xlsx files; energy label is
#                  parsed from each filename.

BASE_PATH = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\Combined Consortium Data\Processed Combined Data"

FILE_FILTER = '*Profile*.xlsx'    # glob used in 'flat' mode only

SHEET1_NAME = "TPS SN 20"             # reference / measured sheet name
SHEET2_NAME = "SN 20"         # comparison sheet name

# Data selection — 'all' means use all common values found in each file
FIELD_SIZES = 'all'                # e.g. [5.0, 10.0, 20.0]  or  'all'
AXES        = 'all'                # e.g. ['X', 'Y']          or  'all'  (Z excluded always)
DEPTHS_CM   = 'all'                # e.g. [1.5, 5.0, 10.0]   or  'all'

DEPTH_ROUND_CM = 0.1               # rounding resolution for depth matching [cm]

ANALYSIS      = 'mppg'             # 'comp', 'dif', 'dist', 'plot', 'gam', 'mppg'
DD_CRITERIA   = 2.0                # dose difference threshold [%]
DTA_CRITERIA  = 0.2                # DTA threshold [cm]  (2 mm)

# Normalization
NORM = 1
# 0 = no normalization (raw dose units)
# 1 = normalize each profile to its CAX  (Off Axis Ratio)
# 2 = normalize each profile to its own Dmax
PLOT_PDD_SCALE = True
CAX_WINDOW_MM = 0
# 0 = use nearest point to axis (recommended for small fields / sparse data)
# N = average within ±N mm of axis
SCALE_TARGET = 0       # 0=none, 1=scale profile 1, 2=scale profile 2, 3=scale both
SCALE_FACTOR = 1.00    # multiplier applied after normalization
SAVE_REPORT  = True    # generate multi-page PDF report (summary + per-energy figures)
REPORT_DPI   = 150     # DPI for PDF report pages (lower = smaller file, faster to open)
# True  = scale the plotted curves by PDD(depth) for visual absolute-dose display
# False = plot in the same units as the normalization above

# Profile centering
CENTER = 3
# 0 = none
# 1 = shift curve 1 to geometric center
# 2 = shift curve 2 to geometric center
# 3 = shift both curves to geometric center
CENTER_THRESHOLD = 0.3             # threshold fraction passed to center()

# Detector convolution
CONV_FWHM_CM = 0.0                 # detector FWHM [cm]; 0 = disabled
CONV_TARGET  = 'none'              # 'none', 'curve1', 'curve2', or 'both'

# MPPG-specific parameters (used only when ANALYSIS == 'mppg')
MPPG_DDTAIL  = 2.0                 # tail dose-difference criterion [% of Dmax]
MPPG_PEN_CM  = 0.25                # penumbra half-width [cm]  (GUI default: pupper_entry="0.5" / 2)
MPPG_OVR_CM  = 1.0                 # overlap buffer zone [cm]  (GUI default: pulower_entry="1")
MPPG_DIAG_FACTOR = 0.80            # diagonal field size factor (GUI default: 0.80)

# XY diagonal scans: only process for the largest field size (matches GUI behavior)
XY_LARGEST_FS_ONLY = True

MARKER_SIZE  = 4                   # plot marker size
DPI          = 600                 # figure output resolution
SAVE_FIGURES = True                # set False to skip figure generation (faster, stats only)

_s1 = SHEET1_NAME.replace(' ', '')
_s2 = SHEET2_NAME.replace(' ', '')

# Build a criteria tag from analysis type and thresholds so runs with different
# criteria land in separate folders automatically (e.g. mppg_2pct_2mm).
_dd  = f"{DD_CRITERIA:.4g}".rstrip('0').rstrip('.')
_dta = f"{DTA_CRITERIA*10:.4g}".rstrip('0').rstrip('.')
_criteria_tag = f"{ANALYSIS}_{_dd}pct_{_dta}mm"

# Output folder hierarchy:
#   Results/
#     <SN1>_vs_<SN2>/
#       Profile/
#         <criteria>/
#           <summary>.xlsx
#           <summary>.pdf          ← all-energy report
#           Individual Figures/    ← per-energy PNGs
#           Individual Reports/    ← per-energy PDFs
COMPARISON_DIR = os.path.join(BASE_PATH, "Results", f"{_s1}_vs_{_s2}")
PROFILE_DIR    = os.path.join(COMPARISON_DIR, "Profile", _criteria_tag)
RESULTS_DIR    = os.path.join(PROFILE_DIR, "Individual Figures")
REPORTS_DIR    = os.path.join(PROFILE_DIR, "Individual Reports")
# ─────────────────────────────────────────────


# ── PDD lookup table (for NORM=2 and MPPG tail scaling) ──────────────────────
PDD_TABLE = {
    "6X":   {1.42: 1.000,  5.0: 0.865, 10.0: 0.668, 20.0: 0.382, 30.0: 0.218},
    "6FFF": {1.32: 1.000,  5.0: 0.848, 10.0: 0.636, 20.0: 0.346, 30.0: 0.190},
    "8FFF": {1.88: 1.000,  5.0: 0.891, 10.0: 0.693, 20.0: 0.406, 30.0: 0.240},
    "10X":  {2.22: 1.000,  5.0: 0.925, 10.0: 0.744, 20.0: 0.469, 30.0: 0.296},
    "10FFF":{2.16: 1.000,  5.0: 0.904, 10.0: 0.708, 20.0: 0.423, 30.0: 0.256},
    "15X":  {2.72: 1.000,  5.0: 0.949, 10.0: 0.774, 20.0: 0.502, 30.0: 0.325},
}


def pdd_lookup_nearest(energy, depth_cm):
    table = PDD_TABLE.get(energy, {})
    if not table:
        return 1.0, None
    keys = np.array(list(table.keys()), dtype=float)
    k = float(keys[np.argmin(np.abs(keys - float(depth_cm)))])
    return float(table[k]), k


# ── helpers ───────────────────────────────────────────────────────────────────

def parse_ssd_energy(path, default_ssd=100.0, default_energy="6X"):
    """Extract SSD (cm) and energy token from file path."""
    s = str(path)
    m_ssd = (re.search(r'ssd[_\-\s]?(\d+(?:\.\d+)?)', s, re.I)
             or re.search(r'(\d+(?:\.\d+)?)\s*ssd', s, re.I))
    ssd_val = float(m_ssd.group(1)) if m_ssd else float(default_ssd)
    s_norm = re.sub(r'[^A-Za-z0-9.]+', ' ', s)
    m_energy = re.search(r'(?:^|\s)(6X|6FFF|8FFF|10X|10FFF|15X)(?=\s|$)', s_norm, re.I)
    energy_val = m_energy.group(1).upper() if m_energy else default_energy
    return ssd_val, energy_val


def apply_detector_convolution(x, dose, fwhm_cm):
    """Gaussian blur along the lateral axis to simulate detector volume averaging."""
    if fwhm_cm <= 0:
        return dose
    sigma_cm  = fwhm_cm / 2.355
    step_cm   = np.mean(np.diff(np.asarray(x, dtype=float)))
    if step_cm <= 0:
        return dose
    sigma_idx = sigma_cm / step_cm
    return gaussian_filter1d(np.asarray(dose, dtype=float), sigma=sigma_idx)


def downsample_to_native(x_interp, vals_interp, x_native):
    """Thin interpolated output back to roughly the native data spacing."""
    x_interp    = np.asarray(x_interp,    dtype=float)
    vals_interp = np.asarray(vals_interp, dtype=float)

    if isinstance(x_native, (list, tuple)):
        natives = [np.asarray(x, dtype=float) for x in x_native if x is not None]
    else:
        natives = [np.asarray(x_native, dtype=float)]

    if x_interp.size < 2 or not natives:
        return x_interp, vals_interp

    def _median_step(x):
        x = x[np.isfinite(x)]
        if x.size < 2:
            return 0.0
        xu = np.unique(x)
        xu.sort()
        diffs = np.diff(xu)
        diffs = diffs[diffs > 0]
        return float(np.median(diffs)) if diffs.size else 0.0

    steps = [s for s in (_median_step(x) for x in natives) if s > 0]
    if not steps:
        return x_interp, vals_interp

    native_step = max(steps)
    xi = np.unique(x_interp[np.isfinite(x_interp)])
    xi.sort()
    if xi.size < 2:
        return x_interp, vals_interp
    interp_step = float(np.median(np.diff(xi)))
    if interp_step <= 0:
        return x_interp, vals_interp

    factor = max(1, int(round(native_step / interp_step)))
    return x_interp[::factor], vals_interp[::factor]


# ── analysis for one file ─────────────────────────────────────────────────────

def run_one_file(xlsx_path, energy_label):
    """
    Run profile analysis on one Excel file.
    Returns (results, skipped) where results is a list of result dicts
    and skipped is a list of warning strings for flagged/skipped combos.
    """
    print(f"\n{'='*60}")
    print(f"  {energy_label}  —  {os.path.basename(xlsx_path)}")
    print(f"{'='*60}")

    # Parse SSD and energy key from filename (energy_label may be a full filename stem)
    ssd_cm, pdd_energy = parse_ssd_energy(xlsx_path, default_ssd=100.0, default_energy=energy_label)
    print(f"  SSD parsed from path: {ssd_cm:.1f} cm  |  PDD energy key: {pdd_energy}")

    df1 = pd.read_excel(xlsx_path, sheet_name=SHEET1_NAME, header=0)
    df2 = pd.read_excel(xlsx_path, sheet_name=SHEET2_NAME, header=0)

    # Standardize depth rounding to avoid floating-point mismatches
    for df in (df1, df2):
        df["Depth"] = ((df["Depth"] / DEPTH_ROUND_CM).round() * DEPTH_ROUND_CM).round(3)

    df1 = df1.sort_values(by=['FS', 'Depth', 'Axis', 'Pos'])
    df2 = df2.sort_values(by=['FS', 'Depth', 'Axis', 'Pos'])

    # ── determine common FS / Axis / Depth (always exclude Z-axis PDDs) ───────
    def _common(col):
        s1 = set(df1.loc[df1['Axis'] != 'Z', col].dropna().unique())
        s2 = set(df2.loc[df2['Axis'] != 'Z', col].dropna().unique())
        return sorted(s1 & s2, key=float)

    fsl_all  = _common('FS')
    dl_all   = _common('Depth')
    axes_all = sorted(
        set(df1.loc[df1['Axis'] != 'Z', 'Axis'].dropna().unique()) &
        set(df2.loc[df2['Axis'] != 'Z', 'Axis'].dropna().unique()),
        key=lambda x: (x in ('XY', 'YX'), x)   # XY/YX processed last
    )

    # Apply config filters
    fsl  = [f for f in fsl_all  if FIELD_SIZES == 'all' or float(f) in [float(x) for x in FIELD_SIZES]]
    axes = [a for a in axes_all if AXES        == 'all' or a in AXES]
    dl   = [d for d in dl_all   if DEPTHS_CM   == 'all' or float(d) in [float(x) for x in DEPTHS_CM]]

    if not fsl or not axes or not dl:
        print("  No common FS / Axis / Depth after applying filters — skipping.")
        return [], [], None

    # Pre-computed thresholds in working units
    dd_frac     = DD_CRITERIA  / 100.0   # fraction  (e.g. 0.03 for 3%)
    ddtail_frac = MPPG_DDTAIL  / 100.0   # fraction of Dmax
    dta_cm      = DTA_CRITERIA            # cm

    stem        = os.path.splitext(os.path.basename(xlsx_path))[0]
    all_results = []
    skipped     = []

    # ── profile panel y-label reflects normalization mode ─────────────────
    profile_ylabel = ('Off Axis Ratio'            if NORM == 1
                      else 'Dose (norm. to Dmax)' if NORM == 2
                      else 'Dose (raw units)')

    # ── single figure for all axes overlaid (matches GUI layout) ──────────
    fig = None
    if SAVE_FIGURES or SAVE_REPORT:
        if ANALYSIS == 'mppg':
            plt.rcParams.update({'font.size': 22})
            fig, (ax0, ax1, ax2, ax3) = plt.subplots(
                4, 1, figsize=(15, 14),
                gridspec_kw={'height_ratios': [1.5, 1, 1, 1]}
            )
            ax1.set_ylabel('DTA [mm]')
            ax2.set_ylabel('$\\Delta$Dose [%]')
            ax3.set_ylabel('$\\Delta$Dose [%Dmax]')
        elif ANALYSIS in ('comp', 'gam'):
            plt.rcParams.update({'font.size': 20})
            fig, (ax0, ax1, ax2) = plt.subplots(
                3, 1, figsize=(15, 11),
                gridspec_kw={'height_ratios': [1.5, 1, 1]}
            )
            if ANALYSIS == 'comp':
                ax1.set_ylabel('DTA [mm]')
                ax2.set_ylabel('Dose Difference [%]')
        elif ANALYSIS in ('dif', 'dist'):
            plt.rcParams.update({'font.size': 18})
            fig, (ax0, ax1) = plt.subplots(
                2, 1, figsize=(15, 9),
                gridspec_kw={'height_ratios': [1.5, 1]}
            )
            ax1.set_ylabel('Dose Difference [%]' if ANALYSIS == 'dif' else 'DTA [mm]')
        else:   # 'plot'
            plt.rcParams.update({'font.size': 16})
            fig, ax0 = plt.subplots(1, 1, figsize=(15, 8))
        ax0.plot([], '+r', ms=10, label=SHEET1_NAME)
        ax0.plot([], '.k', ms=10, label=SHEET2_NAME)
        ax0.legend()
    else:
        class _NullAx:
            def plot(self, *_, **__): pass
            def set_xlim(self, *_, **__): pass
            def set_ylim(self, *_, **__): pass
            def set_ylabel(self, *_, **__): pass
            def set_xlabel(self, *_, **__): pass
            def legend(self, *_, **__): pass
            def hist(self, *_, **__): pass
        fig = None
        ax0 = ax1 = ax2 = ax3 = _NullAx()
    ax0.set_ylabel(profile_ylabel)
    ax0.set_xlabel('Off Axis Position [cm]')

    # Per-file pass/fail accumulators (all axes combined)
    total_pts_axis  = 0
    total_fail_axis = 0
    gvtot = []

    # ── outer loop: all axes overlaid onto the same figure ────────────────
    for axis in axes:
        print(f"\n  ── Axis: {axis} ──")

        # ── loops: field size → depth ──────────────────────────────────────────
        for fs_val in fsl:
            fs_f = float(fs_val)

            # XY / YX diagonal scans: skip if not the largest field size
            if XY_LARGEST_FS_ONLY and axis in ('XY', 'YX') and fs_f != float(max(fsl)):
                continue

            for depth in dl:
                mask1 = (
                    (df1['FS']    == fs_val) &
                    (df1['Axis']  == axis)   &
                    (df1['Depth'] == depth)
                )
                mask2 = (
                    (df2['FS']    == fs_val) &
                    (df2['Axis']  == axis)   &
                    (df2['Depth'] == depth)
                )

                y1 = df1.loc[mask1, 'Pos'].reset_index(drop=True)
                d1 = df1.loc[mask1, 'Dose'].reset_index(drop=True)
                y2 = df2.loc[mask2, 'Pos'].reset_index(drop=True)
                d2 = df2.loc[mask2, 'Dose'].reset_index(drop=True)

                if len(y1) == 0 or len(y2) == 0:
                    continue   # no data for this combo — silent skip
                if len(y1) < 3 or len(y2) < 3:
                    msg = f"{energy_label}  FS {fs_val} depth {depth} [{axis}]: only {min(len(y1),len(y2))} point(s)"
                    print(f"    WARNING — {msg} — skipping.")
                    skipped.append(msg)
                    continue

                # Sort by position (should already be sorted, but enforce it)
                s1 = y1.argsort()
                y1, d1 = y1.iloc[s1].reset_index(drop=True), d1.iloc[s1].reset_index(drop=True)
                s2 = y2.argsort()
                y2, d2 = y2.iloc[s2].reset_index(drop=True), d2.iloc[s2].reset_index(drop=True)

                # Check for duplicate positions (would crash pchip interpolation)
                dup1 = y1[y1.duplicated(keep=False)]
                dup2 = y2[y2.duplicated(keep=False)]
                if not dup1.empty:
                    msg = (f"{energy_label}  FS {fs_val} depth {depth} [{axis}]: "
                           f"{SHEET1_NAME} duplicate Pos {sorted(dup1.unique().tolist())}")
                    print(f"    WARNING — {msg} — skipping.")
                    skipped.append(msg)
                    continue
                if not dup2.empty:
                    msg = (f"{energy_label}  FS {fs_val} depth {depth} [{axis}]: "
                           f"{SHEET2_NAME} duplicate Pos {sorted(dup2.unique().tolist())}")
                    print(f"    WARNING — {msg} — skipping.")
                    skipped.append(msg)
                    continue

                # ── centering ─────────────────────────────────────────────────
                # center() requires a pandas Series with a clean integer index (already done above)
                if CENTER in (1, 3):
                    try:
                        shift1 = center_function(y1, d1, CENTER_THRESHOLD)
                        y1 = y1 + shift1
                    except Exception as e:
                        print(f"    FS {fs_val} depth {depth}: center() curve1 failed ({e})")
                if CENTER in (2, 3):
                    try:
                        shift2 = center_function(y2, d2, CENTER_THRESHOLD)
                        y2 = y2 + shift2
                    except Exception as e:
                        print(f"    FS {fs_val} depth {depth}: center() curve2 failed ({e})")

                # ── detector convolution ───────────────────────────────────────
                if CONV_TARGET in ('curve1', 'both'):
                    d1 = pd.Series(apply_detector_convolution(y1, d1, CONV_FWHM_CM))
                if CONV_TARGET in ('curve2', 'both'):
                    d2 = pd.Series(apply_detector_convolution(y2, d2, CONV_FWHM_CM))

                # ── normalization ──────────────────────────────────────────────
                # d1_analysis/d2_analysis are always PDD-free so 1% = 1% of CAX
                d1_analysis = d1; d2_analysis = d2
                if NORM == 1:
                    if CAX_WINDOW_MM == 0:
                        cax1 = float(d1.iloc[np.abs(np.asarray(y1)).argmin()])
                        cax2 = float(d2.iloc[np.abs(np.asarray(y2)).argmin()])
                    else:
                        cax1 = float(d1[np.abs(np.asarray(y1)) <= CAX_WINDOW_MM / 10.0].mean())
                        cax2 = float(d2[np.abs(np.asarray(y2)) <= CAX_WINDOW_MM / 10.0].mean())
                    if cax1 > 0:
                        d1_analysis = d1 / cax1
                    if cax2 > 0:
                        d2_analysis = d2 / cax2
                elif NORM == 2:
                    d1_analysis = d1 / d1.max(); d2_analysis = d2 / d2.max()
                # profile scaling (after norm, independent of normalization)
                if SCALE_TARGET in (1, 3):
                    d1_analysis = d1_analysis * SCALE_FACTOR
                if SCALE_TARGET in (2, 3):
                    d2_analysis = d2_analysis * SCALE_FACTOR
                # apply PDD scaling for plotting only
                if PLOT_PDD_SCALE:
                    pdd_factor, _ = pdd_lookup_nearest(pdd_energy, depth)
                    d1 = d1_analysis * pdd_factor
                    d2 = d2_analysis * pdd_factor
                else:
                    d1 = d1_analysis; d2 = d2_analysis

                # ── plot the curves ────────────────────────────────────────────
                ax0.plot(y1, d1, '+r', ms=MARKER_SIZE)
                ax0.plot(y2, d2, '.k', ms=MARKER_SIZE)

                if ANALYSIS == 'plot':
                    continue

                # ── gamma analysis ─────────────────────────────────────────────
                elif ANALYSIS == 'gam':
                    from gamma import gamma as g
                    gx, gv = g(y1, d1_analysis, y2, d2_analysis, dd_frac, dta_cm, 1, 0.01, 0.01)
                    gx, gv = downsample_to_native(gx, gv, y1)
                    gv_a   = np.asarray(gv)
                    valid  = np.isfinite(gv_a)
                    total  = int(valid.sum())
                    passed = int(np.count_nonzero(gv_a[valid] <= 1.0))
                    failed = total - passed
                    passrate = passed / total * 100 if total > 0 else 0.0
                    gvmean   = float(np.mean(gv_a[valid])) if total > 0 else float('nan')
                    total_pts_axis  += total
                    total_fail_axis += failed
                    gvtot.extend(gv_a[valid].tolist())

                    ax1.plot(gx[gv_a >  1], gv_a[gv_a >  1], '.r', ms=MARKER_SIZE)
                    ax1.plot(gx[gv_a <= 1], gv_a[gv_a <= 1], '.g', ms=MARKER_SIZE)

                    all_results.append({
                        'Energy': energy_label, 'FS': fs_val, 'Axis': axis,
                        'Depth_cm': depth, 'PassRate_pct': round(passrate, 2),
                        'MeanGamma': round(gvmean, 3),
                        'TotalPoints': total, 'FailPoints': failed,
                    })
                    print(f"    FS {fs_val} depth {depth:.2f}: Pass={passrate:.2f}%  "
                          f"MeanGamma={gvmean:.3f}  Fail={failed}/{total}")

                # ── composite analysis ─────────────────────────────────────────
                elif ANALYSIS == 'comp':
                    dtax, dtav = dtafunc(y1, d1_analysis, y2, d2_analysis, dta_cm)
                    difx, difv = dosedif(y1, d1_analysis, y2, d2_analysis, 0)
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    difx, difv = downsample_to_native(difx, difv, y1)
                    dtafail  = np.where(np.abs(dtav) > dta_cm)[0]
                    ddiffail = np.where(np.abs(difv) > dd_frac)[0]
                    fail_set = set(dtafail) & set(ddiffail)
                    total    = len(dtav)
                    failed   = len(fail_set)
                    passrate = (total - failed) / total * 100 if total > 0 else 0.0
                    mean_dd  = float(np.mean(np.abs(difv))) * 100 if total > 0 else 0.0
                    mean_dta = float(np.mean(np.abs(dtav)))       if total > 0 else 0.0
                    total_pts_axis  += total
                    total_fail_axis += failed

                    ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
                    ax2.plot(difx, difv * 100, '.g', ms=0.5)
                    for n in fail_set:
                        ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)
                        ax2.plot(difx[n], difv[n] * 100, '.r', ms=0.5)

                    all_results.append({
                        'Energy': energy_label, 'FS': fs_val, 'Axis': axis,
                        'Depth_cm': depth, 'PassRate_pct': round(passrate, 2),
                        'MeanDoseDiff_pct': round(mean_dd, 3),
                        'MeanDTA_mm': round(mean_dta * 10, 3),
                        'TotalPoints': total, 'FailPoints': failed,
                    })
                    print(f"    FS {fs_val} depth {depth:.2f}: Pass={passrate:.2f}%  "
                          f"MeanDD={mean_dd:.2f}%  MeanDTA={mean_dta*10:.2f}mm  Fail={failed}/{total}")

                # ── dose difference analysis ───────────────────────────────────
                elif ANALYSIS == 'dif':
                    difx, difv = dosedif(y1, d1_analysis, y2, d2_analysis, 0)
                    difx, difv = downsample_to_native(difx, difv, y1)
                    ddiffail = np.where(np.abs(difv) > dd_frac)[0]
                    total    = len(difv)
                    failed   = len(ddiffail)
                    passrate = (total - failed) / total * 100 if total > 0 else 0.0
                    mean_dd  = float(np.mean(np.abs(difv))) * 100 if total > 0 else 0.0
                    total_pts_axis  += total
                    total_fail_axis += failed

                    ax1.plot(difx, difv * 100, '.g', ms=0.5)
                    for n in ddiffail:
                        ax1.plot(difx[n], difv[n] * 100, '.r', ms=0.5)

                    all_results.append({
                        'Energy': energy_label, 'FS': fs_val, 'Axis': axis,
                        'Depth_cm': depth, 'PassRate_pct': round(passrate, 2),
                        'MeanDoseDiff_pct': round(mean_dd, 3),
                        'TotalPoints': total, 'FailPoints': failed,
                    })
                    print(f"    FS {fs_val} depth {depth:.2f}: Pass={passrate:.2f}%  "
                          f"MeanDD={mean_dd:.2f}%  Fail={failed}/{total}")

                # ── DTA analysis ───────────────────────────────────────────────
                elif ANALYSIS == 'dist':
                    dtax, dtav = dtafunc(y1, d1_analysis, y2, d2_analysis, dta_cm)
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    dtafail  = np.where(np.abs(dtav) > dta_cm)[0]
                    total    = len(dtav)
                    failed   = len(dtafail)
                    passrate = (total - failed) / total * 100 if total > 0 else 0.0
                    mean_dta = float(np.mean(np.abs(dtav))) if total > 0 else 0.0
                    total_pts_axis  += total
                    total_fail_axis += failed

                    ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
                    for n in dtafail:
                        ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)

                    all_results.append({
                        'Energy': energy_label, 'FS': fs_val, 'Axis': axis,
                        'Depth_cm': depth, 'PassRate_pct': round(passrate, 2),
                        'MeanDTA_mm': round(mean_dta * 10, 3),
                        'TotalPoints': total, 'FailPoints': failed,
                    })
                    print(f"    FS {fs_val} depth {depth:.2f}: Pass={passrate:.2f}%  "
                          f"MeanDTA={mean_dta*10:.2f}mm  Fail={failed}/{total}")

                # ── MPPG analysis ──────────────────────────────────────────────
                elif ANALYSIS == 'mppg':
                    pdd_factor, _ = pdd_lookup_nearest(pdd_energy, depth)
                    diag    = np.sqrt(2.0 * MPPG_DIAG_FACTOR) if axis in ('XY', 'YX') else 1.0
                    edge    = (fs_f / 2.0) * diag * (ssd_cm + depth) / 100.0

                    dtax, dtav = dtafunc(y1, d1_analysis, y2, d2_analysis, dta_cm)
                    difx, difv = dosedif(y1, d1_analysis, y2, d2_analysis, 0)
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    difx, difv = downsample_to_native(difx, difv, y1)
                    difv_dmax  = np.asarray(difv) * pdd_factor

                    x_dif = np.asarray(difx)
                    x_dta = np.asarray(dtax)
                    PEN   = MPPG_PEN_CM
                    OVR   = MPPG_OVR_CM
                    left  = -edge
                    right =  edge

                    # Region masks on the ΔDose grid
                    pen_left_x   = (x_dif >= left  - PEN) & (x_dif <= left  + PEN)
                    pen_right_x  = (x_dif >= right - PEN) & (x_dif <= right + PEN)
                    pen_core_x   = pen_left_x  | pen_right_x
                    in_core_x    = (x_dif >  left  + (PEN + OVR)) & (x_dif < right - (PEN + OVR))
                    tail_core_x  = (x_dif <= left  - (PEN + OVR)) | (x_dif >= right + (PEN + OVR))
                    near_hi_x    = (
                        ((x_dif >  left  + PEN) & (x_dif <= left  + PEN + OVR)) |
                        ((x_dif >= right - PEN - OVR) & (x_dif <  right - PEN))
                    )
                    near_lo_x    = (
                        ((x_dif >= left  - PEN - OVR) & (x_dif < left  - PEN)) |
                        ((x_dif >  right + PEN)        & (x_dif <= right + PEN + OVR))
                    )

                    # Region masks on the DTA grid (penumbra)
                    pen_core_dta = (
                        ((x_dta >= left  - PEN) & (x_dta <= left  + PEN)) |
                        ((x_dta >= right - PEN) & (x_dta <= right + PEN))
                    )

                    dtav_a  = np.asarray(dtav)
                    difv_a  = np.asarray(difv)
                    pass_dta    = np.abs(dtav_a)    <= dta_cm
                    pass_dd     = np.abs(difv_a)    <= dd_frac
                    pass_ddtail = np.abs(difv_dmax) <= ddtail_frac

                    # Interpolate DTA onto ΔDose grid for overlap OR-rule
                    if len(dtax) > 3:
                        dta_on_difx      = interp.pchip(dtax, dtav)(difx)
                        pass_dta_on_difx = np.abs(dta_on_difx) <= dta_cm
                    else:
                        pass_dta_on_difx = np.zeros(len(difx), dtype=bool)

                    overlap_hi_pass = near_hi_x & (pass_dd | pass_dta_on_difx)
                    overlap_hi_fail = near_hi_x & ~(pass_dd | pass_dta_on_difx)
                    overlap_lo_pass = near_lo_x & (pass_ddtail | pass_dta_on_difx)
                    overlap_lo_fail = near_lo_x & ~(pass_ddtail | pass_dta_on_difx)

                    domain_mask = in_core_x | tail_core_x | pen_core_x | near_hi_x | near_lo_x
                    pass_any    = (
                        (in_core_x   & pass_dd)          |
                        (tail_core_x & pass_ddtail)       |
                        (pen_core_x  & pass_dta_on_difx)  |
                        overlap_hi_pass                   |
                        overlap_lo_pass
                    )
                    total    = int(domain_mask.sum())
                    passed   = int(pass_any.sum())
                    failed   = total - passed
                    passrate = passed / total * 100 if total > 0 else 0.0
                    mean_dd  = float(np.mean(np.abs(difv_a))) * 100 if total > 0 else 0.0
                    total_pts_axis  += total
                    total_fail_axis += failed

                    # Plotting — ax1 (DTA penumbra)
                    ax1.plot(dtax, dtav_a * 10, '.k', ms=MARKER_SIZE)
                    ax1.plot(dtax[pen_core_dta &  pass_dta],
                             dtav_a[pen_core_dta &  pass_dta] * 10, '.g', ms=MARKER_SIZE)
                    ax1.plot(dtax[pen_core_dta & ~pass_dta],
                             dtav_a[pen_core_dta & ~pass_dta] * 10, '.r', ms=MARKER_SIZE)

                    # ax2 (ΔDose in-field + overlap_hi)
                    ax2.plot(difx, difv_a * 100, '.k', ms=MARKER_SIZE)
                    ax2.plot(difx[in_core_x &  pass_dd],
                             difv_a[in_core_x &  pass_dd] * 100, '.g', ms=MARKER_SIZE)
                    ax2.plot(difx[in_core_x & ~pass_dd],
                             difv_a[in_core_x & ~pass_dd] * 100, '.r', ms=MARKER_SIZE)
                    ax2.plot(difx[overlap_hi_pass], difv_a[overlap_hi_pass] * 100, '.g', ms=MARKER_SIZE)
                    ax2.plot(difx[overlap_hi_fail], difv_a[overlap_hi_fail] * 100, '.r', ms=MARKER_SIZE)

                    # ax3 (ΔDose×PDD tails + overlap_lo)
                    ax3.plot(difx, difv_dmax * 100, '.k', ms=MARKER_SIZE)
                    ax3.plot(difx[tail_core_x &  pass_ddtail],
                             difv_dmax[tail_core_x &  pass_ddtail] * 100, '.g', ms=MARKER_SIZE)
                    ax3.plot(difx[tail_core_x & ~pass_ddtail],
                             difv_dmax[tail_core_x & ~pass_ddtail] * 100, '.r', ms=MARKER_SIZE)
                    ax3.plot(difx[overlap_lo_pass], difv_dmax[overlap_lo_pass] * 100, '.g', ms=MARKER_SIZE)
                    ax3.plot(difx[overlap_lo_fail], difv_dmax[overlap_lo_fail] * 100, '.r', ms=MARKER_SIZE)

                    all_results.append({
                        'Energy': energy_label, 'FS': fs_val, 'Axis': axis,
                        'Depth_cm': depth, 'PassRate_pct': round(passrate, 2),
                        'MeanDoseDiff_pct': round(mean_dd, 3),
                        'TotalPoints': total, 'FailPoints': failed,
                    })
                    print(f"    FS {fs_val} depth {depth:.2f}: MPPG Pass={passrate:.2f}%  Fail={failed}/{total}")

    # ── figure title and pass-rate label (after all axes plotted) ────────
    overall_pass = (
        100.0 * (1 - total_fail_axis / total_pts_axis)
        if total_pts_axis > 0 else 0.0
    )

    if ANALYSIS == 'comp':
        title_tag = f'Composite {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm'
        ax2.set_xlabel(
            f'Points outside {DD_CRITERIA:.1f}% & {DTA_CRITERIA*10:.1f}mm  '
            f'{total_fail_axis}/{total_pts_axis}  Pass Rate: {overall_pass:.2f}%'
        )
    elif ANALYSIS == 'dif':
        title_tag = f'Dose Difference {DD_CRITERIA:.0f}%'
        ax1.set_xlabel(
            f'Points outside {DD_CRITERIA:.1f}%  '
            f'{total_fail_axis}/{total_pts_axis}  Pass Rate: {overall_pass:.2f}%'
        )
    elif ANALYSIS == 'dist':
        title_tag = f'DTA {DTA_CRITERIA*10:.0f}mm'
        ax1.set_xlabel(
            f'Points outside {DTA_CRITERIA*10:.1f}mm  '
            f'{total_fail_axis}/{total_pts_axis}  Pass Rate: {overall_pass:.2f}%'
        )
    elif ANALYSIS == 'gam':
        title_tag = f'Gamma {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm'
        ax1.set_ylabel(f'$\\Gamma$ [{DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm]')
        ax1.set_xlabel(
            f'Points $\\Gamma$ > 1: {total_fail_axis}/{total_pts_axis}  '
            f'Pass Rate: {overall_pass:.2f}%'
        )
        if gvtot:
            bins    = 0.1
            gv_arr  = np.asarray(gvtot)
            weights = np.ones_like(gv_arr) / float(len(gv_arr))
            edges   = np.arange(0, max(gv_arr) + bins, bins)
            ax2.hist(gv_arr, bins=edges, weights=weights)
            ax2.set_xlabel(f'$\\Gamma$ [{DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm]')
            ax2.set_ylabel('Normalized Incidence')
    elif ANALYSIS == 'mppg':
        title_tag = (
            f'MPPG {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm/'
            f'{MPPG_DDTAIL:.0f}%Dmax'
        )
        ax3.set_xlabel(
            f'Outside {DD_CRITERIA:.1f}% in-field | {DTA_CRITERIA*10:.1f}mm pen | '
            f'{MPPG_DDTAIL:.1f}%Dmax tails  '
            f'{total_fail_axis}/{total_pts_axis}  Pass Rate: {overall_pass:.2f}%'
        )
    else:
        title_tag = 'Plots Only'

    safe_tag = title_tag.replace(' ', '_').replace('/', '-').replace('%', 'pct')
    if SAVE_FIGURES or SAVE_REPORT:
        axes_label = '+'.join(axes)
        fig.suptitle(
            f'{energy_label} — {SHEET1_NAME} vs {SHEET2_NAME}  {title_tag}  [{axes_label}]',
            fontsize=14
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.subplots_adjust(hspace=0.45)
    if SAVE_FIGURES:
        s1 = SHEET1_NAME.replace(' ', '')
        s2 = SHEET2_NAME.replace(' ', '')
        out_png = os.path.join(RESULTS_DIR, f"{s1}_{s2}_{energy_label}_{safe_tag}.png")
        fig.savefig(out_png, dpi=DPI, bbox_inches='tight')
        print(f"  Figure saved: {out_png}")

    if not SAVE_REPORT:
        plt.close(fig)
        fig = None

    return all_results, skipped, fig


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    from datetime import datetime
    os.makedirs(RESULTS_DIR, exist_ok=True)   # also creates COMPARISON_DIR and PROFILE_DIR
    os.makedirs(REPORTS_DIR, exist_ok=True)
    timestamp   = datetime.now().strftime("%Y%m%d_%H%M%S")
    s1 = SHEET1_NAME.replace(' ', ''); s2 = SHEET2_NAME.replace(' ', '')
    summary_xlsx = os.path.join(PROFILE_DIR, f"{s1}_{s2}_profile_summary_{timestamp}.xlsx")
    all_results = []

    # ── discover files ────────────────────────────────────────────────────────
    if FILE_MODE == 'flat':
        xlsx_files = sorted(glob.glob(os.path.join(BASE_PATH, FILE_FILTER)))
        file_pairs = []
        for p in xlsx_files:
            stem  = os.path.splitext(os.path.basename(p))[0]
            label = stem
            for suffix in ('_ProfileData', '_Profile', 'Data', '_PDDData', '_PDD'):
                if label.endswith(suffix):
                    label = label[:-len(suffix)]
                    break
            file_pairs.append((p, label))
    else:
        energy_dirs = sorted([
            d for d in os.listdir(BASE_PATH)
            if os.path.isdir(os.path.join(BASE_PATH, d))
            and os.path.join(BASE_PATH, d) != os.path.join(BASE_PATH, "Results")
        ])
        file_pairs = []
        for energy_label in energy_dirs:
            folder     = os.path.join(BASE_PATH, energy_label)
            prof_files = glob.glob(os.path.join(folder, '*Profile*.xlsx'))
            if not prof_files:
                print(f"No Profile xlsx found in {folder} — skipping.")
                continue
            file_pairs.append((prof_files[0], energy_label))

    if not file_pairs:
        print("No files found matching the filter. Check BASE_PATH and FILE_FILTER.")
        return

    # ── process each file ─────────────────────────────────────────────────────
    all_skipped = []
    all_figs = []
    for xlsx_path, energy_label in file_pairs:
        results, skipped, fig = run_one_file(xlsx_path, energy_label)
        all_results.extend(results)
        all_skipped.extend(skipped)
        if fig is not None:
            all_figs.append((energy_label, fig))

    # ── save summary CSV ──────────────────────────────────────────────────────
    if not all_results:
        print("\nNo results to save.")
        return

    # Normalise columns — fill missing metric columns with NaN
    for row in all_results:
        row.setdefault('MeanDoseDiff_pct', float('nan'))
        row.setdefault('MeanDTA_mm',       float('nan'))
        row.setdefault('MeanGamma',        float('nan'))

    df_out = pd.DataFrame(all_results)

    # ── helpers ───────────────────────────────────────────────────────────────
    def weighted_pass(grp):
        t = grp['TotalPoints'].sum()
        f = grp['FailPoints'].sum()
        return (t - f) / t * 100 if t > 0 else float('nan')

    def fmt(v):
        return f"{v:.2f}%" if (not isinstance(v, float) or not np.isnan(v)) else "—"

    _preferred   = ['6X', '6FFF', '8FFF', '10X', '10FFF', '15X']
    _actual      = df_out['Energy'].unique()
    _matched     = [e for e in _preferred if e in _actual]
    energy_order = _matched if _matched else sorted(_actual)

    # ── criteria label (used in console + Excel titles) ───────────────────────
    _analysis_labels = {
        'comp': 'Composite', 'dif': 'Dose Diff', 'dist': 'DTA',
        'gam':  'Gamma',     'mppg': 'MPPG',     'plot': 'Plot only',
    }
    _aname = _analysis_labels.get(ANALYSIS, ANALYSIS)
    if ANALYSIS in ('comp', 'dif', 'dist', 'gam'):
        _criteria_str = f"{_aname}  {DD_CRITERIA:.4g}%/{DTA_CRITERIA*10:.4g}mm"
    elif ANALYSIS == 'mppg':
        _criteria_str = f"MPPG  {DD_CRITERIA:.4g}%/{DTA_CRITERIA*10:.4g}mm/{MPPG_DDTAIL:.4g}%"
    else:
        _criteria_str = _aname

    # ── per-energy Depth × FS matrices ────────────────────────────────────────
    all_matrix_blocks = []   # list of (energy, df_matrix) for CSV output
    energy_texts = {}        # energy -> captured summary text for PDF

    for e in energy_order:
        de = df_out[df_out['Energy'] == e]
        depths = sorted(de['Depth_cm'].unique())
        fss    = sorted(de['FS'].unique())
        fs_col_names = [
            f"{int(fs)} cm" if float(fs).is_integer() else f"{fs} cm"
            for fs in fss
        ]

        rows = []
        for d in depths:
            row = {'Depth': f"{d:g} cm"}
            for fs, col in zip(fss, fs_col_names):
                sub = de[(de['Depth_cm'] == d) & (de['FS'] == fs)]
                row[col] = fmt(weighted_pass(sub)) if len(sub) > 0 else '—'
            row['All FS'] = fmt(weighted_pass(de[de['Depth_cm'] == d]))
            rows.append(row)

        # Footer: per-FS totals
        footer = {'Depth': 'All Depths'}
        for fs, col in zip(fss, fs_col_names):
            sub = de[de['FS'] == fs]
            footer[col] = fmt(weighted_pass(sub)) if len(sub) > 0 else '—'
        footer['All FS'] = fmt(weighted_pass(de))
        rows.append(footer)

        df_mat = pd.DataFrame(rows)
        all_matrix_blocks.append((e, df_mat))

        _buf = io.StringIO()
        with contextlib.redirect_stdout(_buf):
            print(f"\n{'='*60}")
            print(f"  {e}  —  Depth x FS pass rate  ({SHEET1_NAME} vs {SHEET2_NAME})  [{_criteria_str}]")
            print(f"{'='*60}")
            print(df_mat.to_string(index=False))
        energy_texts[e] = _buf.getvalue()
        print(energy_texts[e], end='')

    # ── skipped scan summary ──────────────────────────────────────────────────
    print(f"\n{'='*60}")
    if all_skipped:
        print(f"  SKIPPED SCANS ({len(all_skipped)} total):")
        for msg in all_skipped:
            print(f"    - {msg}")
    else:
        print("  No scans skipped.")
    print(f"{'='*60}")

    # ── write Excel workbook (Summary tab + Detail tab) ──────────────────────
    def _write_xlsx(path):
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            # ── Sheet 1: Summary (Depth × FS matrices, one block per energy) ──
            row_cursor = 0
            summary_entries = []   # (start_row, df_block, write_header)
            for e, df_mat in all_matrix_blocks:
                header_df = pd.DataFrame([[f'{e}  ({SHEET1_NAME} vs {SHEET2_NAME})  [{_criteria_str}]']])
                summary_entries.append((row_cursor, header_df, False))  # no col header row
                row_cursor += 1
                summary_entries.append((row_cursor, df_mat, True))
                row_cursor += len(df_mat) + 2   # +1 for df_mat header row + 1 blank

            for start_row, df_block, with_header in summary_entries:
                df_block.to_excel(writer, sheet_name='Summary',
                                  index=False, header=with_header, startrow=start_row)

            # skipped scans appended below matrices
            if all_skipped:
                skip_df = pd.DataFrame({'Skipped scans': all_skipped})
                skip_df.to_excel(writer, sheet_name='Summary',
                                 index=False, startrow=row_cursor + 1)

            # ── Sheet 2: Detail rows ───────────────────────────────────────────
            df_out.to_excel(writer, sheet_name='Detail', index=False)

    try:
        _write_xlsx(summary_xlsx)
        print(f"\nSummary saved: {summary_xlsx}")
    except PermissionError:
        fallback = summary_xlsx.replace('.xlsx', f'_retry_{timestamp}.xlsx')
        _write_xlsx(fallback)
        print(f"\nFile locked (close it in Excel) — saved fallback: {fallback}")

    # ── PDF report ────────────────────────────────────────────────────────────
    if SAVE_REPORT and all_figs:
        from reportlab.lib.pagesizes import letter as rl_letter
        from reportlab.pdfgen import canvas as rl_canvas
        from reportlab.lib.utils import ImageReader

        pdf_path = os.path.join(PROFILE_DIR, os.path.basename(summary_xlsx).replace('.xlsx', '.pdf'))
        rl_w, rl_h = rl_letter
        c = rl_canvas.Canvas(pdf_path, pagesize=rl_letter)

        # ── summary front page(s) ─────────────────────────────────────────────
        c.setFont("Helvetica-Bold", 14)
        c.drawString(36, rl_h - 48, f"Profile Summary Report — {SHEET1_NAME} vs {SHEET2_NAME}")
        c.setFont("Helvetica", 11)
        c.drawString(36, rl_h - 68, f"Analysis: {_criteria_str}")
        c.setFont("Courier", 7)
        y = rl_h - 100
        for e, df_mat in all_matrix_blocks:
            header = f"  {e}  —  Depth x FS pass rate"
            if y < 80:
                c.showPage()
                c.setFont("Courier", 7)
                y = rl_h - 36
            c.drawString(36, y, header)
            y -= 10
            for line in df_mat.to_string(index=False).split('\n'):
                if y < 80:
                    c.showPage()
                    c.setFont("Courier", 7)
                    y = rl_h - 36
                c.drawString(36, y, line)
                y -= 10
            y -= 8  # gap between energy blocks
        c.showPage()

        for energy_label, fig in all_figs:
            etext = energy_texts.get(energy_label, f'  {energy_label}\n')
            # — text block (top of page) —
            c.setFont("Courier", 9)
            y = rl_h - 36
            for line in etext.split('\n'):
                c.drawString(36, y, line)
                y -= 12
                if y < rl_h * 0.45:
                    break
            # — figure (below text) —
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=REPORT_DPI, bbox_inches='tight')
            buf.seek(0)
            img = ImageReader(buf)
            iw, ih = img.getSize()
            fig_w = rl_w - 72
            fig_h = fig_w * (ih / float(iw))
            fig_top = y - 8
            if fig_top - fig_h < 18:
                fig_h = max(fig_top - 18, 10)
                fig_w = fig_h * (iw / float(ih))
            c.drawImage(img, 36, fig_top - fig_h, width=fig_w, height=fig_h)
            c.showPage()

            # — individual per-energy PDF in Reports/ —
            energy_pdf = os.path.join(REPORTS_DIR, f"{s1}_{s2}_{energy_label}.pdf")
            ec = rl_canvas.Canvas(energy_pdf, pagesize=rl_letter)
            ec.setFont("Courier", 9)
            ey = rl_h - 36
            for line in etext.split('\n'):
                ec.drawString(36, ey, line)
                ey -= 12
                if ey < rl_h * 0.45:
                    break
            buf2 = io.BytesIO()
            fig.savefig(buf2, format='png', dpi=REPORT_DPI, bbox_inches='tight')
            buf2.seek(0)
            img2 = ImageReader(buf2)
            iw2, ih2 = img2.getSize()
            fw2 = rl_w - 72
            fh2 = fw2 * (ih2 / float(iw2))
            ft2 = ey - 8
            if ft2 - fh2 < 18:
                fh2 = max(ft2 - 18, 10)
                fw2 = fh2 * (iw2 / float(ih2))
            ec.drawImage(img2, 36, ft2 - fh2, width=fw2, height=fh2)
            ec.save()

            plt.close(fig)

        c.save()
        print(f"PDF report saved: {pdf_path}")
        if platform.system() == 'Windows':
            os.startfile(pdf_path)
        elif platform.system() == 'Darwin':
            subprocess.call(['open', pdf_path])
        else:
            subprocess.call(['xdg-open', pdf_path])


if __name__ == '__main__':
    main()
