"""
Batch PDD Composite Analysis
-----------------------------
Loops over all energy folders, runs composite analysis on
Measurement vs MC sheets, saves a figure per energy and a
summary CSV of pass-rate statistics.

Parameters are set in the CONFIG section below.
"""

import os
import io
import glob
import contextlib
import platform
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend (no GUI needed)
from matplotlib import pyplot as plt
from scipy import interpolate as interp
from scipy.ndimage import gaussian_filter1d

from comp import dosedif, dta as dtafunc
from gamma import gamma as gammafunc

# ─────────────────────────────────────────────
#  CONFIG  — edit these as needed
# ─────────────────────────────────────────────
FILE_MODE = 'flat'
# 'hierarchical' : BASE_PATH contains energy subfolders (e.g. 6X/, 10FFF/),
#                  each holding one *PDD*.xlsx file.
# 'flat'         : BASE_PATH is a single folder of xlsx files; energy label is
#                  parsed from each filename (e.g. 6X_100SSD_PDDData.xlsx -> 6X_100SSD).
#                  Use FILE_FILTER to restrict which files are processed.

BASE_PATH = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\Combined Consortium Data\Processed Combined Data"


FILE_FILTER = '*PDD*.xlsx'    # glob used in 'flat' mode only (e.g. '*PDD*.xlsx', '*.xlsx')
FS_FILTER   = [1.5, 3.0, 4.5, 6.0, 9.0, 10.5, 12.0, 15.0, 20.0, 30.0, 40.0]      # list of field sizes [cm] to include, e.g. [10.0, 20.0, 30.0]; empty = all

SHEET1_NAME = "TPS SN 21"   # reference / measured sheet name
SHEET2_NAME = "SN 19"            # comparison sheet name

ANALYSIS      = 'comp'        # 'comp', 'dif', 'dist', 'gam', or 'plot'
DD_CRITERIA   = 2           # dose difference threshold [%]
DTA_CRITERIA  = 0.2           # DTA threshold [cm]  (2 mm)
NORM          = 1             # 1 = normalize to dmax, 2 = normalize to fixed depth per energy

# Fixed normalization depth [cm] per energy folder name
NORM_DEPTH = {
    '6X':   1.3,
    '6FFF': 1.3,
    '10X':  2.2,
    '10FFF':2.2,
    '15X':  3.5,
}
DEPTH_SHIFT1  = 0.0           # depth shift for sheet 1 [cm]
DEPTH_SHIFT2  = -0.15           # depth shift for sheet 2 [cm]
CUTOFF_DEPTH  = 0.1           # discard points shallower than this [cm]; GUI uses 0.1 cm
CONV_FWHM_CM  = 0.48          # PTW 31021: 2 × 2.4 mm cavity radius = 4.8 mm = 0.48 cm
CONV_TARGET   = 'none'      # which curve to convolve: 'none', 'curve1', 'curve2', or 'both'
MARKER_SIZE   = 4             # plot marker size
SCALE_TARGET  = 0             # 0=none, 1=scale curve1, 2=scale curve2, 3=scale both
SCALE_FACTOR  = 1.00          # multiplier applied after normalization
DPI           = 600           # figure output resolution
SAVE_FIGURES  = True           # set False to skip individual PNG figure generation
SAVE_REPORT   = True           # generate multi-page PDF report (summary + per-energy figures)
REPORT_DPI    = 150            # DPI for PDF report pages (lower = smaller file, faster to open)

_s1 = SHEET1_NAME.replace(' ', '')
_s2 = SHEET2_NAME.replace(' ', '')

# Build a criteria tag from analysis type and thresholds so runs with different
# criteria land in separate folders automatically (e.g. comp_1pct_1mm).
_dd  = f"{DD_CRITERIA:.4g}".rstrip('0').rstrip('.')
_dta = f"{DTA_CRITERIA*10:.4g}".rstrip('0').rstrip('.')
_criteria_tag = f"{ANALYSIS}_{_dd}pct_{_dta}mm"

# Output folder hierarchy:
#   Results/
#     <SN1>_vs_<SN2>/
#       PDD/
#         <criteria>/
#           <summary>.csv
#           <summary>.pdf          ← all-energy report
#           Individual Figures/    ← per-energy PNGs
#           Individual Reports/    ← per-energy PDFs
COMPARISON_DIR = os.path.join(BASE_PATH, "Results", f"{_s1}_vs_{_s2}")
PDD_DIR        = os.path.join(COMPARISON_DIR, "PDD", _criteria_tag)
RESULTS_DIR    = os.path.join(PDD_DIR, "Individual Figures")
REPORTS_DIR    = os.path.join(PDD_DIR, "Individual Reports")
# ─────────────────────────────────────────────


# ── helpers copied from PDDCompare.py ────────

def apply_detector_convolution(profile_y, profile_dose, fwhm_cm):
    """Gaussian blur along depth axis to simulate detector volume averaging."""
    if fwhm_cm <= 0:
        return profile_dose
    sigma_cm  = fwhm_cm / 2.355
    step_cm   = np.mean(np.diff(np.asarray(profile_y, dtype=float)))
    sigma_idx = sigma_cm / step_cm
    dose_arr  = np.asarray(profile_dose, dtype=float)
    return gaussian_filter1d(dose_arr, sigma=sigma_idx)


def downsample_to_native(x_interp, vals_interp, x_native):
    """Thin interpolated output back to roughly native data spacing."""
    x_interp    = np.asarray(x_interp,    dtype=float)
    vals_interp = np.asarray(vals_interp, dtype=float)
    x_native    = np.asarray(x_native,    dtype=float)

    if x_interp.size < 2 or x_native.size < 2:
        return x_interp, vals_interp

    xu = np.unique(x_native[np.isfinite(x_native)])
    xu.sort()
    diffs = np.diff(xu)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return x_interp, vals_interp

    native_step = float(np.median(diffs))
    xi = np.unique(x_interp[np.isfinite(x_interp)])
    xi.sort()
    interp_step = float(np.median(np.diff(xi))) if xi.size >= 2 else 0.0
    if interp_step <= 0:
        return x_interp, vals_interp

    factor = max(1, int(round(native_step / interp_step)))
    return x_interp[::factor], vals_interp[::factor]


# ── analysis for one file ─────────────────────

def run_one_file(xlsx_path, energy_label):
    """
    Run composite analysis on one PDD Excel file.
    Returns a list of result dicts (one per FS).
    """
    print(f"\n{'='*60}")
    print(f"  {energy_label}  —  {os.path.basename(xlsx_path)}")
    print(f"{'='*60}")

    df1 = pd.read_excel(xlsx_path, sheet_name=SHEET1_NAME, header=0)
    df2 = pd.read_excel(xlsx_path, sheet_name=SHEET2_NAME, header=0)
    df1 = df1.sort_values(by=['FS', 'Axis', 'Pos'])
    df2 = df2.sort_values(by=['FS', 'Axis', 'Pos'])

    # FS values common to both sheets (Z-axis rows only)
    fs1 = set(df1.loc[df1['Axis'] == 'Z', 'FS'].unique())
    fs2 = set(df2.loc[df2['Axis'] == 'Z', 'FS'].unique())
    fsl = sorted(fs1 & fs2)
    if FS_FILTER:
        fsl = [fs for fs in fsl if fs in FS_FILTER]

    if not fsl:
        print("  No common field sizes found — skipping.")
        return []

    # ── figure setup ─────────────────────────
    if SAVE_FIGURES or SAVE_REPORT:
        if ANALYSIS == 'comp':
            fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 11),
                gridspec_kw={'height_ratios': [1.5, 1, 1]})
            plt.rcParams.update({'font.size': 20})
            ax1.set_ylabel('DTA [mm]')
            ax2.set_ylabel('Dose Difference [%]')
        elif ANALYSIS == 'gam':
            fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 11),
                gridspec_kw={'height_ratios': [1.5, 1, 1]})
            plt.rcParams.update({'font.size': 20})
            ax1.set_ylabel(f'$\\Gamma$ [{DD_CRITERIA:.4g}%/{DTA_CRITERIA*10:.4g}mm]')
            ax2.set_ylabel('Normalized Incidence')
        elif ANALYSIS in ('dif', 'dist'):
            fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15, 9),
                gridspec_kw={'height_ratios': [1.5, 1]})
            plt.rcParams.update({'font.size': 18})
            ax1.set_ylabel('Dose Difference [%]' if ANALYSIS == 'dif' else 'DTA [mm]')
        else:  # 'plot'
            fig, ax0 = plt.subplots(1, 1, figsize=(15, 8))
            plt.rcParams.update({'font.size': 16})
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
        ax0 = ax1 = ax2 = _NullAx()

    ax0.plot([], '+r', ms=10, label=SHEET1_NAME)
    ax0.plot([], '.k', ms=10, label=SHEET2_NAME)
    ax0.legend()
    ax0.set_ylabel('Percentage Depth Dose [%]')
    ax0.set_xlabel('Depth [cm]')

    results = []
    total_points_cumulative = 0
    total_fails_cumulative  = 0
    gvtot = []

    for fs_val in fsl:
        # ── extract data ───────────────────
        mask1 = (df1['FS'] == fs_val) & (df1['Axis'] == 'Z')
        mask2 = (df2['FS'] == fs_val) & (df2['Axis'] == 'Z')
        y1 = df1.loc[mask1, 'Pos'].reset_index(drop=True)
        d1 = df1.loc[mask1, 'Dose'].reset_index(drop=True)
        y2 = df2.loc[mask2, 'Pos'].reset_index(drop=True)
        d2 = df2.loc[mask2, 'Dose'].reset_index(drop=True)

        if len(y1) < 2 or len(y2) < 2:
            print(f"  FS {fs_val}: not enough points — skipping.")
            continue

        # ── check for duplicate depth positions ──────────────────────────
        dup1 = y1[y1.duplicated(keep=False)]
        dup2 = y2[y2.duplicated(keep=False)]
        if not dup1.empty:
            print(f"  FS {fs_val}: WARNING — {SHEET1_NAME} has duplicate Pos values "
                  f"{sorted(dup1.unique().tolist())} — skipping this FS.")
            continue
        if not dup2.empty:
            print(f"  FS {fs_val}: WARNING — {SHEET2_NAME} has duplicate Pos values "
                  f"{sorted(dup2.unique().tolist())} — skipping this FS.")
            continue
        # Ensure strictly increasing order
        y1 = y1.reset_index(drop=True); d1 = d1.reset_index(drop=True)
        y2 = y2.reset_index(drop=True); d2 = d2.reset_index(drop=True)
        sort1 = y1.argsort(); y1, d1 = y1.iloc[sort1].reset_index(drop=True), d1.iloc[sort1].reset_index(drop=True)
        sort2 = y2.argsort(); y2, d2 = y2.iloc[sort2].reset_index(drop=True), d2.iloc[sort2].reset_index(drop=True)

        # ── depth shifts ───────────────────
        y1 = y1 + DEPTH_SHIFT1
        y2 = y2 + DEPTH_SHIFT2

        # ── cutoff ─────────────────────────
        if CUTOFF_DEPTH > 0:
            m1 = y1 >= CUTOFF_DEPTH
            m2 = y2 >= CUTOFF_DEPTH
            y1, d1 = y1[m1].reset_index(drop=True), d1[m1].reset_index(drop=True)
            y2, d2 = y2[m2].reset_index(drop=True), d2[m2].reset_index(drop=True)
            if len(y1) < 2 or len(y2) < 2:
                print(f"  FS {fs_val}: cutoff removed too many points — skipping.")
                continue

        # ── normalization ──────────────────
        if NORM == 1:
            d1 = d1 / d1.max() * 100
            d2 = d2 / d2.max() * 100
        elif NORM == 2:
            fixed_depth = NORM_DEPTH.get(energy_label, None)
            if fixed_depth is None:
                print(f"  WARNING: no fixed depth for '{energy_label}', falling back to max norm.")
                d1 = d1 / d1.max() * 100
                d2 = d2 / d2.max() * 100
            else:
                d1_fixed = float(interp.pchip(y1, d1)(fixed_depth))
                d2_fixed = float(interp.pchip(y2, d2)(fixed_depth))
                d1 = d1 / d1_fixed * 100
                d2 = d2 / d2_fixed * 100

        # ── detector convolution ──────────────────
        if CONV_TARGET in ('curve1', 'both'):
            d1 = apply_detector_convolution(y1, d1, CONV_FWHM_CM)
        if CONV_TARGET in ('curve2', 'both'):
            d2 = apply_detector_convolution(y2, d2, CONV_FWHM_CM)

        # ── profile scaling (after normalization) ─────
        if SCALE_TARGET in (1, 3):
            d1 = d1 * SCALE_FACTOR
        if SCALE_TARGET in (2, 3):
            d2 = d2 * SCALE_FACTOR

        # ── plot curves ───────────────────
        ax0.plot(y1, d1, '+r', ms=MARKER_SIZE)
        ax0.plot(y2, d2, '.k', ms=MARKER_SIZE)

        # ── analysis metrics & subplot plots ──
        xlim = (max(float(y1.min()), float(y2.min())),
                min(float(y1.max()), float(y2.max())))
        ax0.set_xlim(*xlim)

        if ANALYSIS == 'plot':
            ax0.set_ylim(0, max(1.05 * np.max(d1), 1.05 * np.max(d2)))

        elif ANALYSIS == 'comp':
            dtax, dtav = dtafunc(y1, d1, y2, d2, DTA_CRITERIA)
            difx, difv = dosedif(y1, d1, y2, d2, 0)
            dtax, dtav = downsample_to_native(dtax, dtav, y1)
            difx, difv = downsample_to_native(difx, difv, y1)
            dtafail  = np.where(np.abs(dtav) > DTA_CRITERIA)[0]
            ddiffail = np.where(np.abs(difv) > DD_CRITERIA)[0]
            fail_set = set(dtafail) & set(ddiffail)
            totalfail = len(fail_set)
            total     = len(dtav)
            mean_dd  = float(np.mean(np.abs(difv))) if total > 0 else 0.0
            mean_dta = float(np.mean(np.abs(dtav))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDoseDiff_pct': round(mean_dd, 3),
                'MeanDTA_mm': round(mean_dta * 10, 3), 'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDD={mean_dd:.2f}%  MeanDTA={mean_dta*10:.2f}mm  Fail={totalfail}/{total}")
            for ax in (ax1, ax2):
                ax.set_xlim(*xlim)
            ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
            ax2.plot(difx, difv,      '.g', ms=0.5)
            for n in fail_set:
                ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)
                ax2.plot(difx[n], difv[n],      '.r', ms=0.5)

        elif ANALYSIS == 'dif':
            difx, difv = dosedif(y1, d1, y2, d2, 0)
            difx, difv = downsample_to_native(difx, difv, y1)
            ddiffail = np.where(np.abs(difv) > DD_CRITERIA)[0]
            totalfail = len(ddiffail)
            total     = len(difv)
            mean_dd  = float(np.mean(np.abs(difv))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDoseDiff_pct': round(mean_dd, 3),
                'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDD={mean_dd:.2f}%  Fail={totalfail}/{total}")
            ax1.set_xlim(*xlim)
            ax1.plot(difx, difv, '.g', ms=0.5)
            for n in ddiffail:
                ax1.plot(difx[n], difv[n], '.r', ms=0.5)

        elif ANALYSIS == 'dist':
            dtax, dtav = dtafunc(y1, d1, y2, d2, DTA_CRITERIA)
            dtax, dtav = downsample_to_native(dtax, dtav, y1)
            dtafail  = np.where(np.abs(dtav) > DTA_CRITERIA)[0]
            totalfail = len(dtafail)
            total     = len(dtav)
            mean_dta = float(np.mean(np.abs(dtav))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDTA_mm': round(mean_dta * 10, 3),
                'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDTA={mean_dta*10:.2f}mm  Fail={totalfail}/{total}")
            ax1.set_xlim(*xlim)
            ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
            for n in dtafail:
                ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)

        elif ANALYSIS == 'gam':
            gx, gv = gammafunc(y1, d1, y2, d2, DD_CRITERIA, DTA_CRITERIA, 0, 0, 0.01)
            gx, gv = downsample_to_native(gx, gv, y1)
            gv_a = np.asarray(gv); gx_a = np.asarray(gx)
            valid  = gv_a >= 0   # match GUI denominator (gv >= 0)
            total  = int(valid.sum())
            passed = int(np.count_nonzero(gv_a[valid] <= 1.0))
            totalfail = total - passed
            passrate  = passed / total * 100 if total > 0 else 0.0
            gvmean    = float(np.mean(gv_a[valid])) if total > 0 else float('nan')
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            gvtot.extend(gv_a[valid].tolist())
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanGamma': round(gvmean, 3),
                'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanGamma={gvmean:.3f}  Fail={totalfail}/{total}")
            ax1.set_xlim(*xlim)
            ax1.plot(gx_a[gv_a >  1], gv_a[gv_a >  1], '.r', ms=MARKER_SIZE)
            ax1.plot(gx_a[gv_a <= 1], gv_a[gv_a <= 1], '.g', ms=MARKER_SIZE)

    # ── figure summary label & title ──────────
    if total_points_cumulative > 0:
        overall_pass = 100 * (1 - total_fails_cumulative / total_points_cumulative)
    else:
        overall_pass = 0.0

    if ANALYSIS == 'comp':
        ax2.set_xlabel(f'Points outside {DD_CRITERIA:.1f}% & {DTA_CRITERIA*10:.1f}mm  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'Composite {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm'
    elif ANALYSIS == 'dif':
        ax1.set_xlabel(f'Points outside {DD_CRITERIA:.1f}%  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'Dose Difference {DD_CRITERIA:.0f}%'
    elif ANALYSIS == 'dist':
        ax1.set_xlabel(f'Points outside {DTA_CRITERIA*10:.1f}mm  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'DTA {DTA_CRITERIA*10:.0f}mm'
    elif ANALYSIS == 'gam':
        title_tag = f'Gamma {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm'
        ax1.set_xlabel(f'Points $\\Gamma$ > 1: {total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        if gvtot:
            gv_arr  = np.asarray(gvtot)
            bins    = 0.1
            weights = np.ones_like(gv_arr) / float(len(gv_arr))
            edges   = np.arange(0, max(gv_arr) + bins, bins)
            ax2.hist(gv_arr, bins=edges, weights=weights)
            ax2.set_xlabel(f'$\\Gamma$ [{DD_CRITERIA:.4g}%/{DTA_CRITERIA*10:.4g}mm]')
            ax2.set_ylabel('Normalized Incidence')
    else:
        title_tag = 'Plots Only'

    stem     = os.path.splitext(os.path.basename(xlsx_path))[0]
    safe_tag = title_tag.replace(' ', '_').replace('/', '-').replace('%', 'pct')
    if SAVE_FIGURES or SAVE_REPORT:
        fig.suptitle(f'{energy_label} — {SHEET1_NAME} vs {SHEET2_NAME}  {title_tag}', fontsize=14)
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

    return results, fig


# ── main ─────────────────────────────────────

def main():
    from datetime import datetime
    os.makedirs(RESULTS_DIR, exist_ok=True)   # also creates COMPARISON_DIR and PDD_DIR
    os.makedirs(REPORTS_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    s1 = SHEET1_NAME.replace(' ', ''); s2 = SHEET2_NAME.replace(' ', '')
    summary_csv = os.path.join(PDD_DIR, f"{s1}_{s2}_pdd_summary_{timestamp}.csv")
    all_results = []

    if FILE_MODE == 'flat':
        # All xlsx files sit directly in BASE_PATH; derive label from filename.
        xlsx_files = sorted(glob.glob(os.path.join(BASE_PATH, FILE_FILTER)))
        file_pairs = []
        for p in xlsx_files:
            stem = os.path.splitext(os.path.basename(p))[0]
            # Strip common trailing suffixes to get a tidy energy label
            label = stem
            for suffix in ('_PDDData', '_ProfileData', '_PDD', '_Profile', 'Data'):
                if label.endswith(suffix):
                    label = label[:-len(suffix)]
                    break
            file_pairs.append((p, label))
    else:
        # Hierarchical: one energy subfolder per energy, each with a *PDD*.xlsx
        energy_dirs = sorted([
            d for d in os.listdir(BASE_PATH)
            if os.path.isdir(os.path.join(BASE_PATH, d))
            and os.path.join(BASE_PATH, d) != os.path.join(BASE_PATH, "Results")
        ])
        file_pairs = []
        for energy_label in energy_dirs:
            folder = os.path.join(BASE_PATH, energy_label)
            pdd_files = glob.glob(os.path.join(folder, '*PDD*.xlsx'))
            if not pdd_files:
                print(f"No PDD xlsx found in {folder} — skipping.")
                continue
            file_pairs.append((pdd_files[0], energy_label))

    all_figs = []
    for xlsx_path, energy_label in file_pairs:
        _buf = io.StringIO()
        with contextlib.redirect_stdout(_buf):
            results, fig = run_one_file(xlsx_path, energy_label)
        energy_text = _buf.getvalue()
        print(energy_text, end='')
        all_results.extend(results)
        if fig is not None:
            all_figs.append((energy_label, energy_text, fig))

    # ── save summary CSV ─────────────────────
    if all_results:
        df_out = pd.DataFrame(all_results)
        for col in ('MeanDoseDiff_pct', 'MeanDTA_mm', 'MeanGamma'):
            if col not in df_out.columns:
                df_out[col] = float('nan')

        # ── criteria label ────────────────────────────────────────────
        _analysis_labels = {
            'comp': 'Composite', 'dif': 'Dose Diff',
            'dist': 'DTA',       'gam': 'Gamma',      'plot': 'Plot only',
        }
        _aname = _analysis_labels.get(ANALYSIS, ANALYSIS)
        if ANALYSIS in ('comp', 'dif', 'dist', 'gam'):
            _criteria_str = f"{_aname}  {DD_CRITERIA:.4g}%/{DTA_CRITERIA*10:.4g}mm"
        else:
            _criteria_str = _aname

        # ── build pivot table (pass rate % by energy × FS) ───────────
        _preferred = ['6X', '6FFF', '8FFF', '10X', '10FFF', '15X']
        _actual = df_out['Energy'].unique()
        _matched = [e for e in _preferred if e in _actual]
        energy_order = _matched if _matched else sorted(_actual)
        fss = sorted(df_out['FS'].unique())

        def weighted_pass(g):
            t = g['TotalPoints'].sum()
            f = g['FailPoints'].sum()
            return (t - f) / t * 100 if t > 0 else 0.0

        pivot = df_out.pivot_table(
            index='Energy', columns='FS',
            values='PassRate_pct', aggfunc='first'
        ).reindex(index=energy_order, columns=fss)

        all_fs   = df_out.groupby('Energy').apply(weighted_pass).reindex(energy_order)
        all_en   = df_out.groupby('FS').apply(weighted_pass).reindex(fss)
        all_data = weighted_pass(df_out)

        def fmt(v):
            return f"{v:.2f}%" if not np.isnan(v) else "—"

        fs_col_names = [f"{int(fs)} cm" if float(fs).is_integer() else f"{fs} cm"
                        for fs in fss]
        table_rows = []
        for e in energy_order:
            row = {'Energy': e}
            for fs, col in zip(fss, fs_col_names):
                row[col] = fmt(pivot.loc[e, fs])
            row['All Field Sizes'] = fmt(all_fs[e])
            table_rows.append(row)

        footer = {'Energy': 'All Energies'}
        for fs, col in zip(fss, fs_col_names):
            footer[col] = fmt(all_en[fs])
        footer['All Field Sizes'] = f"All Data: {fmt(all_data)}"
        table_rows.append(footer)

        df_table = pd.DataFrame(table_rows)

        try:
            df_table.to_csv(summary_csv, index=False)
            print(f"\nSummary CSV saved: {summary_csv}")
        except PermissionError:
            from datetime import datetime
            fallback = summary_csv.replace('.csv', f'_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
            df_table.to_csv(fallback, index=False)
            print(f"\nCSV locked (close it in Excel first next time).")
            print(f"Saved to fallback: {fallback}")
        print(f"\n{'='*60}")
        print(f"  PDD pass rate  ({SHEET1_NAME} vs {SHEET2_NAME})  [{_criteria_str}]")
        print(f"{'='*60}")
        print(df_table.to_string(index=False))

        # ── PDF report ────────────────────────────────────────────────────────
        if SAVE_REPORT and all_figs:
            from reportlab.lib.pagesizes import letter as rl_letter
            from reportlab.pdfgen import canvas as rl_canvas
            from reportlab.lib.utils import ImageReader

            pdf_path = os.path.join(PDD_DIR, os.path.basename(summary_csv).replace('.csv', '.pdf'))
            rl_w, rl_h = rl_letter
            c = rl_canvas.Canvas(pdf_path, pagesize=rl_letter)

            # ── summary front page ────────────────────────────────────────────
            c.setFont("Helvetica-Bold", 14)
            c.drawString(36, rl_h - 48, f"PDD Summary Report — {SHEET1_NAME} vs {SHEET2_NAME}")
            c.setFont("Helvetica", 11)
            c.drawString(36, rl_h - 68, f"Analysis: {_criteria_str}")
            c.setFont("Courier", 7)
            y = rl_h - 100
            for line in df_table.to_string(index=False).split('\n'):
                c.drawString(36, y, line)
                y -= 10
            c.showPage()

            for energy_label, energy_text, fig in all_figs:
                # — text block (top of page) —
                c.setFont("Courier", 9)
                y = rl_h - 36
                for line in energy_text.split('\n'):
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
                for line in energy_text.split('\n'):
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
    else:
        print("\nNo results to save.")


if __name__ == '__main__':
    main()
