"""
Batch PDD Composite Analysis
-----------------------------
Loops over all energy folders, runs composite analysis on
Measurement vs MC sheets, saves a figure per energy and a
summary CSV of pass-rate statistics.

Parameters are set in the CONFIG section below.
"""

import os
import sys
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

# ── EPOM shift helpers (inlined from PDDCompare.py to avoid launching GUI) ───
import enum
from typing import NamedTuple
from scipy.interpolate import UnivariateSpline, BSpline, PPoly

class _IonChamberData(NamedTuple):
    radius_cavity_mm: float
    wall_thickness_mm: float

class _IonChambers(str, enum.Enum):
    PTW_31021 = "PTW31021"
    IBA_CC13  = "IBA CC13"
    PTW_31010 = "PTW31010"
    def _parameters(self):
        if self is _IonChambers.PTW_31021: return _IonChamberData(2.4, 0.57 + 0.09)
        if self is _IonChambers.PTW_31010: return _IonChamberData(2.75, 0.55 + 0.15)
        if self is _IonChambers.IBA_CC13:  return _IonChamberData(3.0, 0.4)
        return _IonChamberData(2.4, 0.66)
    @property
    def outer_radius(self):
        p = self._parameters(); return p.radius_cavity_mm + p.wall_thickness_mm
    @property
    def simple_photon_effective_point_of_measurement(self):
        return 0.6 * self._parameters().radius_cavity_mm
    @property
    def simple_electron_effective_point_of_measurement(self):
        return 0.5 * self._parameters().radius_cavity_mm

class _Modality(str, enum.Enum):
    PHOTON   = "PHOTON"
    ELECTRON = "ELECTRON"

class _PDDCurve:
    def __init__(self, z, dose, unit="cm", lam=0.0):
        z = np.asarray(z, dtype=float); dose = np.asarray(dose, dtype=float)
        idx = np.argsort(z); self.z = z[idx]; self.dose = dose[idx]
        self.z_distance_unit = unit
        N = max(1, self.z.size); s = max(0.0, lam * N * float(np.var(self.dose)))
        k = min(3, max(1, N - 1))
        self.curve = UnivariateSpline(self.z, self.dose, s=s, k=k)
    def _ppoly(self):
        t, c, k = self.curve._eval_args; return PPoly.from_spline(BSpline(t, c, k))
    def _deriv(self, order=1): return self._ppoly().derivative(order)
    def _mm_to_z(self, mm): return mm / (10.0 if self.z_distance_unit == "cm" else 1.0)
    def _surface(self):
        try:
            roots = self._ppoly().derivative(2).roots()
            if roots.size:
                return float(roots[int(np.argmax(self._deriv(1)(roots)))])
        except Exception:
            pass
        zmin, zmax = float(self.z.min()), float(self.z.max())
        win = 2.0 if self.z_distance_unit == "cm" else 20.0
        zfine = np.linspace(zmin, min(zmax, zmin + win), max(1000, 10 * len(self.z)))
        d1 = np.nan_to_num(self._deriv(1)(zfine), nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
        return float(zfine[int(np.argmax(d1))])
    def epom_shift(self, chamber, modality):
        epom_mm = (chamber.simple_photon_effective_point_of_measurement
                   if modality == _Modality.PHOTON
                   else chamber.simple_electron_effective_point_of_measurement)
        return self._mm_to_z(chamber.outer_radius - epom_mm) - self._surface()

def _map_detector(name):
    """Map a raw detector string to an _IonChambers enum, or None for no-shift detectors.
    Canonical names (TN31021, TN31010, TN60019, CC13) are matched first.
    Fuzzy fallbacks handle any remaining legacy variants.
    Returns None for no-shift detectors (TPS, TN60019 microDiamond).
    """
    key = (name or "").strip().lower().replace("-","").replace("_","").replace(" ","")
    # No-shift detectors
    if key.startswith("tps"):                      return None  # TPS M5, TPS M6, etc.
    if key in ("tn60019",) or "microdiamond" in key or "60019" in key: return None
    # PTW 31010
    if key in ("tn31010",) or "31010" in key:     return _IonChambers.PTW_31010
    # IBA CC13
    if key in ("cc13",) or "cc13" in key:         return _IonChambers.IBA_CC13
    # PTW 31021 — canonical TN31021, serial TN31XXX, or legacy name variants
    if key in ("tn31021",) or key.startswith("tn31") or key.startswith("tn32") or "31021" in key or "semiflex" in key:
                                                   return _IonChambers.PTW_31021
    # Unrecognised — warn and default to PTW31021
    print(f"[_map_detector] unrecognised detector '{name}' -- defaulting to PTW31021")
    return _IonChambers.PTW_31021

def compute_epom_shifts(z1, d1, det1, z2, d2, det2, modality="PHOTON", unit="cm"):
    import numpy as _np
    z1a = _np.asarray(z1, dtype=float); z2a = _np.asarray(z2, dtype=float)
    ic1 = _map_detector(det1); ic2 = _map_detector(det2)
    mod = _Modality(modality)
    def _shift(za, d, ic):
        if ic is None:
            return 0.0
        if not _np.any(za < 0):
            return 0.0
        return float(_PDDCurve(za, d, unit=unit).epom_shift(ic, mod))
    s1 = _shift(z1a, d1, ic1)
    s2 = _shift(z2a, d2, ic2)
    return s1, s2
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────
#  CONFIG  — edit these as needed
# ─────────────────────────────────────────────
FILE_MODE = 'flat'
# 'hierarchical' : BASE_PATH contains energy subfolders (e.g. 6X/, 10FFF/),
#                  each holding one *PDD*.xlsx file.
# 'flat'         : BASE_PATH is a single folder of xlsx files; energy label is
#                  parsed from each filename (e.g. 6X_100SSD_PDDData.xlsx -> 6X_100SSD).
#                  Use FILE_FILTER to restrict which files are processed.

BASE_PATH = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\Combined Consortium Data\100 SSD"

# Where outputs (figures, reports, summary csv) get written. Set to None to use
# BASE_PATH (legacy behavior); set to a short path like r"C:\OFresults" to keep
# total file paths under the Windows 260-character limit.
OUTPUT_BASE = None

FILE_FILTER = '*PDD*.xlsx'    # glob used in 'flat' mode only (e.g. '*PDD*.xlsx', '*.xlsx')
FS_FILTER   = []      # list of field sizes [cm] to include, e.g. [10.0, 20.0, 30.0]; empty = all

SHEET1_NAME = "SN 21"   # reference / measured sheet name; set 'none' to skip
SHEET2_NAME = "none"            # comparison sheet name; set 'none' to skip

ANALYSIS      = 'none'        # 'comp', 'dif', 'dist', 'gam', 'plot', or 'none' (plot only, no comparison)
DD_CRITERIA   = 2           # dose difference threshold [%]
DTA_CRITERIA  = 0.2           # DTA threshold [cm]  (2 mm)
NORM          = 1             # 1 = each dmax, 2 = fixed depth per energy, 3 = dmax of curve 1, 4 = dmax of curve 2

# Fixed normalization depth [cm] per energy folder name
NORM_DEPTH = {
    '6X':   1.3,
    '6FFF': 1.3,
    '10X':  2.2,
    '10FFF':2.2,
    '15X':  3.5,
}
DEPTH_SHIFT1  = 0.0           # depth shift for sheet 1 [cm]  (ignored when AUTO_EPOM_SHIFT=True)
DEPTH_SHIFT2  = 0.0          # depth shift for sheet 2 [cm]  (ignored when AUTO_EPOM_SHIFT=True)
AUTO_EPOM_SHIFT = False       # True = compute EPOM depth shift automatically per file (like GUI)
EPOM_DETECTOR1  = 'PTW31021' # detector for sheet 1: 'PTW31021', 'PTW31010', 'IBACC13'
EPOM_DETECTOR2  = 'PTW31021' # detector for sheet 2
EPOM_MODALITY   = 'PHOTON'   # 'PHOTON' or 'ELECTRON'
CUTOFF_DEPTH  = 0.0           # discard points shallower than this [cm]; GUI uses 0.1 cm
CONV_FWHM_CM  = 0.48          # PTW 31021: 2 × 2.4 mm cavity radius = 4.8 mm = 0.48 cm
CONV_TARGET   = 'none'      # which curve to convolve: 'none', 'curve1', 'curve2', or 'both'
MARKER_SIZE   = 4             # plot marker size
SCALE_TARGET  = 0             # 0=none, 1=scale curve1, 2=scale curve2, 3=scale both
SCALE_FACTOR  = 1.00          # multiplier applied after normalization
DPI           = 600           # figure output resolution
SAVE_FIGURES  = True           # set False to skip individual PNG figure generation
SAVE_REPORT   = True           # generate multi-page PDF report (summary + per-energy figures)
REPORT_DPI    = 150            # DPI for PDF report pages (lower = smaller file, faster to open)

def _sheet_active(name):
    """A sheet is skipped if its name is blank or none/na."""
    return str(name).strip().lower() not in ('', 'none', 'na', 'n/a')

USE1 = _sheet_active(SHEET1_NAME)
USE2 = _sheet_active(SHEET2_NAME)
PLOT_ONLY = (ANALYSIS == 'none') or (not USE1) or (not USE2)

if ANALYSIS in ('comp', 'dif', 'dist', 'gam') and not (USE1 and USE2):
    raise SystemExit(f"ANALYSIS='{ANALYSIS}' needs two sheets, but one of "
                     f"SHEET1_NAME/SHEET2_NAME is set to none.")

_s1 = SHEET1_NAME.replace(' ', '') if USE1 else 'none'
_s2 = SHEET2_NAME.replace(' ', '') if USE2 else 'none'

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
_OUT_ROOT      = OUTPUT_BASE if OUTPUT_BASE else BASE_PATH
COMPARISON_DIR = os.path.join(_OUT_ROOT, "Results", f"{_s1}_vs_{_s2}")
PDD_DIR        = os.path.join(COMPARISON_DIR, "PDD", _criteria_tag)
RESULTS_DIR    = os.path.join(PDD_DIR, "Individual Figures")
REPORTS_DIR    = os.path.join(PDD_DIR, "Individual Reports")
# ─────────────────────────────────────────────


# ── PDD lookup table (SSD=100, 10.5 cm field, normalized at dmax) ────────────
# First key of each dict is the dmax depth [cm] — used as buildup/past-dmax region boundary.
PDD_TABLE = {
    "6X":   {1.42: 1.000,  5.0: 0.865, 10.0: 0.668, 20.0: 0.382, 30.0: 0.218},
    "6FFF": {1.32: 1.000,  5.0: 0.848, 10.0: 0.636, 20.0: 0.346, 30.0: 0.190},
    "8FFF": {1.88: 1.000,  5.0: 0.891, 10.0: 0.693, 20.0: 0.406, 30.0: 0.240},
    "10X":  {2.22: 1.000,  5.0: 0.925, 10.0: 0.744, 20.0: 0.469, 30.0: 0.296},
    "10FFF":{2.16: 1.000,  5.0: 0.904, 10.0: 0.708, 20.0: 0.423, 30.0: 0.256},
    "15X":  {2.72: 1.000,  5.0: 0.949, 10.0: 0.774, 20.0: 0.502, 30.0: 0.325},
}


def dmax_for_energy(energy):
    """Return tabled dmax depth [cm] for an energy token, or None if missing."""
    if energy is None:
        return None
    key = str(energy).strip().upper().replace(' ', '').replace('MV', '')
    if key in PDD_TABLE:
        return float(next(iter(PDD_TABLE[key].keys())))
    return None


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

    df1 = None
    df2 = None
    if USE1:
        df1 = pd.read_excel(xlsx_path, sheet_name=SHEET1_NAME, header=0).sort_values(by=['FS', 'Axis', 'Pos'])
    if USE2:
        df2 = pd.read_excel(xlsx_path, sheet_name=SHEET2_NAME, header=0).sort_values(by=['FS', 'Axis', 'Pos'])

    # Energy from whichever sheet is present (for region boundary lookup); falls back to label.
    _energy_token = None
    _ref_df = df1 if df1 is not None else df2
    if _ref_df is not None and 'Energy' in _ref_df.columns:
        _vals = _ref_df['Energy'].dropna().unique()
        if len(_vals) > 0:
            _energy_token = str(_vals[0])
    if _energy_token is None:
        _energy_token = energy_label
    _region_dmax = dmax_for_energy(_energy_token)

    # Detector is resolved per-FS inside the loop (each scan may use a different detector)

    # Field sizes: intersection when comparing two sheets, union of active sheets in plot mode.
    fs1 = set(df1.loc[df1['Axis'] == 'Z', 'FS'].unique()) if df1 is not None else set()
    fs2 = set(df2.loc[df2['Axis'] == 'Z', 'FS'].unique()) if df2 is not None else set()
    if PLOT_ONLY:
        fsl = sorted(fs1 | fs2)
    else:
        fsl = sorted(fs1 & fs2)
    if FS_FILTER:
        fsl = [fs for fs in fsl if fs in FS_FILTER]

    if not fsl:
        print("  No field sizes found — skipping.")
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
        else:  # 'plot' or 'none'
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

    if USE1:
        ax0.plot([], '+r', ms=10, label=SHEET1_NAME)
    if USE2:
        ax0.plot([], '.k', ms=10, label=SHEET2_NAME)
    ax0.legend()
    ax0.set_ylabel('Percentage Depth Dose [%]')
    ax0.set_xlabel('Depth [cm]')

    results = []
    total_points_cumulative = 0
    total_fails_cumulative  = 0
    gvtot = []

    def _norm_at_depth(y, d, depth):
        return float(interp.pchip(y, d)(depth))

    def _extract_curve(df, sheet_name):
        """Pull (y, d, detector) for this FS's Z-axis scan; (None,None,None) if unusable."""
        mask = (df['FS'] == fs_val) & (df['Axis'] == 'Z')
        y = df.loc[mask, 'Pos'].reset_index(drop=True)
        d = df.loc[mask, 'Dose'].reset_index(drop=True)
        if len(y) < 2:
            return None, None, None
        dup = y[y.duplicated(keep=False)]
        if not dup.empty:
            print(f"  FS {fs_val}: WARNING — {sheet_name} has duplicate Pos values "
                  f"{sorted(dup.unique().tolist())} — skipping.")
            return None, None, None
        srt = y.argsort()
        y, d = y.iloc[srt].reset_index(drop=True), d.iloc[srt].reset_index(drop=True)
        det = None
        if 'Detector' in df.columns:
            vals = df.loc[mask, 'Detector'].dropna().unique()
            if len(vals) > 0:
                det = str(vals[0]).strip()
        return y, d, det

    for fs_val in fsl:
        # ── extract each active curve ───────────────────
        y1 = d1 = y2 = d2 = None
        det1 = det2 = None
        if df1 is not None and fs_val in fs1:
            y1, d1, det1 = _extract_curve(df1, SHEET1_NAME)
        if df2 is not None and fs_val in fs2:
            y2, d2, det2 = _extract_curve(df2, SHEET2_NAME)

        # Comparison modes require both curves; plot modes accept either.
        if not PLOT_ONLY and (y1 is None or y2 is None):
            print(f"  FS {fs_val}: missing a curve — skipping.")
            continue
        if y1 is None and y2 is None:
            print(f"  FS {fs_val}: no usable curve — skipping.")
            continue

        # ── depth shifts ───────────────────
        if AUTO_EPOM_SHIFT:
            _d1 = det1 or EPOM_DETECTOR1
            _d2 = det2 or EPOM_DETECTOR2
            s1 = s2 = 0.0
            if y1 is not None and y2 is not None:
                s1, s2 = compute_epom_shifts(y1, d1, _d1, y2, d2, _d2,
                                             modality=EPOM_MODALITY, unit="cm")
            elif y1 is not None:
                s1, _ = compute_epom_shifts(y1, d1, _d1, y1, d1, _d1,
                                            modality=EPOM_MODALITY, unit="cm")
            elif y2 is not None:
                s2, _ = compute_epom_shifts(y2, d2, _d2, y2, d2, _d2,
                                            modality=EPOM_MODALITY, unit="cm")
        else:
            s1, s2 = DEPTH_SHIFT1, DEPTH_SHIFT2
        if y1 is not None: y1 = y1 + s1
        if y2 is not None: y2 = y2 + s2

        # ── cutoff: always applied at CUTOFF_DEPTH ──
        if y1 is not None:
            m1 = y1 >= CUTOFF_DEPTH
            y1, d1 = y1[m1].reset_index(drop=True), d1[m1].reset_index(drop=True)
            if len(y1) < 2:
                y1 = d1 = None
        if y2 is not None:
            m2 = y2 >= CUTOFF_DEPTH
            y2, d2 = y2[m2].reset_index(drop=True), d2[m2].reset_index(drop=True)
            if len(y2) < 2:
                y2 = d2 = None
        if not PLOT_ONLY and (y1 is None or y2 is None):
            print(f"  FS {fs_val}: cutoff removed too many points — skipping.")
            continue
        if y1 is None and y2 is None:
            print(f"  FS {fs_val}: cutoff removed too many points — skipping.")
            continue

        # ── normalization ──────────────────
        # NORM 3/4 reference the other curve's dmax; fall back to each-own when absent.
        ref_depth = None
        eff_norm = NORM
        if NORM == 2:
            ref_depth = NORM_DEPTH.get(energy_label, None)
            if ref_depth is None:
                print(f"  WARNING: no fixed depth for '{energy_label}', using each-dmax.")
                eff_norm = 1
        elif NORM == 3:
            if y1 is not None:
                ref_depth = float(y1.iloc[d1.values.argmax()])
            else:
                eff_norm = 1
        elif NORM == 4:
            if y2 is not None:
                ref_depth = float(y2.iloc[d2.values.argmax()])
            else:
                eff_norm = 1

        def _normalize(y, d):
            if y is None:
                return None
            if eff_norm == 1:
                return d / d.max() * 100
            return d / _norm_at_depth(y, d, ref_depth) * 100

        d1 = _normalize(y1, d1)
        d2 = _normalize(y2, d2)

        # ── detector convolution ──────────────────
        if CONV_TARGET in ('curve1', 'both') and y1 is not None:
            d1 = apply_detector_convolution(y1, d1, CONV_FWHM_CM)
        if CONV_TARGET in ('curve2', 'both') and y2 is not None:
            d2 = apply_detector_convolution(y2, d2, CONV_FWHM_CM)

        # ── profile scaling (after normalization) ─────
        if SCALE_TARGET in (1, 3) and y1 is not None:
            d1 = d1 * SCALE_FACTOR
        if SCALE_TARGET in (2, 3) and y2 is not None:
            d2 = d2 * SCALE_FACTOR

        # ── plot curves ───────────────────
        if y1 is not None:
            ax0.plot(y1, d1, '+r', ms=MARKER_SIZE)
        if y2 is not None:
            ax0.plot(y2, d2, '.k', ms=MARKER_SIZE)

        # ── x-limits from available curve(s) ──
        if y1 is not None and y2 is not None:
            xlim = (max(float(y1.min()), float(y2.min())),
                    min(float(y1.max()), float(y2.max())))
        elif y1 is not None:
            xlim = (float(y1.min()), float(y1.max()))
        else:
            xlim = (float(y2.min()), float(y2.max()))
        ax0.set_xlim(*xlim)

        if PLOT_ONLY:
            _maxes = [1.05 * float(np.max(d)) for d in (d1, d2) if d is not None]
            if _maxes:
                ax0.set_ylim(0, max(_maxes))

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
            # Region breakdown (buildup vs past-dmax) — only how existing comp failures distribute.
            b_fails = b_total = p_fails = p_total = 0
            if _region_dmax is not None:
                _xarr = np.asarray(dtax, dtype=float)
                buildup_mask = _xarr <= _region_dmax
                past_mask    = _xarr >  _region_dmax
                fail_mask_arr = np.zeros(len(_xarr), dtype=bool)
                if fail_set:
                    _idx = np.fromiter(fail_set, dtype=int)
                    _idx = _idx[(_idx >= 0) & (_idx < len(_xarr))]
                    fail_mask_arr[_idx] = True
                b_total = int(buildup_mask.sum())
                p_total = int(past_mask.sum())
                b_fails = int((buildup_mask & fail_mask_arr).sum())
                p_fails = int((past_mask    & fail_mask_arr).sum())

            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDoseDiff_pct': round(mean_dd, 3),
                'MeanDTA_mm': round(mean_dta * 10, 3), 'TotalPoints': total, 'FailPoints': totalfail,
                'BuildupFails': b_fails, 'BuildupTotal': b_total,
                'PastFails': p_fails, 'PastTotal': p_total})
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
    # Title reflects which sheets are present.
    if USE1 and USE2:
        _sheet_label = f'{SHEET1_NAME} vs {SHEET2_NAME}'
    elif USE1:
        _sheet_label = SHEET1_NAME
    else:
        _sheet_label = SHEET2_NAME
    if SAVE_FIGURES or SAVE_REPORT:
        fig.suptitle(f'{energy_label} — {_sheet_label}  {title_tag}', fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.subplots_adjust(hspace=0.45)
    if SAVE_FIGURES:
        out_png = os.path.join(RESULTS_DIR, f"{_s1}_{_s2}_{energy_label}_{safe_tag}.png")
        fig.savefig(out_png, dpi=DPI, bbox_inches='tight')
        # Bypass stdout capture so the path line doesn't land in the per-energy PDF page.
        sys.__stdout__.write(f"  Figure saved: {out_png}\n")

    # ── Per-energy region breakdown (comp mode only) ──
    if ANALYSIS == 'comp' and _region_dmax is not None and results:
        tot_b_f = tot_b_t = tot_p_f = tot_p_t = 0
        print(f"\n  Region breakdown ({_energy_token}, dmax boundary = {_region_dmax:.2f} cm)")
        print(f"  {'FS [cm]':<10}{'Buildup Pass':>16}{'Buildup Fails/Total':>24}{'Past Pass':>14}{'Past Fails/Total':>24}")
        for r in results:
            b_f = r.get('BuildupFails', 0); b_t = r.get('BuildupTotal', 0)
            p_f = r.get('PastFails',    0); p_t = r.get('PastTotal',    0)
            b_pr = f"{(b_t - b_f)/b_t*100:.2f}%" if b_t > 0 else "—"
            p_pr = f"{(p_t - p_f)/p_t*100:.2f}%" if p_t > 0 else "—"
            print(f"  {r['FS']:<10.1f}{b_pr:>16}{f'{b_f}/{b_t}':>24}{p_pr:>14}{f'{p_f}/{p_t}':>24}")
            tot_b_f += b_f; tot_b_t += b_t; tot_p_f += p_f; tot_p_t += p_t
        b_pr_all = f"{(tot_b_t - tot_b_f)/tot_b_t*100:.2f}%" if tot_b_t > 0 else "—"
        p_pr_all = f"{(tot_p_t - tot_p_f)/tot_p_t*100:.2f}%" if tot_p_t > 0 else "—"
        print(f"  {'Overall':<10}{b_pr_all:>16}{f'{tot_b_f}/{tot_b_t}':>24}{p_pr_all:>14}{f'{tot_p_f}/{tot_p_t}':>24}")

    if not SAVE_REPORT:
        plt.close(fig)
        fig = None

    return results, fig


# ── main ─────────────────────────────────────

def _warn_path_length():
    """On Windows, warn if the expected output path may exceed the 260-char MAX_PATH limit."""
    if os.name != 'nt':
        return
    s1 = SHEET1_NAME.replace(' ', ''); s2 = SHEET2_NAME.replace(' ', '')
    safe_tag = _criteria_tag.replace('_', '-').upper()
    sample_name = f"{s1}_{s2}_10FFF_90SSD_{safe_tag}.png"
    sample_path = os.path.join(RESULTS_DIR, sample_name)
    n = len(sample_path)
    if n > 240:
        bar = '!' * 78
        print(f"\n{bar}")
        print(f"  WARNING: predicted output path is {n} chars (Windows limit is 260).")
        if n > 260:
            print(f"  This run WILL fail when saving figures.")
        else:
            print(f"  This is close to the limit — long energy labels may push it over.")
        print(f"  Fixes:")
        print(f"    • Shorten BASE_PATH (or move data out of OneDrive — that prefix alone is ~60 chars)")
        print(f"    • Set OUTPUT_BASE to a short path like  r\"C:\\BeamResults\"")
        print(f"    • Enable Windows long paths via registry (admin + reboot required)")
        print(f"  Sample path ({n} chars):")
        print(f"    {sample_path}")
        print(f"{bar}\n")


def main():
    from datetime import datetime
    print(f"  Input  : {BASE_PATH}")
    print(f"  Output : {_OUT_ROOT}")
    _warn_path_length()
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

        # ── Region breakdown tables (buildup / past-dmax) ─────────────────
        # Only emitted for comp analysis where per-region fail counts were recorded.
        df_buildup = None
        df_past    = None
        if ANALYSIS == 'comp' and all(c in df_out.columns for c in
                                      ('BuildupFails', 'BuildupTotal', 'PastFails', 'PastTotal')):
            def _region_pass(g, fail_col, tot_col):
                t = g[tot_col].sum()
                f = g[fail_col].sum()
                return (t - f) / t * 100 if t > 0 else float('nan')

            def _build_region_table(fail_col, tot_col):
                pivot_r = df_out.pivot_table(
                    index='Energy', columns='FS',
                    values=[fail_col, tot_col], aggfunc='sum'
                ).reindex(index=energy_order, columns=pd.MultiIndex.from_product(
                    [[fail_col, tot_col], fss]))
                rows = []
                for e in energy_order:
                    row = {'Energy': e}
                    for fs, col in zip(fss, fs_col_names):
                        try:
                            f_v = pivot_r.loc[e, (fail_col, fs)]
                            t_v = pivot_r.loc[e, (tot_col,  fs)]
                            row[col] = fmt((t_v - f_v) / t_v * 100) if t_v > 0 else "—"
                        except KeyError:
                            row[col] = "—"
                    e_rows = df_out[df_out['Energy'] == e]
                    row['All Field Sizes'] = fmt(_region_pass(e_rows, fail_col, tot_col))
                    rows.append(row)
                footer_r = {'Energy': 'All Energies'}
                for fs, col in zip(fss, fs_col_names):
                    fs_rows = df_out[df_out['FS'] == fs]
                    footer_r[col] = fmt(_region_pass(fs_rows, fail_col, tot_col))
                footer_r['All Field Sizes'] = f"All Data: {fmt(_region_pass(df_out, fail_col, tot_col))}"
                rows.append(footer_r)
                return pd.DataFrame(rows)

            df_buildup = _build_region_table('BuildupFails', 'BuildupTotal')
            df_past    = _build_region_table('PastFails',    'PastTotal')

        def _write_csv(path):
            with open(path, 'w', newline='', encoding='utf-8') as fh:
                df_table.to_csv(fh, index=False)
                if df_buildup is not None:
                    fh.write('\n')
                    fh.write('Buildup region pass rate (depth <= dmax)\n')
                    df_buildup.to_csv(fh, index=False)
                if df_past is not None:
                    fh.write('\n')
                    fh.write('Past-dmax region pass rate (depth > dmax)\n')
                    df_past.to_csv(fh, index=False)

        try:
            _write_csv(summary_csv)
            print(f"\nSummary CSV saved: {summary_csv}")
        except PermissionError:
            from datetime import datetime
            fallback = summary_csv.replace('.csv', f'_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
            _write_csv(fallback)
            print(f"\nCSV locked (close it in Excel first next time).")
            print(f"Saved to fallback: {fallback}")
        print(f"\n{'='*60}")
        print(f"  PDD pass rate  ({SHEET1_NAME} vs {SHEET2_NAME})  [{_criteria_str}]")
        print(f"{'='*60}")
        print(df_table.to_string(index=False))
        if df_buildup is not None:
            print(f"\nBuildup region pass rate (depth <= dmax)")
            print(df_buildup.to_string(index=False))
        if df_past is not None:
            print(f"\nPast-dmax region pass rate (depth > dmax)")
            print(df_past.to_string(index=False))

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
                if y < 40:
                    c.showPage()
                    c.setFont("Courier", 7)
                    y = rl_h - 36
                c.drawString(36, y, line)
                y -= 10

            def _draw_region_block(title, df_block):
                nonlocal y
                if y < 80:
                    c.showPage()
                    c.setFont("Courier", 7)
                    y = rl_h - 36
                y -= 14
                c.setFont("Helvetica-Bold", 10)
                c.drawString(36, y, title)
                y -= 12
                c.setFont("Courier", 7)
                for line in df_block.to_string(index=False).split('\n'):
                    if y < 40:
                        c.showPage()
                        c.setFont("Courier", 7)
                        y = rl_h - 36
                    c.drawString(36, y, line)
                    y -= 10

            if df_buildup is not None:
                _draw_region_block("Buildup region pass rate (depth <= dmax)", df_buildup)
            if df_past is not None:
                _draw_region_block("Past-dmax region pass rate (depth > dmax)", df_past)
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

    elif SAVE_REPORT and all_figs:
        # Plot-only / single-sheet mode: figures-only report (no stats front page).
        from reportlab.lib.pagesizes import letter as rl_letter
        from reportlab.pdfgen import canvas as rl_canvas
        from reportlab.lib.utils import ImageReader

        pdf_path = os.path.join(PDD_DIR, f"{s1}_{s2}_pdd_plots_{timestamp}.pdf")
        rl_w, rl_h = rl_letter
        c = rl_canvas.Canvas(pdf_path, pagesize=rl_letter)
        for energy_label, energy_text, fig in all_figs:
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=REPORT_DPI, bbox_inches='tight')
            buf.seek(0)
            img = ImageReader(buf)
            iw, ih = img.getSize()
            fig_w = rl_w - 72
            fig_h = fig_w * (ih / float(iw))
            top = rl_h - 36
            if top - fig_h < 18:
                fig_h = max(top - 18, 10)
                fig_w = fig_h * (iw / float(ih))
            c.drawImage(img, 36, top - fig_h, width=fig_w, height=fig_h)
            c.showPage()
            # individual per-energy PDF
            energy_pdf = os.path.join(REPORTS_DIR, f"{s1}_{s2}_{energy_label}.pdf")
            ec = rl_canvas.Canvas(energy_pdf, pagesize=rl_letter)
            ec.drawImage(img, 36, top - fig_h, width=fig_w, height=fig_h)
            ec.save()
            plt.close(fig)
        c.save()
        print(f"PDF (plots only) saved: {pdf_path}")
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
