"""
Beam Specifier Comparison Tool
------------------------------
Extracts single-value beam metrics (starting with PDD at reference depth)
from multiple sheets across one or more xlsx files and plots the
distribution as a box plot grouped by (energy, SSD), with optional
highlighting of one sheet across all groups.

Data format expected: columns FS, Pos, Dose, Axis (same as PDDCompare).
Depth-dose scans use Axis == 'Z'; Pos is depth [cm].

Filename parsing (folder mode): looks for tokens like "6X", "6FFF", "10X",
"10FFF", "15X", "8FFF" for energy and "100SSD" / "90SSD" etc. for SSD.
"""

import os
import re
import glob
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import numpy as np
import pandas as pd
from scipy import interpolate as interp
from matplotlib import pyplot as plt


# ── Metric functions ─────────────────────────────────────────────────────────

def _z_axis_sorted(df, fs):
    """Return (y, d) arrays for the Z-axis PDD at the given FS, sorted by depth."""
    sub = df[(df['FS'] == fs) & (df['Axis'] == 'Z')].sort_values('Pos')
    if sub.empty:
        return None, None
    y = sub['Pos'].to_numpy(dtype=float)
    d = sub['Dose'].to_numpy(dtype=float)
    if len(y) < 4:
        return None, None
    return y, d


def metric_pdd_at_depth(df, fs, depth_cm):
    """PDD [% of dmax] at the given depth for the given FS on the Z axis."""
    y, d = _z_axis_sorted(df, fs)
    if y is None:
        return float('nan')
    dmax = float(d.max())
    if dmax <= 0:
        return float('nan')
    try:
        val = float(interp.pchip(y, d)(depth_cm))
    except Exception:
        return float('nan')
    return val / dmax * 100.0


def metric_dmax_depth(df, fs, _depth_cm):
    """Depth of dmax [cm] via pchip on a fine grid."""
    y, d = _z_axis_sorted(df, fs)
    if y is None:
        return float('nan')
    try:
        pchip = interp.pchip(y, d)
    except Exception:
        return float('nan')
    y_lo = float(y.min())
    y_hi = min(float(y.max()), y_lo + 5.0)
    fine = np.linspace(y_lo, y_hi, 2001)
    vals = pchip(fine)
    return float(fine[int(np.argmax(vals))])


def _profile_sorted(df, fs, depth_cm, axis):
    """Return (pos, dose) arrays for a profile at the given FS, depth, and axis."""
    if 'Depth' not in df.columns:
        return None, None
    sub = df[(df['FS'] == fs) & (df['Axis'] == axis)]
    if sub.empty:
        return None, None
    available = sub['Depth'].dropna().to_numpy(dtype=float)
    if len(available) == 0:
        return None, None
    nearest = available[np.argmin(np.abs(available - depth_cm))]
    sub = sub[sub['Depth'] == nearest].sort_values('Pos')
    if len(sub) < 4:
        return None, None
    return sub['Pos'].to_numpy(dtype=float), sub['Dose'].to_numpy(dtype=float)


def _fwhm(pos, dose):
    """FWHM [same units as pos] via pchip 50% crossing on a fine grid."""
    try:
        pchip = interp.pchip(pos, dose)
    except Exception:
        return float('nan')
    d_max = float(dose.max())
    fine = np.linspace(float(pos.min()), float(pos.max()), 10001)
    fine_dose = pchip(fine)
    above = fine_dose >= d_max * 0.5
    if not above.any() or above.all():
        return float('nan')
    diff = np.diff(above.astype(int))
    up_x   = np.where(diff > 0)[0]
    down_x = np.where(diff < 0)[0]
    if len(up_x) == 0 or len(down_x) == 0:
        return float('nan')
    return float(fine[down_x[-1]]) - float(fine[up_x[0]])


def _inflection_field_size(pos, dose):
    """Field size [same units as pos] via inflection points (peak |derivative|) on each side.
    Works for both FF and FFF profiles since it finds the steepest fall-off edge
    regardless of absolute dose level at the field edge."""
    try:
        pchip = interp.pchip(pos, dose)
    except Exception:
        return float('nan')
    dpchip = pchip.derivative()
    fine = np.linspace(float(pos.min()), float(pos.max()), 10001)
    deriv = dpchip(fine)
    mid = len(fine) // 2
    # Left penumbra: derivative is positive (dose rising); pick its maximum
    left_idx = int(np.argmax(deriv[:mid]))
    # Right penumbra: derivative is negative (dose falling); pick its minimum
    right_idx = mid + int(np.argmin(deriv[mid:]))
    return float(fine[right_idx]) - float(fine[left_idx])


def metric_field_size_inflection(df, fs, depth_cm, axis='X'):
    """Measured field size (inflection-point method) for axis='X', 'Y', or 'Both'."""
    if axis == 'Both':
        fx = metric_field_size_inflection(df, fs, depth_cm, 'X')
        fy = metric_field_size_inflection(df, fs, depth_cm, 'Y')
        valid = [v for v in (fx, fy) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    return _inflection_field_size(pos, dose)


def _flatness(pos, dose):
    """Flatness within 80% of FWHM: (Dmax-Dmin)/(2*D_CAX)*100."""
    if pos is None:
        return float('nan')
    try:
        d_cax = float(interp.pchip(pos, dose)(0.0))
    except Exception:
        return float('nan')
    if d_cax <= 0:
        return float('nan')
    fwhm = _fwhm(pos, dose)
    if np.isnan(fwhm) or fwhm <= 0:
        return float('nan')
    half_eval = fwhm * 0.4
    mask = np.abs(pos) <= half_eval
    if mask.sum() < 2:
        return float('nan')
    d = dose[mask]
    lo, hi = float(d.min()), float(d.max())
    return (hi - lo) / (2.0 * d_cax) * 100.0


def metric_flatness(df, fs, depth_cm, axis='X'):
    """Profile flatness for axis='X', 'Y', or 'Both' (worst of X/Y)."""
    if axis == 'Both':
        fx = _flatness(*_profile_sorted(df, fs, depth_cm, 'X'))
        fy = _flatness(*_profile_sorted(df, fs, depth_cm, 'Y'))
        valid = [v for v in (fx, fy) if not np.isnan(v)]
        return max(valid) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    return _flatness(pos, dose)


def metric_symmetry(df, fs, depth_cm, axis='X'):
    """Symmetry: max|D(+x)-D(-x)|/D_CAX*100 within 80% of FWHM. Both=worst of X/Y."""
    if axis == 'Both':
        sx = metric_symmetry(df, fs, depth_cm, 'X')
        sy = metric_symmetry(df, fs, depth_cm, 'Y')
        valid = [v for v in (sx, sy) if not np.isnan(v)]
        return max(valid) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    try:
        pchip = interp.pchip(pos, dose)
        d_cax = float(pchip(0.0))
    except Exception:
        return float('nan')
    if d_cax <= 0:
        return float('nan')
    fwhm = _fwhm(pos, dose)
    if np.isnan(fwhm) or fwhm <= 0:
        return float('nan')
    half_eval = fwhm * 0.4
    pos_eval = pos[(pos > 0) & (pos <= half_eval)]
    if len(pos_eval) == 0:
        return float('nan')
    d_pos = pchip(pos_eval)
    d_neg = pchip(-pos_eval)
    return float(np.max(np.abs(d_pos - d_neg))) / d_cax * 100.0


def _penumbra_side(pos, dose, side):
    """80-20 penumbra width [cm] on one side. side='left' or 'right'."""
    try:
        pchip = interp.pchip(pos, dose)
    except Exception:
        return float('nan')
    d_max = float(dose.max())
    fine = np.linspace(float(pos.min()), float(pos.max()), 10001)
    fine_dose = pchip(fine)
    if side == 'right':
        f = fine[fine >= 0]
        d = fine_dose[fine >= 0]
    else:
        f = fine[fine <= 0][::-1]       # reverse so index 0 = centre, outward
        d = fine_dose[fine <= 0][::-1]
    d_80, d_20 = d_max * 0.80, d_max * 0.20
    diff_80 = np.diff((d >= d_80).astype(int))
    diff_20 = np.diff((d >= d_20).astype(int))
    down_80 = np.where(diff_80 < 0)[0]
    down_20 = np.where(diff_20 < 0)[0]
    if len(down_80) == 0 or len(down_20) == 0:
        return float('nan')
    return abs(float(f[down_20[0]]) - float(f[down_80[0]]))


def _penumbra_inflection_side(pos, dose, side, geo_half=None):
    """Penumbra width via inflection-point-referenced 80/20 surrogates.
    D_inf (dose at steepest gradient) ≈ 50% of the penumbra transition.
    Upper threshold = 1.6 × D_inf, lower = 0.4 × D_inf.
    geo_half: geometric half-field width [cm] = (FS/2)*(SSD+depth)/100.
    When provided, inflection search is restricted to within 1 cm of the geometric edge
    so FFF central slope can't be mistaken for the penumbra.
    The 80/20 crossing search is capped at 2.5 cm from the inflection to prevent
    far-field artifacts from producing spuriously large values."""
    try:
        pchip = interp.pchip(pos, dose)
    except Exception:
        return float('nan')
    fine = np.linspace(float(pos.min()), float(pos.max()), 10001)
    fine_dose = pchip(fine)
    deriv = pchip.derivative()(fine)
    mid = len(fine) // 2
    step = float(fine[1] - fine[0])
    n_search = max(1, int(2.5 / step))   # 2.5 cm cap on crossing search

    if side == 'left':
        if geo_half is not None:
            edge = -geo_half
            search = (fine >= edge - 1.0) & (fine <= edge + 1.0) & (fine < 0)
            if not search.any():
                search = fine < 0   # fallback
            inf_idx = int(np.argmax(np.where(search, deriv, -np.inf)))
        else:
            inf_idx = int(np.argmax(deriv[:mid]))
        d_inf = float(fine_dose[inf_idx])
        upper, lower = d_inf * 1.6, d_inf * 0.4
        seg = fine_dose[inf_idx : min(inf_idx + n_search, mid)]
        above = np.where(seg >= upper)[0]
        if len(above) == 0:
            return float('nan')
        x_upper = float(fine[inf_idx + above[0]])
        seg2 = fine_dose[max(0, inf_idx - n_search) : inf_idx + 1][::-1]
        below = np.where(seg2 <= lower)[0]
        if len(below) == 0:
            return float('nan')
        x_lower = float(fine[inf_idx - below[0]])
        result = float(x_upper - x_lower)

    else:  # right
        if geo_half is not None:
            edge = geo_half
            search = (fine >= edge - 1.0) & (fine <= edge + 1.0) & (fine > 0)
            if not search.any():
                search = fine > 0   # fallback
            inf_idx = int(np.argmin(np.where(search, deriv, np.inf)))
        else:
            inf_idx = mid + int(np.argmin(deriv[mid:]))
        d_inf = float(fine_dose[inf_idx])
        upper, lower = d_inf * 1.6, d_inf * 0.4
        seg = fine_dose[max(mid, inf_idx - n_search) : inf_idx + 1][::-1]
        above = np.where(seg >= upper)[0]
        if len(above) == 0:
            return float('nan')
        x_upper = float(fine[inf_idx - above[0]])
        seg2 = fine_dose[inf_idx : min(inf_idx + n_search, len(fine))]
        below = np.where(seg2 <= lower)[0]
        if len(below) == 0:
            return float('nan')
        x_lower = float(fine[inf_idx + below[0]])
        result = float(x_lower - x_upper)

    return result if 0 < result < 3.0 else float('nan')


def metric_penumbra(df, fs, depth_cm, axis='X', side='Both'):
    """80-20 penumbra width. axis: X/Y/Both; side: Left/Right/Both (avg of selected)."""
    if axis == 'Both':
        px = metric_penumbra(df, fs, depth_cm, 'X', side)
        py = metric_penumbra(df, fs, depth_cm, 'Y', side)
        valid = [v for v in (px, py) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    if side == 'Both':
        pl = _penumbra_side(pos, dose, 'left')
        pr = _penumbra_side(pos, dose, 'right')
        valid = [v for v in (pl, pr) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    return _penumbra_side(pos, dose, side.lower())


def metric_penumbra_inflection(df, fs, depth_cm, axis='X', side='Both'):
    """Inflection-referenced penumbra. axis: X/Y/Both; side: Left/Right/Both (avg of selected)."""
    ssd_vals = df['SSD'].dropna() if 'SSD' in df.columns else []
    ssd = float(ssd_vals.iloc[0]) if len(ssd_vals) > 0 else None
    geo_half = (fs / 2.0) * (ssd + depth_cm) / 100.0 if ssd is not None else None

    if axis == 'Both':
        px = metric_penumbra_inflection(df, fs, depth_cm, 'X', side)
        py = metric_penumbra_inflection(df, fs, depth_cm, 'Y', side)
        valid = [v for v in (px, py) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    if side == 'Both':
        pl = _penumbra_inflection_side(pos, dose, 'left',  geo_half)
        pr = _penumbra_inflection_side(pos, dose, 'right', geo_half)
        valid = [v for v in (pl, pr) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    return _penumbra_inflection_side(pos, dose, side.lower(), geo_half)


def metric_oar(df, fs, depth_cm, axis='X', side='Both', oar_pos=2.0):
    """Off-axis ratio D(±oar_pos)/D_CAX [%]. Both=average of + and - sides."""
    if axis == 'Both':
        ox = metric_oar(df, fs, depth_cm, 'X', side, oar_pos)
        oy = metric_oar(df, fs, depth_cm, 'Y', side, oar_pos)
        valid = [v for v in (ox, oy) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    try:
        pchip = interp.pchip(pos, dose)
        d_cax = float(pchip(0.0))
    except Exception:
        return float('nan')
    if d_cax <= 0:
        return float('nan')
    if side == 'Both':
        vals = [float(pchip(oar_pos)), float(pchip(-oar_pos))]
        return float(np.mean(vals)) / d_cax * 100.0
    elif side == 'Right':
        return float(pchip(oar_pos)) / d_cax * 100.0
    else:
        return float(pchip(-oar_pos)) / d_cax * 100.0


def metric_field_size(df, fs, depth_cm, axis='X'):
    """Measured field size (FWHM) for axis='X', 'Y', or 'Both' (average of X and Y)."""
    if axis == 'Both':
        fx = metric_field_size(df, fs, depth_cm, 'X')
        fy = metric_field_size(df, fs, depth_cm, 'Y')
        valid = [v for v in (fx, fy) if not np.isnan(v)]
        return float(np.mean(valid)) if valid else float('nan')
    pos, dose = _profile_sorted(df, fs, depth_cm, axis)
    if pos is None:
        return float('nan')
    return _fwhm(pos, dose)


def metric_tmr_20_10(df, fs, depth_cm):
    """TMR(20,10) beam quality index, derived from PDD data.

    PDD20,10 = D(20 cm) / D(10 cm) on the Z axis (the dmax normalization
    cancels in the ratio), converted via the Followill et al. (1998) fit:
        TMR20,10 = 1.2661 * PDD20,10 - 0.0595
    Conventionally a 10x10 field at SSD 100 cm; the GUI field size is used as
    given so it can be checked at other sizes. `depth_cm` is ignored — the
    depths are fixed at 10 and 20 cm.
    """
    y, d = _z_axis_sorted(df, fs)
    if y is None:
        return float('nan')
    # Both reference depths must be inside the scan — don't extrapolate.
    if float(y.min()) > 10.0 or float(y.max()) < 20.0:
        return float('nan')
    try:
        pchip = interp.pchip(y, d)
        d10 = float(pchip(10.0))
        d20 = float(pchip(20.0))
    except Exception:
        return float('nan')
    if d10 <= 0:
        return float('nan')
    return 1.2661 * (d20 / d10) - 0.0595


_ENERGIES = ["6X", "10X", "15X", "6FFF", "8FFF", "10FFF"]

METRICS = {
    "PDD at depth":     {"fn": metric_pdd_at_depth, "file_keyword": "PDD",     "unit": "%",
                         "spec": {"6X": 66.6, "10X": 73.8, "15X": 76.8,
                                  "6FFF": 63.8, "8FFF": 68.9, "10FFF": 70.9}, "tol": 0.5},
    "Depth of dmax":    {"fn": metric_dmax_depth,   "file_keyword": "PDD",     "unit": "cm",
                         "spec": {"6X": 1.43, "10X": 2.37, "15X": 2.70,
                                  "6FFF": 1.36, "8FFF": 1.79, "10FFF": 2.16}, "tol": 0.15},
    "TMR 20/10":        {"fn": metric_tmr_20_10,    "file_keyword": "PDD",     "unit": "",
                         "spec": None, "tol": None},
    "Profile Flatness":   {"fn": metric_flatness,    "file_keyword": "Profile", "unit": "%",
                           "spec": None, "tol": None, "target": 0.0, "needs_axis": True,
                           "spec_type": "energy_fs", "comparison": "max"},
    "Profile Field Size":        {"fn": metric_field_size_inflection, "file_keyword": "Profile", "unit": "cm",
                           "spec": None,
                           "spec_fn": lambda ssd, depth, fs: ((ssd + depth) / 100.0) * fs if ssd is not None else None,
                           "tol": 0.2,  "needs_axis": True, "default_axis": "Both"},
    "Profile Field Size (FWHM)": {"fn": metric_field_size,  "file_keyword": "Profile", "unit": "cm",
                           "spec": None,
                           "spec_fn": lambda ssd, depth, fs: ((ssd + depth) / 100.0) * fs if ssd is not None else None,
                           "tol": 0.2,  "needs_axis": True, "default_axis": "Both"},
    "Profile Symmetry":   {"fn": metric_symmetry,   "file_keyword": "Profile", "unit": "%",
                           "spec": 0.0,  "tol": 1.5,  "needs_axis": True},
    "Profile Penumbra":        {"fn": metric_penumbra_inflection, "file_keyword": "Profile", "unit": "cm",
                                "spec": None, "tol": None, "needs_axis": True, "needs_side": True, "default_axis": "Both"},
    "Profile Penumbra (80-20%)": {"fn": metric_penumbra,          "file_keyword": "Profile", "unit": "cm",
                                "spec": None, "tol": None, "needs_axis": True, "needs_side": True, "default_axis": "Both"},
    "Profile OAR":        {"fn": metric_oar,       "file_keyword": "Profile", "unit": "%",
                           "spec": None, "tol": 1.0, "needs_axis": True, "needs_side": True,
                           "needs_oar_pos": True, "spec_type": "energy_fs_pos"},
}


# ── Filename parsing ─────────────────────────────────────────────────────────

_ENERGY_RE = re.compile(r'(?<![A-Za-z0-9])(6X|6FFF|8FFF|10X|10FFF|15X)(?![A-Za-z0-9])', re.I)
_SSD_RE    = re.compile(r'(?<!\d)(\d{2,3})\s*SSD', re.I)


def parse_energy_ssd(name):
    """Return (energy, ssd) parsed from a filename/path, or (None, None)."""
    e = _ENERGY_RE.search(str(name))
    s = _SSD_RE.search(str(name))
    energy = e.group(1).upper() if e else None
    ssd    = int(s.group(1)) if s else None
    return energy, ssd


# ── File registry: list of {'path', 'energy', 'ssd', 'sheets'} ──────────────

loaded_files = []
_user_specs  = {}   # {metric_name: varies by spec_type}
_tps_values  = {}   # {metric_name: {energy: value}}
_df_cache    = {}   # {path: {sheet: DataFrame}} — populated by Load Data

_known_raw_detectors = set()  # populated from runs

# Default Pion corrections (multiplicative) per canonical detector / energy.
# Blank energies = 1.0 (no correction). FFF beams only — X beams have negligible Pion.
_pion_corrections = {
    "CC 04":    {"6X": 1.0, "10X": 1.0, "15X": 1.0, "6FFF": 1.000, "8FFF": 0.998, "10FFF": 0.998},
    "CC 13":    {"6X": 1.0, "10X": 1.0, "15X": 1.0, "6FFF": 0.995, "8FFF": 0.992, "10FFF": 0.990},
    "TN 31010": {"6X": 1.0, "10X": 1.0, "15X": 1.0, "6FFF": 0.999, "8FFF": 0.998, "10FFF": 0.997},
    "TN 31021": {"6X": 1.0, "10X": 1.0, "15X": 1.0, "6FFF": 0.998, "8FFF": 0.997, "10FFF": 0.996},
}

# Alias map pre-populated with common detector name variants → canonical name.
# Shown and editable in the Detector Aliases dialog. Substring-matched (case-insensitive).
_detector_aliases = {
    "CC 04":         "CC 04",
    "CC04":          "CC 04",
    "CC 13":         "CC 13",
    "CC13":          "CC 13",
    "31010":         "TN 31010",
    "Semiflex (1)":  "TN 31010",
    "31021":         "TN 31021",
    "TN31021":       "TN 31021",
    "Semiflex 3D":   "TN 31021",
}


def _canonical_detector(raw):
    """Return canonical detector name via exact then substring match in _detector_aliases."""
    raw = str(raw).strip()
    if raw in _detector_aliases:
        return _detector_aliases[raw]
    for key, canon in _detector_aliases.items():
        if key.lower() in raw.lower():
            return canon
    return raw

# OAR seed data: {energy: {fs: {oar_pos: spec_%}}}  tol ±1%
_OAR_SEED = {
    "6X":    {10.5: {2.0: 100.0, 4.0: 100.0}, 40.0: {6.0: 100.0, 18.0: 100.0}},
    "10X":   {10.5: {2.0: 100.0, 4.0: 100.0}, 40.0: {6.0: 100.0, 18.0: 100.0}},
    "15X":   {10.5: {2.0: 100.0, 4.0: 100.0}, 40.0: {6.0: 100.0, 18.0: 100.0}},
    "6FFF":  {10.5: {2.0: 97.8,  4.0: 91.1},  40.0: {6.0: 89.8,  18.0: 59.5}},
    "8FFF":  {10.5: {2.0: 96.5,  4.0: 87.7},  40.0: {6.0: 83.4,  18.0: 49.6}},
    "10FFF": {10.5: {2.0: 95.8,  4.0: 85.8},  40.0: {6.0: 80.2,  18.0: 45.6}},
}
_user_specs["Profile OAR"] = {e: {fs: dict(pos_map)
                                  for fs, pos_map in fs_map.items()}
                               for e, fs_map in _OAR_SEED.items()}

# Flatness seed data: {energy: {fs: limit_%}}  — FFF beams have no spec
_FLATNESS_SEED = {
    "6X":  {10.5: 2.75, 40.0: 2.25},
    "10X": {10.5: 2.75, 40.0: 2.25},
    "15X": {10.5: 2.75, 40.0: 2.37},
}
_user_specs["Profile Flatness"] = {e: dict(fs_map) for e, fs_map in _FLATNESS_SEED.items()}

# TPS reference values for PDD at depth — 100 SSD, 10.5×10.5 cm field, 10 cm depth
_tps_values["PDD at depth"] = {
    "6X":    66.5,
    "10X":   73.6,
    "15X":   76.5,
    "6FFF":  63.4,
    "8FFF":  68.2,
    "10FFF": 70.1,
}


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Beam Specifier Compare v0.2")

fs_var          = tk.StringVar(master=root, value="10.5")
depth_var       = tk.StringVar(master=root, value="10")
metric_var      = tk.StringVar(master=root, value="PDD at depth")
highlight_var   = tk.StringVar(master=root)
show_labels_var    = tk.BooleanVar(master=root, value=False)
show_detector_var  = tk.BooleanVar(master=root, value=True)
apply_pion_var     = tk.BooleanVar(master=root, value=False)
color_tps_var      = tk.BooleanVar(master=root, value=False)
show_pion_spec_var = tk.BooleanVar(master=root, value=True)
axis_var        = tk.StringVar(master=root, value="Both")
side_var        = tk.StringVar(master=root, value="Both")
oar_pos_var     = tk.StringVar(master=root, value="2.0")
spec_var        = tk.StringVar(master=root, value="")
tol_var         = tk.StringVar(master=root, value="")
show_spec_var   = tk.BooleanVar(master=root, value=True)
show_tps_var    = tk.BooleanVar(master=root, value=False)
tps_tol_var     = tk.StringVar(master=root, value="0.5")

main = ttk.Frame(root, padding=10)
main.grid(row=0, column=0, sticky="nsew")
root.rowconfigure(0, weight=1)
root.columnconfigure(0, weight=1)
main.rowconfigure(14, weight=1)
main.columnconfigure(1, weight=1)


def _refresh_sheets():
    """Repopulate sheet listbox with union of sheet names across loaded files."""
    sheet_listbox.delete(0, tk.END)
    all_sheets = []
    seen = set()
    for entry in loaded_files:
        for s in entry['sheets']:
            if s not in seen:
                all_sheets.append(s); seen.add(s)
    all_sheets.sort(key=str.casefold)
    for s in all_sheets:
        sheet_listbox.insert(tk.END, s)
    highlight_combo['values'] = [''] + all_sheets
    if highlight_var.get() not in all_sheets:
        highlight_var.set('')


def _refresh_filter_checks(target_frame, target_vars, values, label_fmt=str):
    """Rebuild the energy/SSD filter checkboxes."""
    for w in target_frame.winfo_children():
        w.destroy()
    target_vars.clear()
    for v in values:
        var = tk.BooleanVar(master=root, value=True)
        cb = ttk.Checkbutton(target_frame, text=label_fmt(v), variable=var)
        cb.pack(side='left', padx=4)
        target_vars[v] = var


def _refresh_filters_from_loaded():
    energies = sorted({f['energy'] for f in loaded_files if f['energy'] is not None})
    ssds     = sorted({f['ssd']    for f in loaded_files if f['ssd']    is not None})
    _refresh_filter_checks(energy_frame, energy_vars, energies)
    _refresh_filter_checks(ssd_frame,    ssd_vars,    ssds, label_fmt=lambda v: f"{v} SSD")


def _load_paths(paths):
    """Given a list of xlsx paths, read metadata and populate registry."""
    loaded_files.clear()
    _df_cache.clear()
    load_status_lbl.config(text="")
    file_listbox.delete(0, tk.END)
    for p in paths:
        try:
            xl = pd.ExcelFile(p)
            sheets = xl.sheet_names
        except Exception as e:
            file_listbox.insert(tk.END, f"[ERROR reading] {os.path.basename(p)}: {e}")
            continue
        energy, ssd = parse_energy_ssd(os.path.basename(p))
        loaded_files.append({
            'path':   p,
            'energy': energy,
            'ssd':    ssd,
            'sheets': sheets,
        })
        tag = f"{energy or '?'}  {ssd or '?'} SSD"
        file_listbox.insert(tk.END, f"{tag:<18} {os.path.basename(p)}")
    _refresh_sheets()
    _refresh_filters_from_loaded()


def choose_file():
    p = filedialog.askopenfilename(
        filetypes=[("Excel", "*.xlsx *.xls"), ("All files", "*.*")]
    )
    if not p:
        return
    _load_paths([p])


def choose_folder():
    d = filedialog.askdirectory()
    if not d:
        return
    paths = sorted(glob.glob(os.path.join(d, "*.xlsx")))
    if not paths:
        messagebox.showwarning("No files", "No .xlsx files found in that folder.")
        return
    _load_paths(paths)


# Row 0: file buttons
ttk.Button(main, text="Browse File…",   command=choose_file  ).grid(row=0, column=0, sticky="w")
ttk.Button(main, text="Browse Folder…", command=choose_folder).grid(row=0, column=1, sticky="w", padx=5)

# Row 1: file list
ttk.Label(main, text="Loaded files:").grid(row=1, column=0, sticky="nw", pady=(10, 0))
file_listbox = tk.Listbox(main, height=6, width=70)
file_listbox.grid(row=1, column=1, columnspan=2, sticky="w", pady=(10, 0))

# Row 2: energy filters
ttk.Label(main, text="Energies:").grid(row=2, column=0, sticky="w", pady=(10, 0))
energy_frame = ttk.Frame(main)
energy_frame.grid(row=2, column=1, columnspan=2, sticky="w", pady=(10, 0))
energy_vars = {}

# Row 3: SSD filters
ttk.Label(main, text="SSDs:").grid(row=3, column=0, sticky="w")
ssd_frame = ttk.Frame(main)
ssd_frame.grid(row=3, column=1, columnspan=2, sticky="w")
ssd_vars = {}

# Row 4: sheet multi-select
ttk.Label(main, text="Sheets to compare:").grid(row=4, column=0, sticky="nw", pady=(10, 0))
sheet_frame = ttk.Frame(main)
sheet_frame.grid(row=4, column=1, columnspan=2, sticky="w", pady=(10, 0))
sheet_listbox = tk.Listbox(sheet_frame, selectmode=tk.EXTENDED, height=10, width=40, exportselection=False)
sheet_listbox.grid(row=0, column=0)
sheet_scroll = ttk.Scrollbar(sheet_frame, orient="vertical", command=sheet_listbox.yview)
sheet_scroll.grid(row=0, column=1, sticky="ns")
sheet_listbox.configure(yscrollcommand=sheet_scroll.set)

# Row 5: metric
ttk.Label(main, text="Metric:").grid(row=5, column=0, sticky="w", pady=(10, 0))
metric_combo = ttk.Combobox(main, textvariable=metric_var, values=list(METRICS.keys()),
                            state="readonly", width=20)
metric_combo.grid(row=5, column=1, sticky="w", padx=5, pady=(10, 0))

# Row 6: FS entry + Spec
ttk.Label(main, text="Reference field size [cm]:").grid(row=6, column=0, sticky="w")
ttk.Entry(main, textvariable=fs_var, width=8).grid(row=6, column=1, sticky="w", padx=5)
spec_single_frame = ttk.Frame(main)
spec_single_frame.grid(row=6, column=2, sticky="w", padx=5)
ttk.Label(spec_single_frame, text="Spec:").pack(side='left')
ttk.Entry(spec_single_frame, textvariable=spec_var, width=8).pack(side='left', padx=4)

spec_multi_frame = ttk.Frame(main)
spec_multi_frame.grid(row=6, column=2, sticky="w", padx=5)
ttk.Label(spec_multi_frame, text="Spec:").pack(side='left')
# Button command assigned after _open_spec_editor is defined below
spec_multi_btn = ttk.Button(spec_multi_frame, text="Per energy…")
spec_multi_btn.pack(side='left', padx=4)
spec_multi_frame.grid_remove()

spec_2d_frame = ttk.Frame(main)
spec_2d_frame.grid(row=6, column=2, sticky="w", padx=5)
ttk.Label(spec_2d_frame, text="Spec:").pack(side='left')
# Button command assigned after _open_spec_editor_2d is defined below
spec_2d_btn = ttk.Button(spec_2d_frame, text="Per energy/FS…")
spec_2d_btn.pack(side='left', padx=4)
spec_2d_frame.grid_remove()

spec_ef_frame = ttk.Frame(main)
spec_ef_frame.grid(row=6, column=2, sticky="w", padx=5)
ttk.Label(spec_ef_frame, text="Tolerance:").pack(side='left')
# Button command assigned after _open_spec_editor_ef is defined below
spec_ef_btn = ttk.Button(spec_ef_frame, text="Per energy/FS…")
spec_ef_btn.pack(side='left', padx=4)
spec_ef_frame.grid_remove()


# Row 7: depth entry + Tolerance
ttk.Label(main, text="Reference depth [cm]:").grid(row=7, column=0, sticky="w")
ttk.Entry(main, textvariable=depth_var, width=8).grid(row=7, column=1, sticky="w", padx=5)
tol_frame = ttk.Frame(main)
tol_frame.grid(row=7, column=2, sticky="w", padx=5)
tol_label_w  = ttk.Label(tol_frame, text="Tolerance ±:")
tol_label_w.pack(side='left')
tol_entry_w  = ttk.Entry(tol_frame, textvariable=tol_var, width=8)
tol_entry_w.pack(side='left', padx=4)
tol_spec_lbl = ttk.Label(tol_frame, text="")   # "Spec: 0 %" shown for limit-type metrics
show_spec_chk = ttk.Checkbutton(tol_frame, text="Show spec/tol lines", variable=show_spec_var)
show_spec_chk.pack(side='left', padx=10)
show_tps_chk  = ttk.Checkbutton(tol_frame, text="Show TPS", variable=show_tps_var)
show_tps_chk.pack(side='left', padx=10)
ttk.Label(tol_frame, text="TPS tol ±:").pack(side='left')
ttk.Entry(tol_frame, textvariable=tps_tol_var, width=6).pack(side='left', padx=4)
tps_btn = ttk.Button(tol_frame, text="TPS per energy…")
tps_btn.pack(side='left', padx=4)

# Row 8: profile axis selector (shown only for axis-aware metrics)
axis_label = ttk.Label(main, text="Profile axis:")
axis_label.grid(row=8, column=0, sticky="w")
axis_frame = ttk.Frame(main)
axis_frame.grid(row=8, column=1, columnspan=2, sticky="w")
for _txt, _val in (("X", "X"), ("Y", "Y"), ("Both", "Both")):
    ttk.Radiobutton(axis_frame, text=_txt, variable=axis_var, value=_val).pack(side='left', padx=6)

# Row 9: profile side selector (shown only for metrics with needs_side)
side_label = ttk.Label(main, text="Side:")
side_label.grid(row=9, column=0, sticky="w")
side_frame = ttk.Frame(main)
side_frame.grid(row=9, column=1, columnspan=2, sticky="w")
for _txt, _val in (("Left", "Left"), ("Right", "Right"), ("Both", "Both")):
    ttk.Radiobutton(side_frame, text=_txt, variable=side_var, value=_val).pack(side='left', padx=6)
side_label.grid_remove()
side_frame.grid_remove()

# Row 10: off-axis position entry (shown only for OAR metric)
oar_pos_frame = ttk.Frame(main)
oar_pos_frame.grid(row=10, column=0, columnspan=2, sticky="w")
ttk.Label(oar_pos_frame, text="Off-axis position [cm]:").pack(side='left')
ttk.Entry(oar_pos_frame, textvariable=oar_pos_var, width=8).pack(side='left', padx=5)
oar_pos_frame.grid_remove()


def _open_spec_editor():
    metric = metric_var.get()
    defaults = METRICS[metric].get("spec", {})
    if metric not in _user_specs:
        _user_specs[metric] = dict(defaults)

    dlg = tk.Toplevel(root)
    dlg.title(f"Spec per energy — {metric}")
    dlg.resizable(False, False)
    dlg.grab_set()

    entries = {}
    for i, energy in enumerate(_ENERGIES):
        ttk.Label(dlg, text=energy, width=8).grid(row=i, column=0, padx=12, pady=4, sticky='w')
        val = _user_specs[metric].get(energy)
        var = tk.StringVar(value="" if val is None else str(val))
        ttk.Entry(dlg, textvariable=var, width=10).grid(row=i, column=1, padx=12, pady=4)
        entries[energy] = var

    def _save():
        for e, v in entries.items():
            txt = v.get().strip()
            try:
                _user_specs[metric][e] = float(txt) if txt else None
            except ValueError:
                _user_specs[metric][e] = None
        dlg.destroy()

    n = len(_ENERGIES)
    ttk.Button(dlg, text="OK",     command=_save       ).grid(row=n, column=0, pady=8)
    ttk.Button(dlg, text="Cancel", command=dlg.destroy ).grid(row=n, column=1, pady=8)

spec_multi_btn.configure(command=_open_spec_editor)


def _open_tps_editor():
    metric = metric_var.get()
    if metric not in _tps_values:
        _tps_values[metric] = {}

    dlg = tk.Toplevel(root)
    dlg.title(f"TPS values per energy — {metric}")
    dlg.resizable(False, False)
    dlg.grab_set()

    # For PDD metrics show a derived column: D(SAD=100, d=10) = %dd(10)/100 × (110/100)²
    _is_pdd = "PDD" in metric or "pdd" in metric
    _ISF = (110 / 100) ** 2   # SSD=100, d=10 → SAD=100, d=10

    unit = METRICS.get(metric, {}).get("unit", "")
    header = f"TPS %dd(10)" if _is_pdd else (f"TPS value [{unit}]" if unit else "TPS value")
    ttk.Label(dlg, text="Energy",  width=8 ).grid(row=0, column=0, padx=12, pady=4, sticky='w')
    ttk.Label(dlg, text=header,    width=14).grid(row=0, column=1, padx=12, pady=4, sticky='w')
    if _is_pdd:
        ttk.Label(dlg, text="D(SAD=100, d=10)\n[cGy/MU @ 1cGy/MU cal]",
                  width=22).grid(row=0, column=2, padx=12, pady=4, sticky='w')

    entries = {}
    sad_labels = {}
    for i, energy in enumerate(_ENERGIES):
        ttk.Label(dlg, text=energy, width=8).grid(row=i + 1, column=0, padx=12, pady=4, sticky='w')
        val = _tps_values[metric].get(energy)
        var = tk.StringVar(value="" if val is None else str(val))
        ttk.Entry(dlg, textvariable=var, width=12).grid(row=i + 1, column=1, padx=12, pady=4)
        entries[energy] = var

        if _is_pdd:
            sad_lbl = ttk.Label(dlg, text="", width=22, foreground='steelblue')
            sad_lbl.grid(row=i + 1, column=2, padx=12, pady=4, sticky='w')
            sad_labels[energy] = sad_lbl

            def _update_sad(_, v=var, lbl=sad_lbl):
                try:
                    pdd = float(v.get())
                    lbl.config(text=f"{pdd / 100 * _ISF:.4f}")
                except (ValueError, ZeroDivisionError):
                    lbl.config(text="")
            var.trace_add('write', _update_sad)
            _update_sad(None)   # populate on open

    def _save():
        for e, v in entries.items():
            txt = v.get().strip()
            try:
                _tps_values[metric][e] = float(txt) if txt else None
            except ValueError:
                _tps_values[metric][e] = None
        dlg.destroy()

    n = len(_ENERGIES)
    ttk.Button(dlg, text="OK",     command=_save       ).grid(row=n + 1, column=0, pady=8)
    ttk.Button(dlg, text="Cancel", command=dlg.destroy ).grid(row=n + 1, column=1, pady=8)

tps_btn.configure(command=_open_tps_editor)


def _open_spec_editor_2d():
    """2-D spec editor: energy (rows) × FS (columns) for the current OAR position."""
    metric  = metric_var.get()
    try:
        cur_pos = float(oar_pos_var.get())
    except ValueError:
        cur_pos = 2.0

    if metric not in _user_specs:
        _user_specs[metric] = {}

    dlg = tk.Toplevel(root)
    dlg.title(f"OAR spec at ±{cur_pos:g} cm — per energy / FS")
    dlg.resizable(True, True)
    dlg.grab_set()

    # ── FS selector ────────────────────────────────────────────────────────────
    top = ttk.Frame(dlg, padding=8)
    top.pack(fill='x')
    ttk.Label(top, text="Field sizes (comma-separated cm):").pack(side='left')

    # Gather FS values already stored for this oar_pos
    existing_fs: set[float] = set()
    for e_dict in _user_specs[metric].values():
        for fs, pos_map in e_dict.items():
            if cur_pos in pos_map:
                existing_fs.add(fs)
    fs_var_dlg = tk.StringVar(value=', '.join(str(f) for f in sorted(existing_fs)) or '10.5, 40')
    ttk.Entry(top, textvariable=fs_var_dlg, width=28).pack(side='left', padx=6)

    grid_frame = ttk.Frame(dlg, padding=8)
    grid_frame.pack(fill='both', expand=True)

    entries: dict[tuple, tk.StringVar] = {}

    def rebuild_grid():
        for w in grid_frame.winfo_children():
            w.destroy()
        entries.clear()
        try:
            fs_vals = [float(x.strip()) for x in fs_var_dlg.get().split(',') if x.strip()]
        except ValueError:
            return
        ttk.Label(grid_frame, text="Energy \\ FS", width=10).grid(row=0, column=0, padx=6, pady=3)
        for j, fs in enumerate(fs_vals):
            ttk.Label(grid_frame, text=f"{fs:g} cm", width=9).grid(row=0, column=j+1, padx=6, pady=3)
        for i, energy in enumerate(_ENERGIES):
            ttk.Label(grid_frame, text=energy, width=10).grid(row=i+1, column=0, padx=6, pady=3)
            for j, fs in enumerate(fs_vals):
                existing = (_user_specs[metric].get(energy, {}).get(fs, {}).get(cur_pos))
                var = tk.StringVar(value="" if existing is None else str(existing))
                ttk.Entry(grid_frame, textvariable=var, width=9).grid(row=i+1, column=j+1, padx=6, pady=3)
                entries[(energy, fs)] = var

    ttk.Button(top, text="Update grid", command=rebuild_grid).pack(side='left', padx=4)
    rebuild_grid()

    def _save():
        for (energy, fs), var in entries.items():
            txt = var.get().strip()
            _user_specs[metric].setdefault(energy, {}).setdefault(fs, {})
            try:
                _user_specs[metric][energy][fs][cur_pos] = float(txt) if txt else None
            except ValueError:
                _user_specs[metric][energy][fs][cur_pos] = None
        dlg.destroy()

    btn_row = ttk.Frame(dlg, padding=8)
    btn_row.pack()
    ttk.Button(btn_row, text="OK",     command=_save      ).pack(side='left', padx=8)
    ttk.Button(btn_row, text="Cancel", command=dlg.destroy).pack(side='left', padx=8)


spec_2d_btn.configure(command=_open_spec_editor_2d)


def _open_spec_editor_ef():
    """2-D limit editor: energy (rows) × FS (columns) — no oar_pos dimension."""
    metric = metric_var.get()
    if metric not in _user_specs:
        _user_specs[metric] = {}

    dlg = tk.Toplevel(root)
    dlg.title(f"Tolerance per energy / FS — {metric}")
    dlg.resizable(True, True)
    dlg.grab_set()

    top = ttk.Frame(dlg, padding=8)
    top.pack(fill='x')
    ttk.Label(top, text="Field sizes (comma-separated cm):").pack(side='left')

    existing_fs: set[float] = set()
    for fs_map in _user_specs[metric].values():
        existing_fs.update(fs_map.keys())
    fs_var_dlg = tk.StringVar(value=', '.join(str(f) for f in sorted(existing_fs)) or '10.5, 40')
    ttk.Entry(top, textvariable=fs_var_dlg, width=28).pack(side='left', padx=6)

    grid_frame = ttk.Frame(dlg, padding=8)
    grid_frame.pack(fill='both', expand=True)

    entries: dict[tuple, tk.StringVar] = {}

    def rebuild_grid():
        for w in grid_frame.winfo_children():
            w.destroy()
        entries.clear()
        try:
            fs_vals = [float(x.strip()) for x in fs_var_dlg.get().split(',') if x.strip()]
        except ValueError:
            return
        ttk.Label(grid_frame, text="Energy \\ FS", width=10).grid(row=0, column=0, padx=6, pady=3)
        for j, fs in enumerate(fs_vals):
            ttk.Label(grid_frame, text=f"{fs:g} cm", width=9).grid(row=0, column=j+1, padx=6, pady=3)
        for i, energy in enumerate(_ENERGIES):
            ttk.Label(grid_frame, text=energy, width=10).grid(row=i+1, column=0, padx=6, pady=3)
            for j, fs in enumerate(fs_vals):
                existing = _user_specs[metric].get(energy, {}).get(fs)
                var = tk.StringVar(value="" if existing is None else str(existing))
                ttk.Entry(grid_frame, textvariable=var, width=9).grid(row=i+1, column=j+1, padx=6, pady=3)
                entries[(energy, fs)] = var

    ttk.Button(top, text="Update grid", command=rebuild_grid).pack(side='left', padx=4)
    rebuild_grid()

    def _save():
        for (energy, fs), var in entries.items():
            txt = var.get().strip()
            _user_specs[metric].setdefault(energy, {})
            try:
                _user_specs[metric][energy][fs] = float(txt) if txt else None
            except ValueError:
                _user_specs[metric][energy][fs] = None
        dlg.destroy()

    btn_row = ttk.Frame(dlg, padding=8)
    btn_row.pack()
    ttk.Button(btn_row, text="OK",     command=_save      ).pack(side='left', padx=8)
    ttk.Button(btn_row, text="Cancel", command=dlg.destroy).pack(side='left', padx=8)


spec_ef_btn.configure(command=_open_spec_editor_ef)


def _on_metric_change(*_):
    info = METRICS.get(metric_var.get(), {})
    if info.get("needs_axis"):
        axis_label.grid()
        axis_frame.grid()
        if "default_axis" in info:
            axis_var.set(info["default_axis"])
    else:
        axis_label.grid_remove()
        axis_frame.grid_remove()
    if info.get("needs_side"):
        side_label.grid()
        side_frame.grid()
    else:
        side_label.grid_remove()
        side_frame.grid_remove()
    if info.get("needs_oar_pos"):
        oar_pos_frame.grid()
    else:
        oar_pos_frame.grid_remove()
    spec_info = info.get("spec")
    spec_type = info.get("spec_type")
    for _f in [spec_single_frame, spec_multi_frame, spec_2d_frame, spec_ef_frame]:
        _f.grid_remove()

    if spec_type == "energy_fs_pos":
        spec_2d_frame.grid()
        spec_var.set("")
    elif spec_type == "energy_fs":
        spec_ef_frame.grid()
        spec_var.set("")
    elif info.get("spec_fn") is not None:
        spec_single_frame.grid()
        try:
            _preview = info["spec_fn"](100, float(depth_var.get()), float(fs_var.get()))
            spec_var.set(f"{_preview:.3f}")
        except Exception:
            spec_var.set("")
    elif isinstance(spec_info, dict):
        spec_multi_frame.grid()
        spec_var.set("")
    else:
        spec_single_frame.grid()
        spec_var.set("" if spec_info is None else str(spec_info))
    # Tol row: limit-type metrics show target; abs-type metrics show scalar tolerance entry
    for _w in (tol_label_w, tol_entry_w, tol_spec_lbl, show_spec_chk):
        _w.pack_forget()
    if info.get("comparison") == "max":
        _target = info.get("target")
        tol_spec_lbl.config(
            text=f"Spec: {_target:g} {info.get('unit', '')}" if _target is not None else "Spec: —"
        )
        tol_spec_lbl.pack(side='left', padx=4)
        show_spec_chk.config(text="Show tol lines")
        tol_var.set("")
    else:
        tol_label_w.pack(side='left')
        tol_entry_w.pack(side='left', padx=4)
        show_spec_chk.config(text="Show spec/tol lines")
        tol_var.set("" if info.get("tol") is None else str(info["tol"]))
    show_spec_chk.pack(side='left', padx=10)


def _update_spec_preview(*_):
    fn = METRICS.get(metric_var.get(), {}).get("spec_fn")
    if fn is None:
        return
    try:
        spec_var.set(f"{fn(100, float(depth_var.get()), float(fs_var.get())):.3f}")
    except Exception:
        pass

fs_var.trace_add('write', _update_spec_preview)
depth_var.trace_add('write', _update_spec_preview)

metric_combo.bind("<<ComboboxSelected>>", _on_metric_change)
_on_metric_change()   # initialize visibility

# Row 11: highlight + show all labels
ttk.Label(main, text="Highlight sheet:").grid(row=11, column=0, sticky="w")
highlight_combo = ttk.Combobox(main, textvariable=highlight_var, values=[''],
                               state="readonly", width=20)
highlight_combo.grid(row=11, column=1, sticky="w", padx=5)
ttk.Checkbutton(main, text="Label all points", variable=show_labels_var).grid(
    row=11, column=2, sticky="w", padx=5)
ttk.Checkbutton(main, text="Show detector/SN per sheet", variable=show_detector_var).grid(
    row=11, column=3, sticky="w", padx=5)
ttk.Checkbutton(main, text="Color TPS vs measured", variable=color_tps_var).grid(
    row=11, column=4, sticky="w", padx=5)

# Row 12: Detector aliases + Pion correction
detector_row = ttk.Frame(main)
detector_row.grid(row=12, column=0, columnspan=4, sticky="w", pady=(4, 0))
detector_alias_btn = ttk.Button(detector_row, text="Detector Aliases…")
detector_alias_btn.pack(side="left", padx=(0, 8))
pion_btn = ttk.Button(detector_row, text="Pion Corrections…")
pion_btn.pack(side="left", padx=(0, 8))
ttk.Checkbutton(detector_row, text="Apply Pion correction", variable=apply_pion_var).pack(side="left")
ttk.Checkbutton(detector_row, text="Show pion-scaled spec", variable=show_pion_spec_var).pack(side="left", padx=(8, 0))

# Row 13: Run button (above output so it's always visible)
# Row 13 run button is bound after run_compare is defined, below.

# Row 14: output
output_frame = ttk.Frame(main)
output_frame.grid(row=14, column=0, columnspan=4, sticky="nsew", pady=(6, 0))
output_frame.rowconfigure(0, weight=1)
output_frame.columnconfigure(0, weight=1)
output_text = tk.Text(output_frame, height=15, width=110, wrap='none',
                      font=('Courier', 9))
output_text.grid(row=0, column=0, sticky="nsew")
output_vsb = ttk.Scrollbar(output_frame, orient="vertical",   command=output_text.yview)
output_hsb = ttk.Scrollbar(output_frame, orient="horizontal", command=output_text.xview)
output_text.configure(yscrollcommand=output_vsb.set, xscrollcommand=output_hsb.set)
output_vsb.grid(row=0, column=1, sticky="ns")
output_hsb.grid(row=1, column=0, sticky="ew")


def run_compare():
    if not loaded_files:
        messagebox.showerror("Error", "Load at least one file first.")
        return

    sel_idx = sheet_listbox.curselection()
    if not sel_idx:
        messagebox.showerror("Error", "Select at least one sheet.")
        return
    sheets = [sheet_listbox.get(i) for i in sel_idx]

    try:
        fs_ref = float(fs_var.get())
    except ValueError:
        messagebox.showerror("Error", "Reference field size must be numeric.")
        return
    try:
        depth_ref = float(depth_var.get())
    except ValueError:
        messagebox.showerror("Error", "Reference depth must be numeric.")
        return

    allowed_energies = {e for e, v in energy_vars.items() if v.get()}
    allowed_ssds     = {s for s, v in ssd_vars.items()    if v.get()}
    # If no filters populated (single-file with unknown energy/SSD) — don't filter.
    has_energy_filter = bool(energy_vars)
    has_ssd_filter    = bool(ssd_vars)

    metric_name = metric_var.get()
    metric_info = METRICS[metric_name]
    file_kw     = metric_info.get("file_keyword")
    unit        = metric_info.get("unit", "")
    highlight   = highlight_var.get().strip() or None

    # Build (suffix, callable) combos — one per axis×side combination to run
    axis_list = (["X", "Y"] if metric_info.get("needs_axis") and axis_var.get() == "Both"
                 else [axis_var.get()] if metric_info.get("needs_axis") else [None])
    side_list = (["Left", "Right"] if metric_info.get("needs_side") and side_var.get() == "Both"
                 else [side_var.get()] if metric_info.get("needs_side") else [None])
    try:
        oar_p = float(oar_pos_var.get()) if metric_info.get("needs_oar_pos") else None
    except ValueError:
        oar_p = 2.0

    def _make_combo_fn(ax, sd, p):
        def _fn(df, fs, d):
            if ax is None:
                return metric_info["fn"](df, fs, d)
            if sd is None:
                return metric_info["fn"](df, fs, d, ax)
            if p is None:
                return metric_info["fn"](df, fs, d, ax, sd)
            return metric_info["fn"](df, fs, d, ax, sd, p)
        return _fn

    combos = []
    for ax in axis_list:
        for sd in side_list:
            parts = ([ax] if len(axis_list) > 1 else []) + ([sd] if len(side_list) > 1 else [])
            suffix = ("_" + "_".join(parts)) if parts else ""
            combos.append((suffix, _make_combo_fn(ax, sd, oar_p)))
    expand = len(combos) > 1

    axis_str = ""
    if metric_info.get("needs_axis"):
        axis_str += f", axis={axis_var.get()}"
    if metric_info.get("needs_side"):
        axis_str += f", side={side_var.get()}"
    if oar_p is not None:
        axis_str += f", pos=±{oar_p}cm"

    # results[(energy, ssd)] = { sheet: value }
    results = {}

    spec_info      = metric_info.get("spec")
    spec_type_key  = metric_info.get("spec_type")
    spec_fn        = metric_info.get("spec_fn")
    comparison     = metric_info.get("comparison", "abs")
    use_energy_fs_pos_spec = (spec_type_key == "energy_fs_pos")
    use_energy_fs_spec     = (spec_type_key == "energy_fs")
    use_per_energy_spec    = isinstance(spec_info, dict) and not use_energy_fs_pos_spec and not use_energy_fs_spec
    if use_energy_fs_pos_spec:
        energy_fs_pos_data = _user_specs.get(metric_name, {})
        spec_val = None
    elif use_energy_fs_spec:
        energy_fs_data = _user_specs.get(metric_name, {})
        spec_val = None
    elif use_per_energy_spec:
        energy_specs = _user_specs.get(metric_name, spec_info)
        spec_val = None
    elif spec_fn is not None:
        spec_val = None   # computed per group
    else:
        try:
            spec_val = float(spec_var.get()) if spec_var.get().strip() else None
        except ValueError:
            spec_val = None
    try:
        tol_val = float(tol_var.get()) if tol_var.get().strip() else None
    except ValueError:
        tol_val = None
    try:
        tps_tol_val = float(tps_tol_var.get()) if tps_tol_var.get().strip() else None
    except ValueError:
        tps_tol_val = None

    output_text.delete("1.0", tk.END)
    output_text.insert(tk.END, f"{metric_name}  (FS={fs_ref} cm, depth={depth_ref} cm{axis_str})\n")
    output_text.insert(tk.END, "=" * 70 + "\n")

    # Accumulate sheet → (sn, detectors) for the summary
    _sheet_detector_map = {}   # {sheet: {'sn': str, 'detectors': set}}

    for entry in loaded_files:
        energy = entry['energy']; ssd = entry['ssd']
        if use_energy_fs_pos_spec:
            spec_val = (energy_fs_pos_data.get(energy, {}).get(fs_ref, {}).get(oar_p)
                        if oar_p is not None else None)
        elif use_energy_fs_spec:
            spec_val = energy_fs_data.get(energy, {}).get(fs_ref)
        elif use_per_energy_spec:
            spec_val = energy_specs.get(energy)
        elif spec_fn is not None:
            spec_val = spec_fn(ssd, depth_ref, fs_ref)
        if has_energy_filter and energy is not None and energy not in allowed_energies:
            continue
        if has_ssd_filter and ssd is not None and ssd not in allowed_ssds:
            continue
        # Skip files that don't apply to this metric (e.g. profile files for PDD metric)
        if file_kw and file_kw.lower() not in os.path.basename(entry['path']).lower():
            continue
        group_key = (energy, ssd)
        group_label = f"{energy or '?'}  {ssd or '?'} SSD"
        results.setdefault(group_key, {})

        output_text.insert(tk.END, f"\n{group_label}   ({os.path.basename(entry['path'])})\n")
        output_text.insert(tk.END, "-" * 70 + "\n")
        output_text.update_idletasks()

        wanted = [s for s in sheets if s in entry['sheets']]
        if not wanted:
            continue
        cached = _df_cache.get(entry['path'])
        if cached is not None:
            dfs = {s: cached[s] for s in wanted if s in cached}
        else:
            try:
                dfs = pd.read_excel(entry['path'], sheet_name=wanted)
            except Exception as e:
                output_text.insert(tk.END, f"  ERROR reading file: {e}\n")
                continue

        for s in wanted:
            df = dfs.get(s)
            if df is None:
                continue
            try:
                df['FS'] = df['FS'].astype(float)
            except Exception:
                pass

            if show_detector_var.get():
                sub = df[df['FS'] == fs_ref]
                if file_kw == 'PDD':
                    sub = sub[sub['Axis'] == 'Z']
                elif file_kw == 'Profile' and 'Depth' in sub.columns:
                    avail = sub['Depth'].dropna().to_numpy(dtype=float)
                    if len(avail) > 0:
                        nearest = avail[np.argmin(np.abs(avail - depth_ref))]
                        sub = sub[sub['Depth'] == nearest]
                sn_vals = ([str(v) for v in sub['SN'].dropna().unique()
                            if str(v).strip() not in ("", "nan")]
                           if 'SN' in sub.columns else [])
                det_vals = []
                for col in ("Detector", "Chamber"):
                    if col in sub.columns:
                        vals = [str(v) for v in sub[col].dropna().unique()
                                if str(v).strip() not in ("", "nan")]
                        if vals:
                            det_vals = vals  # use first column that has data
                            break
                _known_raw_detectors.update(det_vals)
                canonicals_for_sheet = list({_canonical_detector(d) for d in det_vals}) if det_vals else []
                _sheet_detector_map[s] = {
                    'sn': sn_vals[0] if sn_vals else None,
                    'detectors': set(det_vals),
                    'canonical': canonicals_for_sheet[0] if canonicals_for_sheet else None,
                }

            for suffix, combo_fn in combos:
                d_arg = (spec_val if metric_name == "Depth of dmax" and spec_val is not None
                         else depth_ref)
                val = combo_fn(df, fs_ref, d_arg)
                if apply_pion_var.get() and not np.isnan(val):
                    canon = _sheet_detector_map.get(s, {}).get('canonical')
                    if canon:
                        pion = _pion_corrections.get(canon, {}).get(energy)
                        if pion is not None:
                            val *= pion
                key = f"{s}{suffix}"
                results[group_key][key] = val
                flag = ""
                if show_spec_var.get() and spec_val is not None and not np.isnan(val):
                    if comparison == "max":
                        out_of_tol = val > spec_val
                    else:
                        out_of_tol = tol_val is not None and abs(val - spec_val) > tol_val
                    if out_of_tol:
                        if metric_name == "Depth of dmax":
                            pdd_at_spec = metric_pdd_at_depth(df, fs_ref, spec_val)
                            flag = (f"  ***  (PDD at spec dmax {spec_val:.2f} cm: {pdd_at_spec:.1f}%)"
                                    if not np.isnan(pdd_at_spec) else "  ***")
                        else:
                            flag = "  ***"
                w = 18 if expand else 15
                tps_str = ""
                if show_tps_var.get() and not np.isnan(val):
                    _tv = _tps_values.get(metric_name, {}).get(energy)
                    if _tv is not None and _tv != 0:
                        tps_str = f"  Δ TPS: {(val - _tv) / _tv * 100:+.2f}%"
                output_text.insert(tk.END, f"  {key:<{w}}  {val:.3f}{tps_str}{flag}\n")

    # Keep a stable order
    def _group_sort_key(k):
        e, s = k
        # Group by energy first (preferred order), then SSD.
        pref = ['6X', '6FFF', '8FFF', '10X', '10FFF', '15X']
        return (pref.index(e) if e in pref else 99, str(e), s if s is not None else -1)
    group_keys = sorted(results.keys(), key=_group_sort_key)

    # Drop empty groups
    group_keys = [k for k in group_keys if any(not np.isnan(v) for v in results[k].values())]
    if not group_keys:
        messagebox.showwarning("No data", "No valid metric values computed.")
        return

    # ── summary stats table (per group) ──
    output_text.insert(tk.END, "\n" + "=" * 104 + "\n")
    output_text.insert(tk.END, "Summary statistics\n")
    output_text.insert(tk.END, "=" * 104 + "\n")
    _show_tps_summary = show_tps_var.get() and bool(_tps_values.get(metric_name))
    _tps_col = f" {'Δ TPS%':>8}" if _show_tps_summary else ""
    _sep_w = 104 + (9 if _show_tps_summary else 0)
    hdr = (f"{'Group':<18} {'n':>3} {'mean':>8} {'median':>8} {'std':>7} "
           f"{'P25':>8} {'P75':>8} {'min':>8} {'max':>8} {'range':>7} {'CoV%':>6}{_tps_col}")
    output_text.insert(tk.END, hdr + "\n")
    output_text.insert(tk.END, "-" * _sep_w + "\n")
    for k in group_keys:
        e, s = k
        label = f"{e or '?'} {s or '?'} SSD"
        vals = np.array([v for v in results[k].values() if not np.isnan(v)], dtype=float)
        if vals.size == 0:
            continue
        mean   = float(vals.mean())
        median = float(np.median(vals))
        std    = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
        p25    = float(np.percentile(vals, 25))
        p75    = float(np.percentile(vals, 75))
        vmin   = float(vals.min()); vmax = float(vals.max())
        rng_v  = vmax - vmin
        cov    = (std / mean * 100) if mean != 0 else float('nan')
        tps_summary_str = ""
        if _show_tps_summary:
            _tv = _tps_values.get(metric_name, {}).get(e)
            if _tv is not None and _tv != 0:
                tps_summary_str = f" {(mean - _tv) / _tv * 100:>+8.2f}"
            else:
                tps_summary_str = f" {'—':>8}"
        output_text.insert(tk.END,
            f"{label:<18} {vals.size:>3d} {mean:>8.3f} {median:>8.3f} {std:>7.3f} "
            f"{p25:>8.3f} {p75:>8.3f} {vmin:>8.3f} {vmax:>8.3f} {rng_v:>7.3f} {cov:>6.2f}"
            f"{tps_summary_str}\n"
        )

    # Highlighted sheet position within each group
    if highlight:
        output_text.insert(tk.END, "\n" + "-" * _sep_w + "\n")
        output_text.insert(tk.END, f"Position of '{highlight}' within each group:\n")
        for k in group_keys:
            vals_all = results[k]
            hl_keys = [key for key in vals_all
                       if key == highlight or key.startswith(highlight + "_")]
            if not hl_keys:
                continue
            e, s = k
            label = f"{e or '?'} {s or '?'} SSD"
            for hl_key in hl_keys:
                v = vals_all[hl_key]
                if np.isnan(v):
                    continue
                others = np.array([x for s_, x in vals_all.items()
                                   if s_ != hl_key and not np.isnan(x)], dtype=float)
                if others.size == 0:
                    continue
                mean  = float(others.mean())
                std   = float(others.std(ddof=1)) if others.size > 1 else 0.0
                z     = (v - mean) / std if std > 0 else float('nan')
                delta = v - mean
                output_text.insert(tk.END,
                    f"  {label:<18} [{hl_key}]  value={v:.3f}  Δ(mean others)={delta:+.3f}  "
                    f"z={z:+.2f}σ\n"
                )

    # ── Detector/SN summary ──
    if show_detector_var.get() and _sheet_detector_map:
        output_text.insert(tk.END, "\n" + "=" * 60 + "\n")
        output_text.insert(tk.END, "Detectors used\n")
        output_text.insert(tk.END, "=" * 60 + "\n")
        sh_w = max(len(s) for s in _sheet_detector_map) + 2
        has_sn = any(v['sn'] for v in _sheet_detector_map.values())
        for sheet, info in _sheet_detector_map.items():
            sn_str = f"SN {info['sn']:<8}" if has_sn else ""
            raw_dets = ", ".join(sorted(info['detectors'])) or "—"
            canon = info.get('canonical')
            canon_str = f"  → {canon}" if canon and canon != raw_dets else ""
            output_text.insert(tk.END, f"  {sheet:<{sh_w}}  {sn_str}{raw_dets}{canon_str}\n")

    # ── plot: one subplot per group, per-group y-axis scaling ──
    import math
    n = len(group_keys)
    ncols = min(n, 3)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.8 * nrows), squeeze=False)
    fig.suptitle(f"{metric_name}{axis_str} at FS {fs_ref:g} cm, depth {depth_ref:g} cm", fontsize=12)

    rng = np.random.default_rng(0)
    for i, k in enumerate(group_keys):
        r, col_ = divmod(i, ncols)
        ax = axes[r][col_]
        e, s = k
        if use_energy_fs_pos_spec:
            spec_val = (energy_fs_pos_data.get(e, {}).get(fs_ref, {}).get(oar_p)
                        if oar_p is not None else None)
        elif use_energy_fs_spec:
            spec_val = energy_fs_data.get(e, {}).get(fs_ref)
        elif use_per_energy_spec:
            spec_val = energy_specs.get(e)
        elif spec_fn is not None:
            spec_val = spec_fn(s, depth_ref, fs_ref)
        sheet_items = [(sh, v) for sh, v in results[k].items() if not np.isnan(v)]
        vals = [v for _, v in sheet_items]

        ax.boxplot(vals, positions=[0], widths=0.45, patch_artist=True,
                   boxprops=dict(facecolor='lightgray'))
        xs = (rng.random(len(sheet_items)) - 0.5) * 0.15
        label_all = show_labels_var.get()
        do_color_tps = color_tps_var.get()
        for (sh, v), x in zip(sheet_items, xs):
            if highlight and (sh == highlight or sh.startswith(highlight + "_")):
                ax.plot(x, v, 'o', color='red', ms=11, zorder=5)
                ax.annotate(sh, (x, v), xytext=(8, 0), textcoords='offset points',
                            color='red', fontweight='bold', fontsize=9, va='center')
            else:
                is_tps = do_color_tps and "TPS" in sh.upper()
                pt_color = 'darkorange' if is_tps else 'black'
                ax.plot(x, v, 'o', color=pt_color, ms=6)
                if label_all:
                    ax.annotate(sh, (x, v), xytext=(8, 0), textcoords='offset points',
                                color=pt_color, fontsize=7, va='center')

        ax.set_xticks([0])
        ax.set_xticklabels([f"{e or '?'}  {s or '?'} SSD"])
        ax.set_ylabel(f"{metric_name} [{unit}]" if unit else metric_name)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        # Pad y-range a bit around the data so the box isn't hugging the edges.
        # Floor scales with the value magnitude (not an absolute 0.05) so small
        # dimensionless metrics like TMR 20/10 (~0.7) aren't dwarfed by the pad.
        # Include spec/tol values in the range so lines are always visible.
        tps_val = _tps_values.get(metric_name, {}).get(e) if show_tps_var.get() else None

        cc13_pion = (_pion_corrections.get("CC 13", {}).get(e)
                     if apply_pion_var.get() and show_pion_spec_var.get() else None)
        spec_corr = spec_val * cc13_pion if (spec_val is not None and cc13_pion is not None) else None

        y_range = list(vals)
        if show_spec_var.get() and spec_val is not None:
            y_range.append(spec_val)
            if tol_val is not None:
                y_range.extend([spec_val + tol_val, spec_val - tol_val])
        if spec_corr is not None:
            y_range.append(spec_corr)
            if tol_val is not None:
                y_range.extend([spec_corr + tol_val, spec_corr - tol_val])
        if tps_val is not None:
            y_range.append(tps_val)
            if tps_tol_val is not None:
                y_range.extend([tps_val + tps_tol_val, tps_val - tps_tol_val])

        vmin, vmax = min(y_range), max(y_range)
        pad = max((vmax - vmin) * 0.3, abs(vmax) * 0.01)
        ax.set_ylim(vmin - pad, vmax + pad)

        if show_spec_var.get() and spec_val is not None:
            _line_lbl = 'Tol. limit' if comparison == "max" else 'Spec'
            ax.axhline(spec_val, color='green', linewidth=1.2, alpha=0.8, label=_line_lbl)
            if tol_val is not None:
                ax.axhline(spec_val + tol_val, color='red', linestyle='--', linewidth=1.0, alpha=0.8, label='Tol')
                ax.axhline(spec_val - tol_val, color='red', linestyle='--', linewidth=1.0, alpha=0.8)
        if spec_corr is not None:
            ax.axhline(spec_corr, color='green', linestyle='--', linewidth=1.2, alpha=0.8,
                       label='Spec (CC13 pion)')
            if tol_val is not None:
                ax.axhline(spec_corr + tol_val, color='red', linestyle=':', linewidth=1.0, alpha=0.8)
                ax.axhline(spec_corr - tol_val, color='red', linestyle=':', linewidth=1.0, alpha=0.8)
        if tps_val is not None:
            ax.axhline(tps_val, color='steelblue', linestyle='--', linewidth=1.4, alpha=0.9, label='TPS')
            if tps_tol_val is not None:
                ax.axhline(tps_val + tps_tol_val, color='steelblue', linestyle=':', linewidth=1.0, alpha=0.7,
                           label=f'TPS ±{tps_tol_val:g}')
                ax.axhline(tps_val - tps_tol_val, color='steelblue', linestyle=':', linewidth=1.0, alpha=0.7)

    # Hide any unused axes (when n isn't a multiple of ncols)
    for j in range(n, nrows * ncols):
        r, col_ = divmod(j, ncols)
        axes[r][col_].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


def _open_detector_alias_dialog():
    dlg = tk.Toplevel(root)
    dlg.title("Detector Aliases")
    dlg.resizable(True, True)
    dlg.grab_set()

    ttk.Label(dlg, text="Raw name (in data)", width=36, anchor="w").grid(
        row=0, column=0, padx=8, pady=4, sticky="w")
    ttk.Label(dlg, text="Canonical name", width=20, anchor="w").grid(
        row=0, column=1, padx=8, pady=4, sticky="w")

    # Only actual raw names seen in data; canonical pre-filled via full resolution
    all_raw = sorted(_known_raw_detectors)
    entries = {}
    for i, raw in enumerate(all_raw):
        ttk.Label(dlg, text=raw, width=36, anchor="w").grid(
            row=i + 1, column=0, padx=8, pady=2, sticky="w")
        resolved = _detector_aliases.get(raw) or _canonical_detector(raw)
        var = tk.StringVar(value=resolved if resolved != raw else "")
        ttk.Entry(dlg, textvariable=var, width=22).grid(
            row=i + 1, column=1, padx=8, pady=2)
        entries[raw] = var

    # blank rows for adding new entries
    extra_vars = []
    base = len(all_raw) + 1
    for j in range(3):
        raw_var = tk.StringVar()
        can_var = tk.StringVar()
        ttk.Entry(dlg, textvariable=raw_var, width=36).grid(row=base + j, column=0, padx=8, pady=2)
        ttk.Entry(dlg, textvariable=can_var, width=22).grid(row=base + j, column=1, padx=8, pady=2)
        extra_vars.append((raw_var, can_var))

    def _save():
        _detector_aliases.clear()
        for raw, var in entries.items():
            txt = var.get().strip()
            if txt:
                _detector_aliases[raw] = txt
        for rv, cv in extra_vars:
            r, c = rv.get().strip(), cv.get().strip()
            if r and c:
                _detector_aliases[r] = c
        dlg.destroy()

    n = base + 3
    ttk.Button(dlg, text="OK",     command=_save       ).grid(row=n, column=0, pady=8)
    ttk.Button(dlg, text="Cancel", command=dlg.destroy ).grid(row=n, column=1, pady=8)
    dlg.bind("<Return>", lambda *_: _save())

detector_alias_btn.configure(command=_open_detector_alias_dialog)


def _open_pion_dialog():
    # Canonical detector names = alias values, plus any already in the pion table
    canonicals = sorted(set(_detector_aliases.values()) | set(_pion_corrections.keys()))
    if not canonicals:
        messagebox.showinfo("Pion Corrections",
                            "No canonical detector names yet.\n"
                            "Set Detector Aliases first, then run once to discover names.")
        return

    dlg = tk.Toplevel(root)
    dlg.title("Pion Corrections (multiplicative, per detector/energy)")
    dlg.resizable(True, True)
    dlg.grab_set()

    ttk.Label(dlg, text="Detector \\ Energy", width=18).grid(
        row=0, column=0, padx=6, pady=4, sticky="w")
    for j, e in enumerate(_ENERGIES):
        ttk.Label(dlg, text=e, width=9).grid(row=0, column=j + 1, padx=4, pady=4)

    entries = {}
    for i, det in enumerate(canonicals):
        ttk.Label(dlg, text=det, width=18, anchor="w").grid(
            row=i + 1, column=0, padx=6, pady=3, sticky="w")
        entries[det] = {}
        for j, e in enumerate(_ENERGIES):
            cur = _pion_corrections.get(det, {}).get(e, "")
            var = tk.StringVar(value="" if cur == "" else str(cur))
            ttk.Entry(dlg, textvariable=var, width=9).grid(
                row=i + 1, column=j + 1, padx=4, pady=3)
            entries[det][e] = var

    ttk.Label(dlg, text="Leave blank = 1.0 (no correction)",
              foreground="gray").grid(row=len(canonicals) + 1, column=0,
                                     columnspan=len(_ENERGIES) + 1, padx=6, pady=4)

    def _save():
        _pion_corrections.clear()
        for det, evars in entries.items():
            for e, var in evars.items():
                txt = var.get().strip()
                if txt:
                    try:
                        f = float(txt)
                        _pion_corrections.setdefault(det, {})[e] = f
                    except ValueError:
                        pass
        dlg.destroy()

    def _clear():
        for evars in entries.values():
            for var in evars.values():
                var.set("")

    n = len(canonicals) + 2
    btn_row = ttk.Frame(dlg)
    btn_row.grid(row=n, column=0, columnspan=len(_ENERGIES) + 1, pady=8)
    ttk.Button(btn_row, text="OK",     command=_save       ).pack(side="left", padx=6)
    ttk.Button(btn_row, text="Clear",  command=_clear      ).pack(side="left", padx=6)
    ttk.Button(btn_row, text="Cancel", command=dlg.destroy ).pack(side="left", padx=6)

pion_btn.configure(command=_open_pion_dialog)


def load_data():
    if not loaded_files:
        messagebox.showerror("Error", "Load at least one file first.")
        return
    sel_idx = sheet_listbox.curselection()
    if not sel_idx:
        messagebox.showerror("Error", "Select at least one sheet.")
        return
    sheets = [sheet_listbox.get(i) for i in sel_idx]

    metric_name = metric_var.get()
    file_kw = METRICS[metric_name].get("file_keyword")

    _df_cache.clear()
    n_loaded = 0
    load_status_lbl.config(text="Loading…", foreground="darkorange")
    root.update_idletasks()
    for entry in loaded_files:
        if file_kw and file_kw.lower() not in os.path.basename(entry['path']).lower():
            continue
        wanted = [s for s in sheets if s in entry['sheets']]
        if not wanted:
            continue
        load_status_lbl.config(text=f"Loading {os.path.basename(entry['path'])}…")
        root.update_idletasks()
        try:
            dfs = pd.read_excel(entry['path'], sheet_name=wanted)
            _df_cache[entry['path']] = dfs
            n_loaded += len(dfs)
            try:
                fs_ref = float(fs_var.get())
                depth_ref = float(depth_var.get())
            except ValueError:
                fs_ref = None
                depth_ref = None
            for df in dfs.values():
                sub = df
                if fs_ref is not None and 'FS' in df.columns:
                    try:
                        sub = df[df['FS'].astype(float) == fs_ref]
                    except Exception:
                        sub = df
                if file_kw == 'PDD' and 'Axis' in sub.columns:
                    sub = sub[sub['Axis'] == 'Z']
                elif file_kw == 'Profile' and depth_ref is not None and 'Depth' in sub.columns:
                    avail = sub['Depth'].dropna().to_numpy(dtype=float)
                    if len(avail) > 0:
                        nearest = avail[np.argmin(np.abs(avail - depth_ref))]
                        sub = sub[sub['Depth'] == nearest]
                for col in ("Detector", "Chamber"):
                    if col in sub.columns:
                        vals = [str(v) for v in sub[col].dropna().unique()
                                if str(v).strip() not in ("", "nan")]
                        if vals:
                            _known_raw_detectors.update(vals)
                            break
        except Exception as e:
            messagebox.showerror("Error", f"Could not read {os.path.basename(entry['path'])}:\n{e}")
            load_status_lbl.config(text="Load failed", foreground="red")
            return
    load_status_lbl.config(text=f"{n_loaded} sheet(s) loaded", foreground="green")


load_btn_frame = ttk.Frame(main)
load_btn_frame.grid(row=13, column=0, columnspan=3, pady=(10, 0))
ttk.Button(load_btn_frame, text="Load Data", command=load_data).pack(side="left", padx=8)
ttk.Button(load_btn_frame, text="Run",       command=run_compare).pack(side="left", padx=8)
load_status_lbl = ttk.Label(load_btn_frame, text="")
load_status_lbl.pack(side="left", padx=8)

root.mainloop()
