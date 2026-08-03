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


def metric_dmax_depth(df, fs, depth_cm):
    """Depth of dmax [cm] via pchip interpolation on a fine grid."""
    y, d = _z_axis_sorted(df, fs)
    if y is None:
        return float('nan')
    try:
        pchip = interp.pchip(y, d)
    except Exception:
        return float('nan')
    # Search on a fine grid over the shallow half of the scan where dmax should live.
    y_lo = float(y.min())
    y_hi = min(float(y.max()), y_lo + 5.0)   # dmax is always within top ~5 cm of clinical beams
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
                         "spec": {e: None for e in _ENERGIES}, "tol": 0.5},
    "Depth of dmax":    {"fn": metric_dmax_depth,   "file_keyword": "PDD",     "unit": "cm",
                         "spec": {e: None for e in _ENERGIES}, "tol": 0.5},
    "TMR 20/10":        {"fn": metric_tmr_20_10,    "file_keyword": "PDD",     "unit": "",
                         "spec": None, "tol": None},
    "Profile Flatness":   {"fn": metric_flatness,    "file_keyword": "Profile", "unit": "%",
                           "spec": 0.0, "tol": 3.0,  "needs_axis": True},
    "Profile Field Size": {"fn": metric_field_size,  "file_keyword": "Profile", "unit": "cm",
                           "spec": None, "tol": 0.2,  "needs_axis": True},
    "Profile Symmetry":   {"fn": metric_symmetry,   "file_keyword": "Profile", "unit": "%",
                           "spec": 0.0,  "tol": 3.0,  "needs_axis": True},
    "Profile Penumbra":   {"fn": metric_penumbra,  "file_keyword": "Profile", "unit": "cm",
                           "spec": None, "tol": None, "needs_axis": True, "needs_side": True},
    "Profile OAR":        {"fn": metric_oar,       "file_keyword": "Profile", "unit": "%",
                           "spec": None, "tol": None, "needs_axis": True, "needs_side": True,
                           "needs_oar_pos": True},
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
_user_specs  = {}   # {metric_name: {energy: float or None}}


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Beam Specifier Compare v0.2")

fs_var          = tk.StringVar(master=root, value="10.5")
depth_var       = tk.StringVar(master=root, value="10")
metric_var      = tk.StringVar(master=root, value="PDD at depth")
highlight_var   = tk.StringVar(master=root)
show_labels_var = tk.BooleanVar(master=root, value=False)
axis_var        = tk.StringVar(master=root, value="X")
side_var        = tk.StringVar(master=root, value="Both")
oar_pos_var     = tk.StringVar(master=root, value="2.0")
spec_var        = tk.StringVar(master=root, value="")
tol_var         = tk.StringVar(master=root, value="")

main = ttk.Frame(root, padding=10)
main.grid(row=0, column=0, sticky="nsew")


def _refresh_sheets():
    """Repopulate sheet listbox with union of sheet names across loaded files."""
    sheet_listbox.delete(0, tk.END)
    all_sheets = []
    seen = set()
    for entry in loaded_files:
        for s in entry['sheets']:
            if s not in seen:
                all_sheets.append(s); seen.add(s)
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

# Row 7: depth entry + Tolerance
ttk.Label(main, text="Reference depth [cm]:").grid(row=7, column=0, sticky="w")
ttk.Entry(main, textvariable=depth_var, width=8).grid(row=7, column=1, sticky="w", padx=5)
tol_frame = ttk.Frame(main)
tol_frame.grid(row=7, column=2, sticky="w", padx=5)
ttk.Label(tol_frame, text="Tolerance ±:").pack(side='left')
ttk.Entry(tol_frame, textvariable=tol_var, width=8).pack(side='left', padx=4)

# Row 8: profile axis selector (shown only for axis-aware metrics)
axis_label = ttk.Label(main, text="Profile axis:")
axis_label.grid(row=8, column=0, sticky="w")
axis_frame = ttk.Frame(main)
axis_frame.grid(row=8, column=1, columnspan=2, sticky="w")
for _txt, _val in (("X", "X"), ("Y", "Y"), ("Both (worst)", "Both")):
    ttk.Radiobutton(axis_frame, text=_txt, variable=axis_var, value=_val).pack(side='left', padx=6)

# Row 9: profile side selector (shown only for metrics with needs_side)
side_label = ttk.Label(main, text="Side:")
side_label.grid(row=9, column=0, sticky="w")
side_frame = ttk.Frame(main)
side_frame.grid(row=9, column=1, columnspan=2, sticky="w")
for _txt, _val in (("Left", "Left"), ("Right", "Right"), ("Both (avg)", "Both")):
    ttk.Radiobutton(side_frame, text=_txt, variable=side_var, value=_val).pack(side='left', padx=6)
side_label.grid_remove()
side_frame.grid_remove()


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


def _on_metric_change(*_):
    info = METRICS.get(metric_var.get(), {})
    if info.get("needs_axis"):
        axis_label.grid()
        axis_frame.grid()
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
    if isinstance(spec_info, dict):
        spec_single_frame.grid_remove()
        spec_multi_frame.grid()
        spec_var.set("")
    else:
        spec_multi_frame.grid_remove()
        spec_single_frame.grid()
        spec_var.set("" if spec_info is None else str(spec_info))
    tol_var.set("" if info.get("tol") is None else str(info["tol"]))


metric_combo.bind("<<ComboboxSelected>>", _on_metric_change)
_on_metric_change()   # initialize visibility

# Row 10: off-axis position entry (shown only for OAR metric)
oar_pos_frame = ttk.Frame(main)
oar_pos_frame.grid(row=10, column=0, columnspan=2, sticky="w")
ttk.Label(oar_pos_frame, text="Off-axis position [cm]:").pack(side='left')
ttk.Entry(oar_pos_frame, textvariable=oar_pos_var, width=8).pack(side='left', padx=5)
oar_pos_frame.grid_remove()

# Row 11: highlight + show all labels
ttk.Label(main, text="Highlight sheet:").grid(row=11, column=0, sticky="w")
highlight_combo = ttk.Combobox(main, textvariable=highlight_var, values=[''],
                               state="readonly", width=20)
highlight_combo.grid(row=11, column=1, sticky="w", padx=5)
ttk.Checkbutton(main, text="Label all points", variable=show_labels_var).grid(
    row=11, column=2, sticky="w", padx=5)

# Row 12: Run button (above output so it's always visible)
# Row 12 run button is bound after run_compare is defined, below.

# Row 13: output
output_frame = ttk.Frame(main)
output_frame.grid(row=13, column=0, columnspan=3, sticky="w")
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

    if metric_info.get("needs_axis"):
        _axis = axis_var.get()
        if metric_info.get("needs_side"):
            _side = side_var.get()
            if metric_info.get("needs_oar_pos"):
                try:
                    _oar = float(oar_pos_var.get())
                except ValueError:
                    _oar = 2.0
                metric_fn = lambda df, fs, d, __ax=_axis, __sd=_side, __p=_oar: metric_info["fn"](df, fs, d, __ax, __sd, __p)
                axis_str = f", axis={_axis}, side={_side}, pos=±{_oar}cm"
            else:
                metric_fn = lambda df, fs, d, __ax=_axis, __sd=_side: metric_info["fn"](df, fs, d, __ax, __sd)
                axis_str = f", axis={_axis}, side={_side}"
        else:
            metric_fn = lambda df, fs, d, __ax=_axis: metric_info["fn"](df, fs, d, __ax)
            axis_str = f", axis={_axis}"
    else:
        metric_fn = metric_info["fn"]
        axis_str = ""

    # results[(energy, ssd)] = { sheet: value }
    results = {}

    spec_info = metric_info.get("spec")
    use_per_energy_spec = isinstance(spec_info, dict)
    if use_per_energy_spec:
        energy_specs = _user_specs.get(metric_name, spec_info)
        spec_val = None   # resolved per group below
    else:
        try:
            spec_val = float(spec_var.get()) if spec_var.get().strip() else None
        except ValueError:
            spec_val = None
    try:
        tol_val = float(tol_var.get()) if tol_var.get().strip() else None
    except ValueError:
        tol_val = None

    output_text.delete("1.0", tk.END)
    output_text.insert(tk.END, f"{metric_name}  (FS={fs_ref} cm, depth={depth_ref} cm{axis_str})\n")
    output_text.insert(tk.END, "=" * 70 + "\n")

    for entry in loaded_files:
        energy = entry['energy']; ssd = entry['ssd']
        if use_per_energy_spec:
            spec_val = energy_specs.get(energy)
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

        # Read only the needed sheets in ONE open (much faster than per-sheet read_excel)
        wanted = [s for s in sheets if s in entry['sheets']]
        if not wanted:
            continue
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
            val = metric_fn(df, fs_ref, depth_ref)
            results[group_key][s] = val
            flag = ""
            if spec_val is not None and tol_val is not None and not np.isnan(val):
                if abs(val - spec_val) > tol_val:
                    flag = "  ***"
            output_text.insert(tk.END, f"  {s:<15}  {val:.3f}{flag}\n")

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
    hdr = (f"{'Group':<18} {'n':>3} {'mean':>8} {'median':>8} {'std':>7} "
           f"{'P25':>8} {'P75':>8} {'min':>8} {'max':>8} {'range':>7} {'CoV%':>6}")
    output_text.insert(tk.END, hdr + "\n")
    output_text.insert(tk.END, "-" * 104 + "\n")
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
        output_text.insert(tk.END,
            f"{label:<18} {vals.size:>3d} {mean:>8.3f} {median:>8.3f} {std:>7.3f} "
            f"{p25:>8.3f} {p75:>8.3f} {vmin:>8.3f} {vmax:>8.3f} {rng_v:>7.3f} {cov:>6.2f}\n"
        )

    # Highlighted sheet position within each group
    if highlight:
        output_text.insert(tk.END, "\n" + "-" * 104 + "\n")
        output_text.insert(tk.END, f"Position of '{highlight}' within each group:\n")
        for k in group_keys:
            vals_all = results[k]
            if highlight not in vals_all or np.isnan(vals_all[highlight]):
                continue
            v = vals_all[highlight]
            others = np.array([x for s_, x in vals_all.items()
                               if s_ != highlight and not np.isnan(x)], dtype=float)
            if others.size == 0:
                continue
            mean = float(others.mean())
            std  = float(others.std(ddof=1)) if others.size > 1 else 0.0
            z    = (v - mean) / std if std > 0 else float('nan')
            delta = v - mean
            e, s = k
            label = f"{e or '?'} {s or '?'} SSD"
            output_text.insert(tk.END,
                f"  {label:<18} value={v:.3f}  Δ(mean others)={delta:+.3f}  "
                f"z={z:+.2f}σ\n"
            )

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
        if use_per_energy_spec:
            spec_val = energy_specs.get(e)
        sheet_items = [(sh, v) for sh, v in results[k].items() if not np.isnan(v)]
        vals = [v for _, v in sheet_items]

        ax.boxplot(vals, positions=[0], widths=0.45, patch_artist=True,
                   boxprops=dict(facecolor='lightgray'))
        xs = (rng.random(len(sheet_items)) - 0.5) * 0.15
        label_all = show_labels_var.get()
        for (sh, v), x in zip(sheet_items, xs):
            if highlight and sh == highlight:
                ax.plot(x, v, 'o', color='red', ms=11, zorder=5)
                ax.annotate(sh, (x, v), xytext=(8, 0), textcoords='offset points',
                            color='red', fontweight='bold', fontsize=9, va='center')
            else:
                ax.plot(x, v, 'o', color='black', ms=6)
                if label_all:
                    ax.annotate(sh, (x, v), xytext=(8, 0), textcoords='offset points',
                                color='gray', fontsize=7, va='center')

        ax.set_xticks([0])
        ax.set_xticklabels([f"{e or '?'}  {s or '?'} SSD"])
        ax.set_ylabel(f"{metric_name} [{unit}]" if unit else metric_name)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        # Pad y-range a bit around the data so the box isn't hugging the edges.
        # Floor scales with the value magnitude (not an absolute 0.05) so small
        # dimensionless metrics like TMR 20/10 (~0.7) aren't dwarfed by the pad.
        if spec_val is not None:
            ax.axhline(spec_val, color='green', linewidth=1.2, alpha=0.8, label='Spec')
            if tol_val is not None:
                ax.axhline(spec_val + tol_val, color='red', linestyle='--', linewidth=1.0, alpha=0.8, label='Tol')
                ax.axhline(spec_val - tol_val, color='red', linestyle='--', linewidth=1.0, alpha=0.8)

        vmin, vmax = min(vals), max(vals)
        pad = max((vmax - vmin) * 0.3, abs(vmax) * 0.01)
        ax.set_ylim(vmin - pad, vmax + pad)

    # Hide any unused axes (when n isn't a multiple of ncols)
    for j in range(n, nrows * ncols):
        r, col_ = divmod(j, ncols)
        axes[r][col_].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


ttk.Button(main, text="Run", command=run_compare).grid(row=12, column=0, columnspan=3, pady=10)

root.mainloop()
