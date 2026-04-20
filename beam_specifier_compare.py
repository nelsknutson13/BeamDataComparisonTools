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


METRICS = {
    "PDD at depth": {"fn": metric_pdd_at_depth, "file_keyword": "PDD", "unit": "%"},
    "Depth of dmax": {"fn": metric_dmax_depth,  "file_keyword": "PDD", "unit": "cm"},
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


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Beam Specifier Compare v0.2")

fs_var        = tk.StringVar(master=root, value="10.5")
depth_var     = tk.StringVar(master=root, value="10")
metric_var    = tk.StringVar(master=root, value="PDD at depth")
highlight_var = tk.StringVar(master=root)
show_labels_var = tk.BooleanVar(master=root, value=False)

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

# Row 6: FS entry
ttk.Label(main, text="Reference field size [cm]:").grid(row=6, column=0, sticky="w")
ttk.Entry(main, textvariable=fs_var, width=8).grid(row=6, column=1, sticky="w", padx=5)

# Row 7: depth entry
ttk.Label(main, text="Reference depth [cm]:").grid(row=7, column=0, sticky="w")
ttk.Entry(main, textvariable=depth_var, width=8).grid(row=7, column=1, sticky="w", padx=5)

# Row 8: highlight + show all labels
ttk.Label(main, text="Highlight sheet:").grid(row=8, column=0, sticky="w")
highlight_combo = ttk.Combobox(main, textvariable=highlight_var, values=[''],
                               state="readonly", width=20)
highlight_combo.grid(row=8, column=1, sticky="w", padx=5)
ttk.Checkbutton(main, text="Label all points", variable=show_labels_var).grid(
    row=8, column=2, sticky="w", padx=5)

# Row 9: Run button (above output so it's always visible)
# Row 9 run button is bound after run_compare is defined, below.

# Row 10: output
output_frame = ttk.Frame(main)
output_frame.grid(row=10, column=0, columnspan=3, sticky="w")
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
    metric_fn   = metric_info["fn"]
    file_kw     = metric_info.get("file_keyword")
    unit        = metric_info.get("unit", "")
    highlight   = highlight_var.get().strip() or None

    # results[(energy, ssd)] = { sheet: value }
    results = {}

    output_text.delete("1.0", tk.END)
    output_text.insert(tk.END, f"{metric_name}  (FS={fs_ref} cm, depth={depth_ref} cm)\n")
    output_text.insert(tk.END, "=" * 70 + "\n")

    for entry in loaded_files:
        energy = entry['energy']; ssd = entry['ssd']
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
            output_text.insert(tk.END, f"  {s:<15}  {val:.3f}\n")

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
    fig.suptitle(f"{metric_name} at FS {fs_ref:g} cm, depth {depth_ref:g} cm", fontsize=12)

    rng = np.random.default_rng(0)
    for i, k in enumerate(group_keys):
        r, col_ = divmod(i, ncols)
        ax = axes[r][col_]
        e, s = k
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
        vmin, vmax = min(vals), max(vals)
        pad = max((vmax - vmin) * 0.3, 0.05)
        ax.set_ylim(vmin - pad, vmax + pad)

    # Hide any unused axes (when n isn't a multiple of ncols)
    for j in range(n, nrows * ncols):
        r, col_ = divmod(j, ncols)
        axes[r][col_].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


ttk.Button(main, text="Run", command=run_compare).grid(row=9, column=0, columnspan=3, pady=10)

root.mainloop()
