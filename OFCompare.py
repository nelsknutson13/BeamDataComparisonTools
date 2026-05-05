"""
Output Factor Compare
---------------------
Load a single long-format combined OF dataset and compare two groupings of
that data (e.g. SN 21 vs SN 19, or TPS vs measured). Built to consume the
file produced by OFDataConverter.

Long-format columns expected:
    Energy, SSD, Depth, SN, Detector, FS_X, FS_Y, Scp

Views:
    - 2D heatmap of Scp (or % diff between A and B)
    - Diagonal trace (X=Y, square fields) as Scp vs FS
    - Row trace (fixed Y, varying X)
    - Column trace (fixed X, varying Y)
"""

import os
import re
import tkinter as tk
from tkinter import filedialog, ttk, messagebox, simpledialog

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


LONG_REQUIRED = {'FS_X', 'FS_Y', 'Scp'}
FILTER_COLS   = ['Energy', 'SSD', 'Depth', 'SN', 'Detector']
DEST_SHEET_DEFAULT = 'OutputFactors'
ENERGY_TOKENS = ['6X', '10X', '15X', '6FFF', '8FFF', '10FFF']


def _is_numeric(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return False
    return f == f


def _wide_to_long(df_raw):
    """Convert a wide-format 2D Y×X sheet into a long DataFrame with FS_X, FS_Y, Scp."""
    arr = df_raw.values
    n_rows, n_cols = arr.shape
    x_header_row = None
    for r in range(min(n_rows, 6)):
        if n_cols < 2:
            continue
        numeric_count = sum(_is_numeric(v) for v in arr[r, 1:])
        if numeric_count >= max(3, int(0.6 * (n_cols - 1))):
            x_header_row = r
            break
    if x_header_row is None:
        raise RuntimeError("Could not find a numeric X-axis header row.")
    x_vals, x_cols = [], []
    for c in range(1, n_cols):
        v = arr[x_header_row, c]
        if _is_numeric(v):
            x_vals.append(float(v))
            x_cols.append(c)
    rows = []
    for r in range(x_header_row + 1, n_rows):
        y = arr[r, 0]
        if not _is_numeric(y):
            continue
        for c, x in zip(x_cols, x_vals):
            v = arr[r, c]
            if _is_numeric(v):
                rows.append({'FS_X': x, 'FS_Y': float(y), 'Scp': float(v)})
    if not rows:
        raise RuntimeError("No data rows below the X header.")
    return pd.DataFrame(rows)


def _guess_energy(text):
    if not text:
        return None
    s = str(text).upper().replace(' ', '').replace('_', '')
    for tok in ENERGY_TOKENS:
        if tok in s:
            return tok
    return None


def _guess_ssd(text):
    if not text:
        return None
    m = re.search(r'(\d{2,3})\s*SSD|SSD\s*(\d{2,3})', str(text).upper())
    if m:
        return int(m.group(1) or m.group(2))
    return None


def _prompt_wide_metadata(parent, default_ssd=None):
    """Pop up a single-shot dialog to gather SSD/Depth/SN/Detector for a wide-format file."""
    dlg = tk.Toplevel(parent)
    dlg.title("Metadata required for wide-format file")
    dlg.transient(parent); dlg.grab_set()
    ttk.Label(dlg, text="This file isn't in long format — please provide metadata.\n"
                        "Energy will be auto-detected per sheet.",
              padding=8, justify='left').grid(row=0, column=0, columnspan=2)
    vars_ = {
        'SSD':      tk.StringVar(master=dlg, value=str(default_ssd) if default_ssd is not None else '90'),
        'Depth':    tk.StringVar(master=dlg, value='10'),
        'SN':       tk.StringVar(master=dlg, value=''),
        'Detector': tk.StringVar(master=dlg, value=''),
    }
    for i, k in enumerate(['SSD', 'Depth', 'SN', 'Detector']):
        ttk.Label(dlg, text=f"{k}:").grid(row=i + 1, column=0, sticky="e", padx=8, pady=2)
        ttk.Entry(dlg, textvariable=vars_[k], width=20).grid(row=i + 1, column=1, sticky="w", padx=8, pady=2)
    result = {'value': None}
    def _ok():
        try:
            ssd   = float(vars_['SSD'].get())
            depth = float(vars_['Depth'].get())
        except ValueError:
            messagebox.showerror("Error", "SSD and Depth must be numeric.", parent=dlg)
            return
        sn  = vars_['SN'].get().strip()
        det = vars_['Detector'].get().strip()
        if not sn or not det:
            messagebox.showerror("Error", "SN and Detector are required.", parent=dlg)
            return
        result['value'] = {'SSD': ssd, 'Depth': depth, 'SN': sn, 'Detector': det}
        dlg.destroy()
    def _cancel():
        dlg.destroy()
    btns = ttk.Frame(dlg, padding=8)
    btns.grid(row=99, column=0, columnspan=2)
    ttk.Button(btns, text="OK",     command=_ok    ).pack(side='left', padx=4)
    ttk.Button(btns, text="Cancel", command=_cancel).pack(side='left', padx=4)
    dlg.wait_window()
    return result['value']


# ── Loader ───────────────────────────────────────────────────────────────────

def load_long(path, sheet=None, prompt_meta_cb=None):
    """Load a long-format combined OF file. If no sheet contains the required
    long-format columns, fall back to wide-format parsing of every sheet,
    using `prompt_meta_cb(default_ssd)` to gather SSD/Depth/SN/Detector once.
    Energy is auto-detected per sheet.
    """
    xl = pd.ExcelFile(path)
    if sheet is None:
        candidates = ([DEST_SHEET_DEFAULT] if DEST_SHEET_DEFAULT in xl.sheet_names else []) + \
                     [s for s in xl.sheet_names if s != DEST_SHEET_DEFAULT]
        for sh in candidates:
            df = pd.read_excel(xl, sheet_name=sh)
            if LONG_REQUIRED.issubset({str(c).strip() for c in df.columns}):
                return df, sh
        # No long-format sheet — try wide format on each sheet
        if prompt_meta_cb is None:
            raise RuntimeError("No long-format sheet found and no metadata prompt available.")
        meta = prompt_meta_cb(_guess_ssd(os.path.basename(path)))
        if meta is None:
            raise RuntimeError("Wide-format load cancelled.")
        long_parts = []
        for sh in xl.sheet_names:
            try:
                raw = pd.read_excel(xl, sheet_name=sh, header=None)
                wdf = _wide_to_long(raw)
            except Exception:
                continue
            wdf['Energy']   = _guess_energy(sh) or sh
            wdf['SSD']      = meta['SSD']
            wdf['Depth']    = meta['Depth']
            wdf['SN']       = meta['SN']
            wdf['Detector'] = meta['Detector']
            long_parts.append(wdf)
        if not long_parts:
            raise RuntimeError("No sheets could be parsed as wide-format OF data.")
        return pd.concat(long_parts, ignore_index=True), '<wide-format>'
    df = pd.read_excel(xl, sheet_name=sheet)
    if not LONG_REQUIRED.issubset({str(c).strip() for c in df.columns}):
        raise RuntimeError(f"Sheet '{sheet}' is not in long format.")
    return df, sheet


def _coerce_long(df):
    """Coerce types so dropdowns and filtering behave consistently."""
    out = df.copy()
    for c in ('FS_X', 'FS_Y', 'Scp', 'SSD', 'Depth'):
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors='coerce')
    for c in ('Energy', 'SN', 'Detector'):
        if c in out.columns:
            out[c] = out[c].astype(str).str.strip()
    out = out.dropna(subset=['FS_X', 'FS_Y', 'Scp']).reset_index(drop=True)
    return out


# ── Comparison helpers ───────────────────────────────────────────────────────

def _matrix(df):
    pivot = df.pivot_table(index='FS_Y', columns='FS_X', values='Scp', aggfunc='first')
    pivot = pivot.sort_index().sort_index(axis=1)
    return pivot.columns.values.astype(float), pivot.index.values.astype(float), pivot.values


def _common_axes(dfA, dfB):
    xA, yA, Ma = _matrix(dfA)
    xB, yB, Mb = _matrix(dfB)
    x_common = np.array(sorted(set(xA).intersection(xB)))
    y_common = np.array(sorted(set(yA).intersection(yB)))
    if len(x_common) == 0 or len(y_common) == 0:
        return None, None, None, None
    def _select(x_all, y_all, M):
        ix = np.array([list(x_all).index(x) for x in x_common])
        iy = np.array([list(y_all).index(y) for y in y_common])
        return M[np.ix_(iy, ix)]
    return x_common, y_common, _select(xA, yA, Ma), _select(xB, yB, Mb)


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Output Factor Compare v0.2")
# Size the window to ~85% of screen height so it always fits on the visible display.
_sw = root.winfo_screenwidth(); _sh = root.winfo_screenheight()
_w = min(1100, int(_sw * 0.9))
_h = int(_sh * 0.85)
root.geometry(f"{_w}x{_h}")
root.minsize(800, 500)

mode_var    = tk.StringVar(master=root, value='single')   # 'single' or 'dual'
file_a_var  = tk.StringVar(master=root)
file_b_var  = tk.StringVar(master=root)
status_var  = tk.StringVar(master=root, value="No file loaded.")

# Per-group filter vars
def _make_group():
    return {c: tk.StringVar(master=root, value='') for c in FILTER_COLS}

group_a = _make_group()
group_b = _make_group()

view_var    = tk.StringVar(master=root, value='diagonal')
compare_var = tk.StringVar(master=root, value='diff')
row_y_var   = tk.StringVar(master=root, value='10.5')
col_x_var   = tk.StringVar(master=root, value='10.5')

_df_a     = None     # long-format DataFrame backing Group A
_df_b     = None     # long-format DataFrame backing Group B (same as _df_a in single mode)
_combos = {'A': {}, 'B': {}}  # holds combobox widgets so we can repopulate


# Scrollable container so the GUI always fits regardless of screen size.
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

_outer = ttk.Frame(root)
_outer.grid(row=0, column=0, sticky="nsew")
_outer.grid_rowconfigure(0, weight=1)
_outer.grid_columnconfigure(0, weight=1)

_canvas_outer = tk.Canvas(_outer, highlightthickness=0)
_canvas_outer.grid(row=0, column=0, sticky="nsew")
_vsb_outer = ttk.Scrollbar(_outer, orient="vertical",   command=_canvas_outer.yview)
_hsb_outer = ttk.Scrollbar(_outer, orient="horizontal", command=_canvas_outer.xview)
_vsb_outer.grid(row=0, column=1, sticky="ns")
_hsb_outer.grid(row=1, column=0, sticky="ew")
_canvas_outer.configure(yscrollcommand=_vsb_outer.set, xscrollcommand=_hsb_outer.set)

main = ttk.Frame(_canvas_outer, padding=8)
_main_window_id = _canvas_outer.create_window((0, 0), window=main, anchor="nw")

def _on_main_resize(_event=None):
    _canvas_outer.configure(scrollregion=_canvas_outer.bbox("all"))
def _on_outer_resize(event):
    # Make the inner frame width track the canvas width so horizontal scroll only kicks in if needed.
    _canvas_outer.itemconfigure(_main_window_id, width=event.width)
main.bind("<Configure>", _on_main_resize)
_canvas_outer.bind("<Configure>", _on_outer_resize)

# Mouse-wheel scrolling on Windows
def _on_mousewheel(event):
    _canvas_outer.yview_scroll(int(-event.delta / 120), "units")
_canvas_outer.bind_all("<MouseWheel>", _on_mousewheel)


def _fmt_val(v):
    if isinstance(v, float) and v.is_integer():
        return f"{int(v)}"
    return str(v)


ALL_TOKEN = 'All'   # special value to mean "iterate over every unique value"
_ALL_COLS = {'SN', 'Detector', 'Energy'}


def _populate_dropdowns_for(grp_letter, df):
    """Fill the comboboxes for one group from a DataFrame's unique values.
    Auto-selects when only one value exists; SN and Detector get an "All" option
    when more than one value is present.
    """
    if df is None:
        return
    for col in FILTER_COLS:
        vals = sorted(df[col].dropna().unique(), key=lambda v: (str(type(v)), v))
        opts = [_fmt_val(v) for v in vals]
        if col in _ALL_COLS and len(opts) > 1:
            opts = [ALL_TOKEN] + opts
        cb = _combos[grp_letter].get(col)
        if cb is not None:
            cb['values'] = opts
            if col not in _ALL_COLS and len(opts) == 1:
                cb.set(opts[0])


def _refilter_options(grp_letter, group_vars, df_source):
    """Cascading filter: each combobox's options reflect what's possible given the OTHER selections.
    Bidirectional — picking detector narrows SN list, picking SN narrows detector list, etc.
    'All' on SN/Detector counts as no constraint for the cascade.
    """
    if df_source is None:
        return
    for col in FILTER_COLS:
        sub = df_source
        for other in FILTER_COLS:
            if other == col:
                continue
            v = group_vars[other].get().strip()
            if not v or v == ALL_TOKEN:
                continue
            if other in ('SSD', 'Depth'):
                try:
                    target = float(v)
                except ValueError:
                    continue
                sub = sub[np.isclose(sub[other].astype(float), target)]
            else:
                sub = sub[sub[other].astype(str) == v]
        opts = sorted(sub[col].dropna().unique(), key=lambda v: (str(type(v)), v))
        opts_str = [_fmt_val(o) for o in opts]
        if col in _ALL_COLS and len(opts_str) > 1:
            opts_str = [ALL_TOKEN] + opts_str
        cb = _combos[grp_letter].get(col)
        if cb is None:
            continue
        cb['values'] = opts_str
        cur = cb.get().strip()
        if cur and cur not in opts_str:
            cb.set('')


def _choose_file(var):
    p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls")])
    if p:
        var.set(p)


def _load_path(p):
    df, sh = load_long(p, prompt_meta_cb=lambda default_ssd: _prompt_wide_metadata(root, default_ssd))
    return _coerce_long(df), sh


def _clear_group_filters(group_vars):
    for v in group_vars.values():
        v.set('')


def do_load():
    """Load file(s) based on current mode, populate group dropdowns."""
    global _df_a, _df_b
    pa = file_a_var.get().strip()
    if not pa:
        messagebox.showerror("Error", "Pick file A first.")
        return
    try:
        _df_a, sha = _load_path(pa)
    except Exception as e:
        messagebox.showerror("Error loading A", str(e))
        return

    if mode_var.get() == 'dual':
        pb = file_b_var.get().strip()
        if not pb:
            messagebox.showerror("Error", "Pick file B (or switch to Single).")
            return
        try:
            _df_b, shb = _load_path(pb)
        except Exception as e:
            messagebox.showerror("Error loading B", str(e))
            return
        status_var.set(f"A: {len(_df_a)} rows from '{sha}'  |  B: {len(_df_b)} rows from '{shb}'.")
    else:
        _df_b = _df_a
        status_var.set(f"Loaded {len(_df_a)} rows from sheet '{sha}'.")

    # Clear stale selections so dropdowns repopulate with first available value
    _clear_group_filters(group_a)
    _clear_group_filters(group_b)
    _populate_dropdowns_for('A', _df_a)
    _populate_dropdowns_for('B', _df_b)


def _filter_group(df_source, group_vars):
    """Apply concrete filter values. 'All' on SN/Detector is treated as no constraint
    (the plot loop handles iteration over those columns)."""
    if df_source is None:
        return None
    df = df_source.copy()
    for col in FILTER_COLS:
        v = group_vars[col].get().strip()
        if not v or v == ALL_TOKEN:
            continue
        if col in ('SSD', 'Depth'):
            try:
                target = float(v)
            except ValueError:
                continue
            df = df[np.isclose(df[col], target)]
        else:
            df = df[df[col].astype(str) == v]
    return df


def _group_active(group_vars):
    """Group is 'active' if any filter is non-empty."""
    return any(v.get().strip() for v in group_vars.values())


def _clear_group(grp_letter):
    group_vars = group_a if grp_letter == 'A' else group_b
    df_src     = _df_a    if grp_letter == 'A' else _df_b
    _populate_dropdowns_for(grp_letter, df_src)   # restore full option lists
    _clear_group_filters(group_vars)              # then force every value back to blank


def _clear_all_filters():
    _clear_group('A')
    _clear_group('B')


def _plot():
    if _df_a is None or _df_b is None:
        messagebox.showinfo("Load first", "Load a file first.")
        return

    active_a = _group_active(group_a)
    active_b = _group_active(group_b)
    if not active_a and not active_b:
        messagebox.showerror("No selection",
                             "Set at least one filter on Group A or Group B.")
        return

    dfA = _filter_group(_df_a, group_a) if active_a else None
    dfB = _filter_group(_df_b, group_b) if active_b else None
    if active_a and (dfA is None or dfA.empty):
        messagebox.showerror("Empty A", "Group A has no rows after filtering.")
        return
    if active_b and (dfB is None or dfB.empty):
        messagebox.showerror("Empty B", "Group B has no rows after filtering.")
        return

    # short labels for legend / titles
    def _label(group_vars):
        parts = []
        for col in ('Energy', 'SN', 'Detector'):
            v = group_vars[col].get().strip()
            if v:
                parts.append(v)
        return ' / '.join(parts) if parts else 'group'

    label_a = _label(group_a)
    label_b = _label(group_b)

    view = view_var.get()
    mode = compare_var.get()

    fig = plt.figure(figsize=(10, 6))

    # ── single-group plotting: just show the data, no diff ──
    if active_a ^ active_b:
        group_vars = group_a if active_a else group_b
        df_only    = dfA if active_a else dfB
        base_label = label_a if active_a else label_b

        # Determine which All-columns to iterate over for this group
        en_all  = group_vars['Energy'].get().strip()   == ALL_TOKEN
        sn_all  = group_vars['SN'].get().strip()       == ALL_TOKEN
        det_all = group_vars['Detector'].get().strip() == ALL_TOKEN

        # Build list of (label, subset_df) traces
        traces = []
        if not en_all and not sn_all and not det_all:
            traces.append((base_label, df_only))
        else:
            keys = []
            if en_all:  keys.append('Energy')
            if sn_all:  keys.append('SN')
            if det_all: keys.append('Detector')
            combos = df_only[keys].drop_duplicates().sort_values(keys).values
            for combo in combos:
                sub = df_only.copy()
                tag_parts = []
                for k, v in zip(keys, combo):
                    sub = sub[sub[k].astype(str) == str(v)]
                    tag_parts.append(str(v))
                # Build label: replace 'All' tokens in base_label with concrete values
                base_parts = []
                for col in ('Energy', 'SN', 'Detector'):
                    v_set = group_vars[col].get().strip()
                    if v_set and v_set != ALL_TOKEN:
                        base_parts.append(v_set)
                lbl = ' / '.join(base_parts + tag_parts) if base_parts else ' / '.join(tag_parts)
                traces.append((lbl, sub))
            if not traces:
                messagebox.showerror("No combos", "No (SN, Detector) combinations found.")
                return

        if view == 'heatmap':
            if len(traces) > 1:
                messagebox.showinfo("Heatmap not multi-trace",
                                    "Heatmap supports a single trace; switch to Diagonal/Row/Col, "
                                    "or pick specific values instead of 'All'.")
                return
            x_vals, y_vals, M = _matrix(traces[0][1])
            ax = fig.add_subplot(111)
            im = ax.imshow(M, origin='lower', aspect='auto',
                           extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
            fig.colorbar(im, ax=ax, label='Scp')
            ax.set_title(traces[0][0])
            ax.set_xlabel('FS_X [cm]'); ax.set_ylabel('FS_Y [cm]')
        else:
            ax = fig.add_subplot(111)
            xlabel_set = None
            for tlabel, tdf in traces:
                x_vals, y_vals, M = _matrix(tdf)
                if view == 'diagonal':
                    common = sorted(set(x_vals).intersection(y_vals))
                    if not common:
                        continue
                    xs = np.array(common)
                    ys = np.array([M[list(y_vals).index(v), list(x_vals).index(v)] for v in common])
                    xlabel_set = 'FS  [cm]  (square fields, X=Y)'
                elif view == 'row':
                    try:
                        y_pick = float(row_y_var.get())
                    except ValueError:
                        messagebox.showerror("Bad Y", "Row Y= must be numeric."); return
                    iy = np.argmin(np.abs(y_vals - y_pick))
                    if not np.isclose(y_vals[iy], y_pick, atol=0.01):
                        continue
                    xs = x_vals; ys = M[iy, :]
                    xlabel_set = f'FS_X [cm]  (Y = {y_vals[iy]:g})'
                else:  # col
                    try:
                        x_pick = float(col_x_var.get())
                    except ValueError:
                        messagebox.showerror("Bad X", "Col X= must be numeric."); return
                    ix = np.argmin(np.abs(x_vals - x_pick))
                    if not np.isclose(x_vals[ix], x_pick, atol=0.01):
                        continue
                    xs = y_vals; ys = M[:, ix]
                    xlabel_set = f'FS_Y [cm]  (X = {x_vals[ix]:g})'
                ax.plot(xs, ys, 'o-', label=tlabel)
            ax.set_xlabel(xlabel_set or 'FS [cm]')
            ax.set_ylabel('Scp')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8, loc='best')
        fig.tight_layout()
        plt.show(block=False)
        plt.pause(0.001)
        return

    # ── both groups active: comparison view ──
    if view == 'heatmap':
        # Heatmap diff still needs overlapping axes.
        x_vals, y_vals, Ma, Mb = _common_axes(dfA, dfB)
        if x_vals is None:
            messagebox.showerror("No overlap",
                                 "Heatmap diff needs overlapping FS_X / FS_Y values.")
            return
        ax = fig.add_subplot(111)
        if mode == 'diff':
            pct = 100.0 * (Mb - Ma) / Ma
            vmax = max(abs(np.nanmin(pct)), abs(np.nanmax(pct)))
            im = ax.imshow(pct, origin='lower', aspect='auto',
                           cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                           extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
            fig.colorbar(im, ax=ax, label='% diff (B − A) / A')
            ax.set_title(f'% diff   {label_b}  vs  {label_a}')
        else:
            im = ax.imshow(Ma, origin='lower', aspect='auto',
                           extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
            fig.colorbar(im, ax=ax, label='Scp')
            ax.set_title(f'A: {label_a}')
        ax.set_xlabel('FS_X [cm]'); ax.set_ylabel('FS_Y [cm]')
        fig.tight_layout()
        plt.show(block=False)
        plt.pause(0.001)
        return

    # ── 1-D views (diagonal/row/col): build each group's data independently ──
    xA, yA, Ma = _matrix(dfA)
    xB, yB, Mb = _matrix(dfB)

    if view == 'diagonal':
        diag_a = sorted(set(xA).intersection(yA))
        diag_b = sorted(set(xB).intersection(yB))
        if not diag_a and not diag_b:
            messagebox.showerror("No diagonal", "No FS values shared between X and Y in either group.")
            return
        xs_a = np.array(diag_a, dtype=float)
        ys_a = np.array([Ma[list(yA).index(v), list(xA).index(v)] for v in diag_a])
        xs_b = np.array(diag_b, dtype=float)
        ys_b = np.array([Mb[list(yB).index(v), list(xB).index(v)] for v in diag_b])
        xlabel = 'FS  [cm]  (square fields, X=Y)'
    elif view == 'row':
        try:
            y_pick = float(row_y_var.get())
        except ValueError:
            messagebox.showerror("Bad Y", "Row Y= must be numeric."); return
        def _row_at(yvals, xvals, M, target):
            if len(yvals) == 0:
                return None, None
            iy = int(np.argmin(np.abs(yvals - target)))
            if not np.isclose(yvals[iy], target, atol=0.01):
                return None, None
            return xvals, M[iy, :]
        xs_a, ys_a = _row_at(yA, xA, Ma, y_pick)
        xs_b, ys_b = _row_at(yB, xB, Mb, y_pick)
        if xs_a is None and xs_b is None:
            messagebox.showerror("Y not found",
                                 f"Y={y_pick} not in either group's Y values."); return
        xlabel = f'FS_X [cm]  (Y = {y_pick:g})'
    else:  # col
        try:
            x_pick = float(col_x_var.get())
        except ValueError:
            messagebox.showerror("Bad X", "Col X= must be numeric."); return
        def _col_at(xvals, yvals, M, target):
            if len(xvals) == 0:
                return None, None
            ix = int(np.argmin(np.abs(xvals - target)))
            if not np.isclose(xvals[ix], target, atol=0.01):
                return None, None
            return yvals, M[:, ix]
        xs_a, ys_a = _col_at(xA, yA, Ma, x_pick)
        xs_b, ys_b = _col_at(xB, yB, Mb, x_pick)
        if xs_a is None and xs_b is None:
            messagebox.showerror("X not found",
                                 f"X={x_pick} not in either group's X values."); return
        xlabel = f'FS_Y [cm]  (X = {x_pick:g})'

    # Diff is computed only at FS values shared between groups.
    diff_xs = diff_ys = None
    if xs_a is not None and xs_b is not None and len(xs_a) and len(xs_b):
        common = sorted(set(np.asarray(xs_a, dtype=float)).intersection(
                        np.asarray(xs_b, dtype=float)))
        if common:
            map_a = dict(zip(np.asarray(xs_a, dtype=float), ys_a))
            map_b = dict(zip(np.asarray(xs_b, dtype=float), ys_b))
            diff_xs = np.array(common, dtype=float)
            ya_c = np.array([map_a[v] for v in common])
            yb_c = np.array([map_b[v] for v in common])
            with np.errstate(divide='ignore', invalid='ignore'):
                diff_ys = 100.0 * (yb_c - ya_c) / ya_c

    if mode == 'overlay' or diff_xs is None:
        ax = fig.add_subplot(111)
        if xs_a is not None and len(xs_a):
            ax.plot(xs_a, ys_a, 'o-', label=label_a)
        if xs_b is not None and len(xs_b):
            ax.plot(xs_b, ys_b, 's-', label=label_b)
        ax.set_xlabel(xlabel); ax.set_ylabel('Scp')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    else:
        ax1 = fig.add_subplot(2, 1, 1)
        if xs_a is not None and len(xs_a):
            ax1.plot(xs_a, ys_a, 'o-', label=label_a)
        if xs_b is not None and len(xs_b):
            ax1.plot(xs_b, ys_b, 's-', label=label_b)
        ax1.set_ylabel('Scp')
        ax1.grid(True, alpha=0.3); ax1.legend(fontsize=8)
        ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
        ax2.plot(diff_xs, diff_ys, 'o-', color='tab:red')
        ax2.axhline(0, color='k', lw=0.5)
        ax2.set_xlabel(xlabel); ax2.set_ylabel('% diff (B − A) / A')
        ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.show(block=False)
    plt.pause(0.001)


# ── layout ───────────────────────────────────────────────────────────────────

# Mode toggle
mode_frame = ttk.Frame(main)
mode_frame.grid(row=0, column=0, columnspan=4, sticky="w")
ttk.Label(mode_frame, text="Mode:").pack(side='left')
ttk.Radiobutton(mode_frame, text="Single file", variable=mode_var, value='single',
                command=lambda: _update_mode_visibility()).pack(side='left', padx=4)
ttk.Radiobutton(mode_frame, text="Dual file",   variable=mode_var, value='dual',
                command=lambda: _update_mode_visibility()).pack(side='left', padx=4)

# File A row
file_a_label = ttk.Label(main, text="File A:")
file_a_label.grid(row=1, column=0, sticky="w")
ttk.Entry(main, textvariable=file_a_var, width=70).grid(row=1, column=1, padx=4)
ttk.Button(main, text="Browse…", command=lambda: _choose_file(file_a_var)).grid(row=1, column=2)
ttk.Button(main, text="Load",    command=do_load).grid(row=1, column=3, padx=4)

# File B row (hidden in single mode)
file_b_label = ttk.Label(main, text="File B:")
file_b_entry = ttk.Entry(main, textvariable=file_b_var, width=70)
file_b_browse = ttk.Button(main, text="Browse…", command=lambda: _choose_file(file_b_var))

status_label = ttk.Label(main, textvariable=status_var, foreground='gray')
status_label.grid(row=3, column=0, columnspan=4, sticky="w")


def _update_mode_visibility():
    if mode_var.get() == 'dual':
        file_a_label.config(text="File A:")
        file_b_label.grid(row=2, column=0, sticky="w")
        file_b_entry.grid(row=2, column=1, padx=4)
        file_b_browse.grid(row=2, column=2)
    else:
        file_a_label.config(text="File:")
        file_b_label.grid_forget()
        file_b_entry.grid_forget()
        file_b_browse.grid_forget()
_update_mode_visibility()


def _make_filter_frame(label, group_vars, grp_letter, row_idx):
    frame = ttk.LabelFrame(main, text=f"Group {label}", padding=6)
    frame.grid(row=row_idx, column=0, columnspan=4, sticky="we", pady=(6, 0))
    cols = FILTER_COLS
    def _on_change(_evt=None):
        df_src = _df_a if grp_letter == 'A' else _df_b
        _refilter_options(grp_letter, group_vars, df_src)
    for i, c in enumerate(cols):
        ttk.Label(frame, text=f"{c}:").grid(row=0, column=2*i, sticky="e", padx=(4, 0))
        cb = ttk.Combobox(frame, textvariable=group_vars[c], width=12, state='normal')
        cb.grid(row=0, column=2*i+1, sticky="w", padx=(0, 4))
        cb.bind('<<ComboboxSelected>>', _on_change)
        cb.bind('<FocusOut>', _on_change)
        _combos[grp_letter][c] = cb


_make_filter_frame('A', group_a, 'A', 4)
_make_filter_frame('B', group_b, 'B', 5)

# View
view_frame = ttk.LabelFrame(main, text="View", padding=4)
view_frame.grid(row=6, column=0, columnspan=4, sticky="we", pady=(6, 0))
ttk.Radiobutton(view_frame, text="2D heatmap", variable=view_var, value='heatmap').pack(side='left', padx=4)
ttk.Radiobutton(view_frame, text="Diagonal (X=Y)", variable=view_var, value='diagonal').pack(side='left', padx=4)
ttk.Radiobutton(view_frame, text="Row (fixed Y)", variable=view_var, value='row').pack(side='left', padx=4)
ttk.Entry(view_frame, textvariable=row_y_var, width=6).pack(side='left')
ttk.Radiobutton(view_frame, text="Col (fixed X)", variable=view_var, value='col').pack(side='left', padx=(12, 4))
ttk.Entry(view_frame, textvariable=col_x_var, width=6).pack(side='left')

# Compare
cmp_frame = ttk.LabelFrame(main, text="Compare", padding=4)
cmp_frame.grid(row=7, column=0, columnspan=4, sticky="we", pady=(4, 0))
ttk.Radiobutton(cmp_frame, text="% difference (B vs A)", variable=compare_var, value='diff').pack(side='left', padx=4)
ttk.Radiobutton(cmp_frame, text="Overlay only",          variable=compare_var, value='overlay').pack(side='left', padx=4)
ttk.Button(cmp_frame, text="Plot", command=_plot).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear all", command=_clear_all_filters).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear B",   command=lambda: _clear_group('B')).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear A",   command=lambda: _clear_group('A')).pack(side='right', padx=4)

# Plots open in their own matplotlib window — see _plot() for figure creation.

main.grid_rowconfigure(8, weight=1)
main.grid_columnconfigure(1, weight=1)


root.mainloop()
