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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


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


main = ttk.Frame(root, padding=8)
main.grid(row=0, column=0, sticky="nsew")


def _fmt_val(v):
    if isinstance(v, float) and v.is_integer():
        return f"{int(v)}"
    return str(v)


def _populate_dropdowns_for(grp_letter, df):
    """Fill the comboboxes for one group from a DataFrame's unique values."""
    if df is None:
        return
    for col in FILTER_COLS:
        vals = sorted(df[col].dropna().unique(), key=lambda v: (str(type(v)), v))
        opts = [_fmt_val(v) for v in vals]
        cb = _combos[grp_letter].get(col)
        if cb is not None:
            cb['values'] = opts
            cur = cb.get()
            if not cur and opts:
                cb.set(opts[0])


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
    if df_source is None:
        return None
    df = df_source.copy()
    for col in FILTER_COLS:
        v = group_vars[col].get().strip()
        if not v:
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


def _plot():
    if _df_a is None or _df_b is None:
        messagebox.showinfo("Load first", "Load a file first.")
        return
    dfA = _filter_group(_df_a, group_a)
    dfB = _filter_group(_df_b, group_b)
    if dfA is None or dfA.empty:
        messagebox.showerror("Empty A", "Group A has no rows after filtering.")
        return
    if dfB is None or dfB.empty:
        messagebox.showerror("Empty B", "Group B has no rows after filtering.")
        return

    x_vals, y_vals, Ma, Mb = _common_axes(dfA, dfB)
    if x_vals is None:
        messagebox.showerror("No overlap", "Groups share no common FS_X / FS_Y values.")
        return

    pct = 100.0 * (Mb - Ma) / Ma   # B vs A

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

    fig.clear()

    if view == 'heatmap':
        ax = fig.add_subplot(111)
        if mode == 'diff':
            data = pct
            vmax = max(abs(np.nanmin(data)), abs(np.nanmax(data)))
            im = ax.imshow(data, origin='lower', aspect='auto',
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

    else:
        if view == 'diagonal':
            common = sorted(set(x_vals).intersection(y_vals))
            if not common:
                messagebox.showerror("No diagonal", "No FS values shared between X and Y.")
                return
            xs = np.array(common)
            ya = np.array([Ma[list(y_vals).index(v), list(x_vals).index(v)] for v in common])
            yb = np.array([Mb[list(y_vals).index(v), list(x_vals).index(v)] for v in common])
            xlabel = 'FS  [cm]  (square fields, X=Y)'
        elif view == 'row':
            try:
                y_pick = float(row_y_var.get())
            except ValueError:
                messagebox.showerror("Bad Y", "Row Y= must be numeric."); return
            iy = np.argmin(np.abs(y_vals - y_pick))
            if not np.isclose(y_vals[iy], y_pick, atol=0.01):
                messagebox.showerror("Y not found",
                                     f"Y={y_pick} not in common Y values: {list(y_vals)}"); return
            xs = x_vals
            ya = Ma[iy, :]; yb = Mb[iy, :]
            xlabel = f'FS_X [cm]  (Y = {y_vals[iy]:g})'
        else:  # col
            try:
                x_pick = float(col_x_var.get())
            except ValueError:
                messagebox.showerror("Bad X", "Col X= must be numeric."); return
            ix = np.argmin(np.abs(x_vals - x_pick))
            if not np.isclose(x_vals[ix], x_pick, atol=0.01):
                messagebox.showerror("X not found",
                                     f"X={x_pick} not in common X values: {list(x_vals)}"); return
            xs = y_vals
            ya = Ma[:, ix]; yb = Mb[:, ix]
            xlabel = f'FS_Y [cm]  (X = {x_vals[ix]:g})'

        if mode == 'overlay':
            ax = fig.add_subplot(111)
            ax.plot(xs, ya, 'o-', label=label_a)
            ax.plot(xs, yb, 's-', label=label_b)
            ax.set_xlabel(xlabel); ax.set_ylabel('Scp')
            ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        else:
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(xs, ya, 'o-', label=label_a)
            ax1.plot(xs, yb, 's-', label=label_b)
            ax1.set_ylabel('Scp')
            ax1.grid(True, alpha=0.3); ax1.legend(fontsize=8)
            ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
            pct_1d = 100.0 * (yb - ya) / ya
            ax2.plot(xs, pct_1d, 'o-', color='tab:red')
            ax2.axhline(0, color='k', lw=0.5)
            ax2.set_xlabel(xlabel); ax2.set_ylabel('% diff (B − A) / A')
            ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    canvas.draw()


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
    for i, c in enumerate(cols):
        ttk.Label(frame, text=f"{c}:").grid(row=0, column=2*i, sticky="e", padx=(4, 0))
        cb = ttk.Combobox(frame, textvariable=group_vars[c], width=12, state='normal')
        cb.grid(row=0, column=2*i+1, sticky="w", padx=(0, 4))
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

# Plot canvas
fig = plt.Figure(figsize=(9, 6), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=main)
canvas.get_tk_widget().grid(row=8, column=0, columnspan=4, sticky="nsew", pady=(8, 0))
toolbar = NavigationToolbar2Tk(canvas, main, pack_toolbar=False)
toolbar.grid(row=9, column=0, columnspan=4, sticky="we")

main.grid_rowconfigure(8, weight=1)
main.grid_columnconfigure(1, weight=1)


root.mainloop()
