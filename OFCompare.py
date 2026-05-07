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
DEFAULT_FILE  = (r"C:\Users\nknutson\OneDrive - Washington University in St. Louis"
                 r"\NGDS QA Consortium\Combined Consortium Data\Processed Combined Data"
                 r"\90 SSD\OF_Data\NGDSOutputFactorData.xlsx")


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


def _common_axes_maybe_interp(dfA, dfB):
    """Like _common_axes, but uses the global interpolate flag. When on,
    builds the union grid and interpolates each group's surface onto it
    (NaN outside each group's measured range). When off, uses intersection.
    """
    if not interp_var.get():
        return _common_axes(dfA, dfB)
    xA, yA, _ = _matrix(dfA)
    xB, yB, _ = _matrix(dfB)
    x_union = np.array(sorted(set(xA.tolist()) | set(xB.tolist())), dtype=float)
    y_union = np.array(sorted(set(yA.tolist()) | set(yB.tolist())), dtype=float)
    if len(x_union) == 0 or len(y_union) == 0:
        return None, None, None, None
    method = interp_method_var.get()
    Ma = _interp_to_grid(dfA, x_union, y_union, method=method)
    Mb = _interp_to_grid(dfB, x_union, y_union, method=method)
    return x_union, y_union, Ma, Mb


def _stack_traces(traces_xy, interpolate=False, method='pchip'):
    """Given a list of (xs, ys) arrays for individual traces, build a stacked
    matrix at the union of x values. When interpolate=False, missing cells
    stay NaN. When interpolate=True, each trace is interpolated onto the
    union grid (within its measured x range; NaN outside).
    Returns (xs_common, Y_stack) where Y_stack is shape (n_traces, n_x).
    """
    all_x = sorted(set().union(*[set(np.asarray(xs).tolist()) for xs, _ in traces_xy]))
    xs_common = np.array(all_x, dtype=float)
    Y = np.full((len(traces_xy), len(xs_common)), np.nan)
    if interpolate:
        from scipy.interpolate import PchipInterpolator
    for i, (xs, ys) in enumerate(traces_xy):
        if xs is None or len(xs) == 0:
            continue
        xs_arr = np.asarray(xs, dtype=float)
        ys_arr = np.asarray(ys, dtype=float)
        ok = ~np.isnan(ys_arr)
        if ok.sum() == 0:
            continue
        xs_arr, ys_arr = xs_arr[ok], ys_arr[ok]
        if interpolate and len(xs_arr) >= 2:
            order = np.argsort(xs_arr)
            xs_arr, ys_arr = xs_arr[order], ys_arr[order]
            # de-duplicate any repeated x's by averaging
            uniq_x, inv = np.unique(xs_arr, return_inverse=True)
            uniq_y = np.zeros_like(uniq_x, dtype=float)
            counts = np.zeros_like(uniq_x, dtype=float)
            np.add.at(uniq_y, inv, ys_arr)
            np.add.at(counts, inv, 1.0)
            uniq_y /= counts
            try:
                if method == 'pchip':
                    f = PchipInterpolator(uniq_x, uniq_y, extrapolate=False)
                else:
                    from scipy.interpolate import interp1d
                    f = interp1d(uniq_x, uniq_y, kind='linear',
                                 bounds_error=False, fill_value=np.nan)
                Y[i, :] = f(xs_common)
            except Exception:
                # Fall back to exact-match placement
                for x, y in zip(uniq_x, uniq_y):
                    j = int(np.argmin(np.abs(xs_common - float(x))))
                    if np.isclose(xs_common[j], float(x), atol=1e-6):
                        Y[i, j] = y
        else:
            for x, y in zip(xs_arr, ys_arr):
                j = int(np.argmin(np.abs(xs_common - float(x))))
                if np.isclose(xs_common[j], float(x), atol=1e-6):
                    Y[i, j] = y
    return xs_common, Y


def _interp_to_grid(df_trace, x_grid, y_grid, method='pchip'):
    """Interpolate one trace's Scp surface onto (x_grid × y_grid). Returns 2D
    array of shape (len(y_grid), len(x_grid)). NaN where (x, y) falls outside
    the trace's measured range — no extrapolation.
    """
    pivot = df_trace.pivot_table(index='FS_Y', columns='FS_X', values='Scp', aggfunc='mean')
    pivot = pivot.sort_index().sort_index(axis=1)
    xs = pivot.columns.values.astype(float)
    ys = pivot.index.values.astype(float)
    M  = pivot.values.astype(float)
    if M.size == 0 or len(xs) < 2 or len(ys) < 2:
        return np.full((len(y_grid), len(x_grid)), np.nan)

    from scipy.interpolate import RegularGridInterpolator
    # Drop rows / columns that are entirely NaN to keep the interpolator happy
    row_ok = ~np.all(np.isnan(M), axis=1)
    col_ok = ~np.all(np.isnan(M), axis=0)
    if not row_ok.any() or not col_ok.any():
        return np.full((len(y_grid), len(x_grid)), np.nan)
    ys2, xs2, M2 = ys[row_ok], xs[col_ok], M[np.ix_(row_ok, col_ok)]
    # Fill any remaining isolated NaN by row/column nearest fill (RGI doesn't accept NaN)
    if np.isnan(M2).any():
        df_M = pd.DataFrame(M2, index=ys2, columns=xs2)
        df_M = df_M.interpolate(axis=1, limit_direction='both').interpolate(axis=0, limit_direction='both')
        M2 = df_M.values

    try:
        rgi = RegularGridInterpolator((ys2, xs2), M2, method=method, bounds_error=False, fill_value=np.nan)
    except (ValueError, KeyError):
        rgi = RegularGridInterpolator((ys2, xs2), M2, method='linear', bounds_error=False, fill_value=np.nan)

    YY, XX = np.meshgrid(y_grid, x_grid, indexing='ij')
    pts = np.stack([YY.ravel(), XX.ravel()], axis=1)
    out = rgi(pts).reshape(len(y_grid), len(x_grid))
    # Mask points outside the original measured range
    out[(YY < ys2.min()) | (YY > ys2.max()) | (XX < xs2.min()) | (XX > xs2.max())] = np.nan
    return out


def _aggregate_traces_2d(traces, interpolate, method='pchip'):
    """Stack all expanded traces onto a common (FS_X × FS_Y) grid and return
    (x_grid, y_grid, stack) where stack has shape (n_traces, n_y, n_x). When
    `interpolate` is False, just place each trace's values on the union grid
    with NaN elsewhere.
    """
    all_x = sorted({float(v) for _, sub in traces for v in sub['FS_X'].unique()})
    all_y = sorted({float(v) for _, sub in traces for v in sub['FS_Y'].unique()})
    x_grid = np.array(all_x, dtype=float)
    y_grid = np.array(all_y, dtype=float)
    stack = np.full((len(traces), len(y_grid), len(x_grid)), np.nan)
    for i, (_, sub) in enumerate(traces):
        if interpolate:
            stack[i] = _interp_to_grid(sub, x_grid, y_grid, method=method)
        else:
            piv = sub.pivot_table(index='FS_Y', columns='FS_X', values='Scp', aggfunc='mean')
            for ry, yv in enumerate(y_grid):
                if yv not in piv.index:
                    continue
                for cx, xv in enumerate(x_grid):
                    if xv not in piv.columns:
                        continue
                    stack[i, ry, cx] = piv.loc[yv, xv]
    return x_grid, y_grid, stack


def _stats_per_cell(stack):
    """Compute per-cell stats across axis 0 of stack. Returns dict of 2D arrays."""
    Y = stack
    with np.errstate(all='ignore'):
        mean   = np.nanmean(Y, axis=0)
        median = np.nanmedian(Y, axis=0)
        std    = np.nanstd(Y, axis=0, ddof=1)
        cv_pct = 100.0 * std / np.where(np.abs(mean) > 0, mean, np.nan)
        q25    = np.nanpercentile(Y, 25, axis=0)
        q75    = np.nanpercentile(Y, 75, axis=0)
        iqr    = q75 - q25
        mad    = np.nanmedian(np.abs(Y - median[None, :, :]), axis=0)
        ymin   = np.nanmin(Y, axis=0)
        ymax   = np.nanmax(Y, axis=0)
        n      = np.sum(~np.isnan(Y), axis=0)
    return {'Mean': mean, 'Median': median, 'Std': std, 'CV_pct': cv_pct,
            'Q25': q25, 'Q75': q75, 'IQR': iqr, 'MAD': mad,
            'Min': ymin, 'Max': ymax, 'N': n}


def _print_diff_matrix(label_a, label_b, x_vals, y_vals, Ma, Mb):
    """Print A, B and % diff matrices to stdout (Y rows × X columns)."""
    pct = 100.0 * (Mb - Ma) / Ma
    print(f"\n{'='*78}")
    print(f"  Difference matrix    B ({label_b})  vs  A ({label_a})")
    print(f"{'='*78}")
    fmt_scp  = lambda v: f"{v:.4f}"
    fmt_diff = lambda v: f"{v:+.2f}"
    for sname, mat, fmt in (('A (Scp)', Ma, fmt_scp),
                            ('B (Scp)', Mb, fmt_scp),
                            ('% diff (B-A)/A', pct, fmt_diff)):
        df_out = pd.DataFrame(mat, index=y_vals, columns=x_vals)
        df_out.index.name = 'FS_Y \\ FS_X'
        print(f"\n-- {sname} --")
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None,
                               'display.width', 250):
            print(df_out.to_string(float_format=fmt))
    print('=' * 78)


def _print_stats_2d(traces, group_label=''):
    """Print per-stat 2D matrices (Y rows × X columns) to stdout for the given
    list of expanded traces. Honors the global Interpolate-to-common-grid flag.
    """
    if not traces or len(traces) < 2:
        return
    interpolate = interp_var.get()
    method      = interp_method_var.get()
    x_grid, y_grid, stack = _aggregate_traces_2d(traces, interpolate, method=method)
    stats = _stats_per_cell(stack)
    interp_tag = (f"PCHIP" if method == 'pchip' else "Linear") if interpolate else "no interp"
    print(f"\n{'='*78}")
    print(f"  Stats — 2D matrices   ({group_label}, {len(traces)} traces, {interp_tag})")
    print(f"{'='*78}")
    stat_order = ['Mean', 'Std', 'CV_pct', 'Median', 'IQR', 'MAD', 'Min', 'Max', 'N']
    for sname in stat_order:
        df_out = pd.DataFrame(stats[sname], index=y_grid, columns=x_grid)
        df_out.index.name = 'FS_Y \\ FS_X'
        print(f"\n-- {sname} --")
        fmt = (lambda v: f"{int(v):d}") if sname == 'N' else (lambda v: f"{v:.4f}")
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None,
                               'display.width', 250):
            print(df_out.to_string(float_format=fmt))
    print('=' * 78)


def _compute_stats(Y_stack):
    """Compute per-column stats on Y_stack (rows = traces, cols = x).
    Returns dict of arrays keyed by name. NaNs are ignored.
    """
    Y = np.asarray(Y_stack, dtype=float)
    with np.errstate(all='ignore'):
        mean   = np.nanmean(Y, axis=0)
        median = np.nanmedian(Y, axis=0)
        std    = np.nanstd(Y, axis=0, ddof=1)
        cv_pct = 100.0 * std / np.where(np.abs(mean) > 0, mean, np.nan)
        q25    = np.nanpercentile(Y, 25, axis=0)
        q75    = np.nanpercentile(Y, 75, axis=0)
        iqr    = q75 - q25
        # MAD around the median
        mad    = np.nanmedian(np.abs(Y - median[None, :]), axis=0)
        ymin   = np.nanmin(Y, axis=0)
        ymax   = np.nanmax(Y, axis=0)
        n      = np.sum(~np.isnan(Y), axis=0)
    return {'mean': mean, 'median': median, 'std': std, 'cv_pct': cv_pct,
            'q25': q25, 'q75': q75, 'iqr': iqr, 'mad': mad,
            'min': ymin, 'max': ymax, 'n': n}


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
file_a_var  = tk.StringVar(master=root, value=DEFAULT_FILE)
file_b_var  = tk.StringVar(master=root, value=DEFAULT_FILE)
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
diff_lim_var = tk.StringVar(master=root, value='6')   # ± y-axis / colorbar limit for % diff plots

# Stat overlay flags (apply when ≥2 traces exist)
stat_mean    = tk.BooleanVar(master=root, value=False)
stat_median  = tk.BooleanVar(master=root, value=False)
stat_std     = tk.BooleanVar(master=root, value=False)
stat_cv      = tk.BooleanVar(master=root, value=False)
stat_iqr     = tk.BooleanVar(master=root, value=False)
stat_minmax  = tk.BooleanVar(master=root, value=False)
stat_curves  = tk.BooleanVar(master=root, value=True)   # show individual traces too
std_k_var    = tk.StringVar(master=root, value='2')     # ± k × std band width
iqr_k_var    = tk.StringVar(master=root, value='1.5')   # k × IQR Tukey fence (Q25 - k*IQR, Q75 + k*IQR)

# Interpolation to common grid for stats aggregation
interp_var   = tk.BooleanVar(master=root, value=False)   # off → use raw cell values, NaN where missing
interp_method_var = tk.StringVar(master=root, value='pchip')   # 'pchip' or 'linear'

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


class MultiSelectCombo(ttk.Frame):
    """Read-only entry + dropdown button that opens a checkbox popup.
    Mimics the ttk.Combobox API for cb['values'] = list, cb.get(), cb.set().
    Stores selection as a comma-separated string in the bound StringVar:
      ''       → no constraint (treated like All in filters)
      'All'    → all current values selected
      'A'      → single value
      'A, B'   → multi-select
    """
    def __init__(self, parent, textvariable, width=12):
        super().__init__(parent)
        self._var = textvariable
        self._values = []
        self.entry = ttk.Entry(self, textvariable=self._var, width=max(4, width - 2), state='readonly')
        self.entry.pack(side='left', fill='x', expand=True)
        self.btn = ttk.Button(self, text='▼', width=2, command=self._open_popup)
        self.btn.pack(side='left')
        self.entry.bind('<Button-1>', lambda e: self._open_popup())

    def __setitem__(self, key, value):
        if key == 'values':
            self._values = list(value)
        else:
            super().__setitem__(key, value)

    def __getitem__(self, key):
        if key == 'values':
            return self._values
        return super().__getitem__(key)

    def get(self):
        return self._var.get()

    def set(self, text):
        self._var.set(text)

    def _parse_current(self):
        s = self._var.get().strip()
        opts = [v for v in self._values if v != ALL_TOKEN]
        if not s:
            return set()
        if s == ALL_TOKEN:
            return set(opts)
        return set(t.strip() for t in s.split(','))

    def _open_popup(self):
        opts = [v for v in self._values if v != ALL_TOKEN]
        if not opts:
            return
        dlg = tk.Toplevel(self.winfo_toplevel())
        dlg.title('Select')
        dlg.transient(self.winfo_toplevel())
        dlg.grab_set()
        current = self._parse_current()
        check_vars = {}
        container = ttk.Frame(dlg, padding=4)
        container.pack(fill='both', expand=True)
        for v in opts:
            cv = tk.BooleanVar(value=(v in current))
            check_vars[v] = cv
            ttk.Checkbutton(container, text=str(v), variable=cv).pack(anchor='w', padx=4, pady=1)
        btns = ttk.Frame(dlg, padding=4)
        btns.pack(fill='x')
        def _all():
            for cv in check_vars.values(): cv.set(True)
        def _none():
            for cv in check_vars.values(): cv.set(False)
        def _ok():
            sel = [v for v in opts if check_vars[v].get()]
            if not sel:
                self._var.set('')
            elif len(opts) > 1 and len(sel) == len(opts):
                self._var.set(ALL_TOKEN)
            elif len(sel) == 1:
                self._var.set(sel[0])
            else:
                self._var.set(', '.join(sel))
            dlg.destroy()
        ttk.Button(btns, text='All',  command=_all ).pack(side='left',  padx=2)
        ttk.Button(btns, text='None', command=_none).pack(side='left',  padx=2)
        ttk.Button(btns, text='OK',   command=_ok  ).pack(side='right', padx=2)
        self.update_idletasks()
        dlg.geometry(f'+{self.winfo_rootx()}+{self.winfo_rooty() + self.winfo_height()}')


def _split_filter_value(v):
    """Return list of values from a filter StringVar string (comma-separated)."""
    return [s.strip() for s in v.split(',') if s.strip()]


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


def _apply_filter(df, col, var_value):
    """Apply a single filter column to df given the StringVar string value.
    Empty or 'All' = no constraint. Comma-separated = multi-select.
    """
    v = (var_value or '').strip()
    if not v or v == ALL_TOKEN:
        return df
    values = _split_filter_value(v)
    if not values:
        return df
    if col in ('SSD', 'Depth'):
        try:
            targets = [float(s) for s in values]
        except ValueError:
            return df
        mask = pd.Series([False] * len(df), index=df.index)
        for t in targets:
            mask |= np.isclose(df[col].astype(float), t)
        return df[mask]
    return df[df[col].astype(str).isin(values)]


def _refilter_options(grp_letter, group_vars, df_source):
    """Cascading filter: each combobox's options reflect what's possible given the OTHER selections.
    Multi-select supported via comma-separated values; 'All' counts as no constraint.
    """
    if df_source is None:
        return
    for col in FILTER_COLS:
        sub = df_source
        for other in FILTER_COLS:
            if other == col:
                continue
            sub = _apply_filter(sub, other, group_vars[other].get())
        opts = sorted(sub[col].dropna().unique(), key=lambda v: (str(type(v)), v))
        opts_str = [_fmt_val(o) for o in opts]
        if col in _ALL_COLS and len(opts_str) > 1:
            opts_str = [ALL_TOKEN] + opts_str
        cb = _combos[grp_letter].get(col)
        if cb is None:
            continue
        cb['values'] = opts_str
        cur_vals = _split_filter_value(cb.get().strip())
        # Drop any selected values that are no longer in the option list
        valid = [v for v in cur_vals if v in opts_str]
        new_text = cb.get().strip()
        if new_text == ALL_TOKEN:
            pass  # keep
        elif valid != cur_vals:
            if not valid:
                cb.set('')
            elif len(valid) == 1:
                cb.set(valid[0])
            else:
                cb.set(', '.join(valid))


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
        df = _apply_filter(df, col, group_vars[col].get())
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

    # Surface views get their own square-ish figure; diff mode opens a 2nd window for the diff.
    if view == 'surface':
        figsize = (9, 7)
    else:
        figsize = (10, 6)
    fig = plt.figure(figsize=figsize)

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
        elif view == 'surface':
            if len(traces) > 1:
                messagebox.showinfo("Surface not multi-trace",
                                    "3D surface supports a single trace; pick specific values "
                                    "instead of 'All'.")
                return
            tdf = traces[0][1]
            x_vals, y_vals, M = _matrix(tdf)
            ax = fig.add_subplot(111, projection='3d')
            X, Y = np.meshgrid(x_vals, y_vals)
            ax.plot_surface(X, Y, M, cmap='viridis', alpha=0.85,
                            edgecolor='none', antialiased=True)
            # Scatter actual measured points
            ax.scatter(tdf['FS_X'], tdf['FS_Y'], tdf['Scp'],
                       c='black', s=12, depthshade=True)
            ax.set_xlabel('FS_X [cm]'); ax.set_ylabel('FS_Y [cm]'); ax.set_zlabel('Scp')
            ax.set_title(traces[0][0])
        else:
            ax = fig.add_subplot(111)
            xlabel_set = None
            traces_xy = []   # collect (xs, ys) for stat aggregation
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
                traces_xy.append((xs, ys, tlabel))
                if stat_curves.get():
                    ax.plot(xs, ys, 'o-', label=tlabel, alpha=0.45, linewidth=1)

            # Optional stat overlays when ≥2 traces
            if len(traces_xy) >= 2 and any(v.get() for v in
                    (stat_mean, stat_median, stat_std, stat_cv, stat_iqr, stat_minmax)):
                xs_common, Y = _stack_traces(
                    [(xs, ys) for xs, ys, _ in traces_xy],
                    interpolate=interp_var.get(),
                    method=interp_method_var.get())
                stats = _compute_stats(Y)
                # Print stats table to console
                _print_stats_2d(traces, base_label)
                try:
                    k_std = abs(float(std_k_var.get()))
                except ValueError:
                    k_std = 1.0
                try:
                    k_iqr = abs(float(iqr_k_var.get()))
                except ValueError:
                    k_iqr = 1.0
                if stat_std.get():
                    ax.fill_between(xs_common, stats['mean'] - k_std * stats['std'],
                                    stats['mean'] + k_std * stats['std'],
                                    color='gray', alpha=0.25, label=f'±{k_std:g} std')
                if stat_cv.get():
                    band = k_std * stats['mean'] * stats['cv_pct'] / 100.0
                    ax.fill_between(xs_common, stats['mean'] - band, stats['mean'] + band,
                                    color='tab:orange', alpha=0.18, label=f'±{k_std:g} CV%')
                if stat_iqr.get():
                    lo = stats['q25'] - k_iqr * stats['iqr'] if k_iqr != 1.0 else stats['q25']
                    hi = stats['q75'] + k_iqr * stats['iqr'] if k_iqr != 1.0 else stats['q75']
                    # When k_iqr == 1, just show the IQR (Q25..Q75); otherwise show Tukey fences.
                    label = 'IQR' if k_iqr == 1.0 else f'Q25-{k_iqr:g}·IQR .. Q75+{k_iqr:g}·IQR'
                    ax.fill_between(xs_common, lo, hi,
                                    color='tab:blue', alpha=0.18, label=label)
                if stat_minmax.get():
                    ax.plot(xs_common, stats['min'], ':', color='gray', linewidth=1, label='min/max')
                    ax.plot(xs_common, stats['max'], ':', color='gray', linewidth=1)
                if stat_mean.get():
                    ax.plot(xs_common, stats['mean'], '-', color='black', linewidth=2.2, label='mean')
                if stat_median.get():
                    ax.plot(xs_common, stats['median'], '--', color='tab:purple', linewidth=2, label='median')

            ax.set_xlabel(xlabel_set or 'FS [cm]')
            ax.set_ylabel('Scp')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8, loc='best')
        if view == 'surface':
            fig.subplots_adjust(left=0.02, right=0.96, top=0.95, bottom=0.04)
        else:
            fig.tight_layout()
        plt.show(block=False)
        plt.pause(0.001)
        return

    # ── both groups active: comparison view ──
    if view == 'heatmap':
        # Heatmap diff still needs overlapping axes.
        x_vals, y_vals, Ma, Mb = _common_axes_maybe_interp(dfA, dfB)
        if x_vals is None:
            messagebox.showerror("No overlap",
                                 "Heatmap diff needs overlapping FS_X / FS_Y values.")
            return
        ax = fig.add_subplot(111)
        if mode == 'diff':
            pct = 100.0 * (Mb - Ma) / Ma
            try:
                vmax = abs(float(diff_lim_var.get()))
            except ValueError:
                vmax = max(abs(np.nanmin(pct)), abs(np.nanmax(pct)))
            im = ax.imshow(pct, origin='lower', aspect='auto',
                           cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                           extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
            fig.colorbar(im, ax=ax, label='% diff (B − A) / A')
            ax.set_title(f'% diff   {label_b}  vs  {label_a}')
            _print_diff_matrix(label_a, label_b, x_vals, y_vals, Ma, Mb)
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

    if view == 'surface':
        x_vals, y_vals, Ma, Mb = _common_axes_maybe_interp(dfA, dfB)
        if x_vals is None:
            messagebox.showerror("No overlap",
                                 "Surface needs overlapping FS_X / FS_Y values.")
            return
        X, Y = np.meshgrid(x_vals, y_vals)
        if mode == 'diff':
            # Repurpose the original `fig` for the A&B overlay
            ax_top = fig.add_subplot(111, projection='3d')
            ax_top.plot_surface(X, Y, Ma, color='tab:blue',  alpha=0.55,
                                edgecolor='none', antialiased=True)
            ax_top.plot_surface(X, Y, Mb, color='tab:orange', alpha=0.55,
                                edgecolor='none', antialiased=True)
            ax_top.scatter(dfA['FS_X'], dfA['FS_Y'], dfA['Scp'], c='tab:blue',  s=10, depthshade=True)
            ax_top.scatter(dfB['FS_X'], dfB['FS_Y'], dfB['Scp'], c='tab:orange', s=10, depthshade=True)
            ax_top.set_xlabel('FS_X'); ax_top.set_ylabel('FS_Y'); ax_top.set_zlabel('Scp')
            ax_top.set_title(f'A (blue): {label_a}   |   B (orange): {label_b}')
            try:
                fig.canvas.manager.set_window_title('OF Surfaces — A & B')
            except Exception:
                pass

            # Separate figure for the % diff surface
            pct = 100.0 * (Mb - Ma) / Ma
            try:
                zlim = abs(float(diff_lim_var.get()))
            except ValueError:
                zlim = None
            fig2 = plt.figure(figsize=(9, 7))
            try:
                fig2.canvas.manager.set_window_title('OF % Difference')
            except Exception:
                pass
            ax_bot = fig2.add_subplot(111, projection='3d')
            surf = ax_bot.plot_surface(X, Y, pct, cmap='RdBu_r', alpha=1.0,
                                       edgecolor='none', antialiased=True,
                                       vmin=-zlim if zlim else None,
                                       vmax= zlim if zlim else None)
            ax_bot.plot_surface(X, Y, np.zeros_like(pct), color='black',
                                alpha=0.05, edgecolor='none')
            fig2.colorbar(surf, ax=ax_bot, shrink=0.7, pad=0.08, label='% diff')
            if zlim:
                ax_bot.set_zlim(-zlim, zlim)
            ax_bot.set_xlabel('FS_X'); ax_bot.set_ylabel('FS_Y'); ax_bot.set_zlabel('% diff')
            ax_bot.set_title(f'% diff   {label_b}  vs  {label_a}')
            fig2.subplots_adjust(left=0.04, right=0.92, top=0.94, bottom=0.06)
            _print_diff_matrix(label_a, label_b, x_vals, y_vals, Ma, Mb)
        else:  # overlay only — both surfaces on one plot
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, Ma, color='tab:blue',  alpha=0.5,
                            edgecolor='none', antialiased=True)
            ax.plot_surface(X, Y, Mb, color='tab:orange', alpha=0.5,
                            edgecolor='none', antialiased=True)
            ax.scatter(dfA['FS_X'], dfA['FS_Y'], dfA['Scp'], c='tab:blue',  s=10, depthshade=True)
            ax.scatter(dfB['FS_X'], dfB['FS_Y'], dfB['Scp'], c='tab:orange', s=10, depthshade=True)
            ax.set_xlabel('FS_X'); ax.set_ylabel('FS_Y'); ax.set_zlabel('Scp')
            ax.set_title(f'A (blue): {label_a}   |   B (orange): {label_b}')
        # Skip tight_layout for 3D — use subplots_adjust to push axes to fill the figure
        fig.subplots_adjust(left=0.02, right=0.96, top=0.96, bottom=0.04, hspace=0.05)
        plt.show(block=False)
        plt.pause(0.001)
        return

    # ── 1-D views (diagonal/row/col): build each group's data independently ──
    if interp_var.get():
        # Interpolate each group's surface onto the union grid before slicing
        xA_raw, yA_raw, _ = _matrix(dfA)
        xB_raw, yB_raw, _ = _matrix(dfB)
        x_union = np.array(sorted(set(xA_raw.tolist()) | set(xB_raw.tolist())), dtype=float)
        y_union = np.array(sorted(set(yA_raw.tolist()) | set(yB_raw.tolist())), dtype=float)
        method = interp_method_var.get()
        Ma = _interp_to_grid(dfA, x_union, y_union, method=method)
        Mb = _interp_to_grid(dfB, x_union, y_union, method=method)
        xA, yA = x_union, y_union
        xB, yB = x_union, y_union
    else:
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
        try:
            lim = abs(float(diff_lim_var.get()))
            ax2.set_ylim(-lim, lim)
        except ValueError:
            pass

    fig.tight_layout()
    plt.show(block=False)
    plt.pause(0.001)


def _expand_traces(df_source, group_vars):
    """Expand the group into one DataFrame per (Energy/SN/Detector) combo when
    those fields are 'All' or multi-selected (comma-separated). Single
    selections collapse to one trace. Returns a list of (label, sub_df).
    """
    df = _filter_group(df_source, group_vars)
    if df is None or df.empty:
        return []
    iter_cols = []
    for col in ('Energy', 'SN', 'Detector'):
        v = group_vars[col].get().strip()
        if v == ALL_TOKEN or ',' in v:
            iter_cols.append(col)
    if not iter_cols:
        return [('group', df)]
    out = []
    for combo in df[iter_cols].drop_duplicates().sort_values(iter_cols).values:
        sub = df.copy()
        for k, v in zip(iter_cols, combo):
            sub = sub[sub[k].astype(str) == str(v)]
        lbl = ' / '.join(str(v) for v in combo)
        out.append((lbl, sub))
    return out


def _export_stats():
    """Compute per-cell stats matrices across expanded traces for active groups.
    Writes one sheet per (group, stat) to xlsx, prints matrices to terminal,
    and shows them in a popup window.
    """
    if _df_a is None:
        messagebox.showinfo("Load first", "Load a file first.")
        return

    interpolate = interp_var.get()
    method      = interp_method_var.get()

    out_groups = []   # list of (letter, x_grid, y_grid, stats_dict)
    for letter, gv, df_src in (('A', group_a, _df_a), ('B', group_b, _df_b)):
        if not _group_active(gv) or df_src is None:
            continue
        traces = _expand_traces(df_src, gv)
        if not traces:
            continue
        x_grid, y_grid, stack = _aggregate_traces_2d(traces, interpolate, method=method)
        stats = _stats_per_cell(stack)
        out_groups.append((letter, x_grid, y_grid, stats, len(traces)))

    if not out_groups:
        messagebox.showinfo("Nothing to export", "No active groups with data.")
        return

    p = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                     filetypes=[("Excel", "*.xlsx")],
                                     initialfile="OF_stats.xlsx")
    if not p:
        return

    stat_order = ['Mean', 'Std', 'CV_pct', 'Median', 'IQR', 'MAD', 'Min', 'Max', 'N']

    def _matrix_df(x_grid, y_grid, mat):
        df = pd.DataFrame(mat, index=y_grid, columns=x_grid)
        df.index.name = 'FS_Y \\ FS_X'
        return df

    def _slice_df(x_grid, y_grid, stats):
        """Extract a 1D slice from each stat matrix based on the current view.
        Returns (slice_label, DataFrame) or (None, None) if view is heatmap.
        """
        view = view_var.get()
        if view == 'heatmap':
            return None, None
        if view == 'diagonal':
            common = [v for v in x_grid if v in y_grid]
            if not common:
                return None, None
            xs = np.array(common, dtype=float)
            data = {}
            for sname in stat_order:
                M = stats[sname]
                vals = np.array([M[list(y_grid).index(v), list(x_grid).index(v)] for v in common])
                data[sname] = vals
            df = pd.DataFrame(data, index=xs)
            df.index.name = 'FS  (X=Y)'
            return f"diagonal", df
        if view == 'row':
            try:
                y_pick = float(row_y_var.get())
            except ValueError:
                return None, None
            iy = int(np.argmin(np.abs(np.asarray(y_grid) - y_pick)))
            if not np.isclose(y_grid[iy], y_pick, atol=0.01):
                return None, None
            data = {sname: stats[sname][iy, :] for sname in stat_order}
            df = pd.DataFrame(data, index=np.asarray(x_grid, dtype=float))
            df.index.name = f'FS_X  (Y={y_grid[iy]:g})'
            return f"row_Y{y_grid[iy]:g}", df
        # col
        try:
            x_pick = float(col_x_var.get())
        except ValueError:
            return None, None
        ix = int(np.argmin(np.abs(np.asarray(x_grid) - x_pick)))
        if not np.isclose(x_grid[ix], x_pick, atol=0.01):
            return None, None
        data = {sname: stats[sname][:, ix] for sname in stat_order}
        df = pd.DataFrame(data, index=np.asarray(y_grid, dtype=float))
        df.index.name = f'FS_Y  (X={x_grid[ix]:g})'
        return f"col_X{x_grid[ix]:g}", df

    # Compute diff matrices when both groups are present (uses interpolation flag)
    diff_block = None
    if len(out_groups) == 2:
        x_d, y_d, Ma, Mb = _common_axes_maybe_interp(
            _filter_group(_df_a, group_a),
            _filter_group(_df_b, group_b))
        if x_d is not None:
            diff_pct = 100.0 * (Mb - Ma) / Ma
            diff_block = {'x': x_d, 'y': y_d, 'A': Ma, 'B': Mb, 'pct': diff_pct}

    try:
        with pd.ExcelWriter(p, engine='openpyxl') as w:
            for letter, x_grid, y_grid, stats, _ in out_groups:
                for sname in stat_order:
                    df_out = _matrix_df(x_grid, y_grid, stats[sname])
                    df_out.to_excel(w, sheet_name=f"{letter}_{sname}"[:31])
                slice_tag, slice_df = _slice_df(x_grid, y_grid, stats)
                if slice_df is not None:
                    slice_df.to_excel(w, sheet_name=f"{letter}_{slice_tag}"[:31])
            if diff_block is not None:
                _matrix_df(diff_block['x'], diff_block['y'], diff_block['A']  ).to_excel(w, sheet_name='Diff_A')
                _matrix_df(diff_block['x'], diff_block['y'], diff_block['B']  ).to_excel(w, sheet_name='Diff_B')
                _matrix_df(diff_block['x'], diff_block['y'], diff_block['pct']).to_excel(w, sheet_name='Diff_pct')
    except Exception as e:
        messagebox.showerror("Export failed", str(e))
        return

    # Build a single text blob for both terminal and popup
    lines = []
    interp_tag = (f"PCHIP" if method == 'pchip' else "Linear") if interpolate else "no interp"
    for letter, x_grid, y_grid, stats, n_traces in out_groups:
        lines.append("=" * 78)
        lines.append(f"  Group {letter}   ({n_traces} traces, {interp_tag})")
        lines.append("=" * 78)
        for sname in stat_order:
            df_out = _matrix_df(x_grid, y_grid, stats[sname])
            lines.append(f"\n-- {sname} --")
            fmt = (lambda v: f"{int(v):d}") if sname == 'N' else (lambda v: f"{v:.4f}")
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.width', 250):
                lines.append(df_out.to_string(float_format=fmt))
        # 1D slice for current view
        slice_tag, slice_df = _slice_df(x_grid, y_grid, stats)
        if slice_df is not None:
            lines.append(f"\n-- Slice: {slice_tag} --")
            def _fmt(v):
                # CV column gets a different format; N is integer; rest 4dp
                return f"{v:.4f}" if not (np.isnan(v)) else "nan"
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.width', 250):
                lines.append(slice_df.to_string(float_format=_fmt))
        lines.append("")
    if diff_block is not None:
        lines.append("=" * 78)
        lines.append("  Difference matrices  B vs A")
        lines.append("=" * 78)
        for sname, mat, fmt in (('A (Scp)', diff_block['A'], lambda v: f"{v:.4f}"),
                                ('B (Scp)', diff_block['B'], lambda v: f"{v:.4f}"),
                                ('% diff (B-A)/A', diff_block['pct'], lambda v: f"{v:+.2f}")):
            df_out = _matrix_df(diff_block['x'], diff_block['y'], mat)
            lines.append(f"\n-- {sname} --")
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None,
                                   'display.width', 250):
                lines.append(df_out.to_string(float_format=fmt))
        lines.append("")
    text = "\n".join(lines)
    print("\n" + text)

    # Popup window
    win = tk.Toplevel(root)
    win.title("OF stats")
    txt = tk.Text(win, wrap='none', font=('Courier', 9), width=110, height=35)
    txt.grid(row=0, column=0, sticky='nsew')
    vsb = ttk.Scrollbar(win, orient='vertical',   command=txt.yview)
    hsb = ttk.Scrollbar(win, orient='horizontal', command=txt.xview)
    txt.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
    vsb.grid(row=0, column=1, sticky='ns')
    hsb.grid(row=1, column=0, sticky='ew')
    win.grid_rowconfigure(0, weight=1); win.grid_columnconfigure(0, weight=1)
    txt.insert('1.0', text)
    txt.configure(state='disabled')

    messagebox.showinfo("Exported", f"Wrote stats matrices to {p}")


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
        cb = MultiSelectCombo(frame, textvariable=group_vars[c], width=14)
        cb.grid(row=0, column=2*i+1, sticky="w", padx=(0, 4))
        # Refilter cascade when the StringVar changes (popup OK closes, sets var)
        group_vars[c].trace_add('write', lambda *_a: _on_change())
        _combos[grp_letter][c] = cb


_make_filter_frame('A', group_a, 'A', 4)
_make_filter_frame('B', group_b, 'B', 5)

# View
view_frame = ttk.LabelFrame(main, text="View", padding=4)
view_frame.grid(row=6, column=0, columnspan=4, sticky="we", pady=(6, 0))
ttk.Radiobutton(view_frame, text="2D heatmap", variable=view_var, value='heatmap').pack(side='left', padx=4)
ttk.Radiobutton(view_frame, text="3D surface", variable=view_var, value='surface').pack(side='left', padx=4)
ttk.Radiobutton(view_frame, text="Diagonal (X=Y)", variable=view_var, value='diagonal').pack(side='left', padx=4)
ttk.Radiobutton(view_frame, text="Row (fixed Y)", variable=view_var, value='row').pack(side='left', padx=4)
ttk.Entry(view_frame, textvariable=row_y_var, width=6).pack(side='left')
ttk.Radiobutton(view_frame, text="Col (fixed X)", variable=view_var, value='col').pack(side='left', padx=(12, 4))
ttk.Entry(view_frame, textvariable=col_x_var, width=6).pack(side='left')

# Stats frame (overlay across traces when "All" expansion produces ≥2 subgroups)
stats_frame = ttk.LabelFrame(main, text="Stats overlay (multi-trace)", padding=4)
stats_frame.grid(row=7, column=0, columnspan=4, sticky="we", pady=(4, 0))
ttk.Checkbutton(stats_frame, text="Curves",   variable=stat_curves ).pack(side='left', padx=4)
ttk.Checkbutton(stats_frame, text="Mean",     variable=stat_mean   ).pack(side='left', padx=4)
ttk.Checkbutton(stats_frame, text="Median",   variable=stat_median ).pack(side='left', padx=4)
ttk.Checkbutton(stats_frame, text="±Std",     variable=stat_std    ).pack(side='left', padx=4)
ttk.Checkbutton(stats_frame, text="±CV%",     variable=stat_cv     ).pack(side='left', padx=4)
ttk.Label(stats_frame, text="k=").pack(side='left')
ttk.Entry(stats_frame, textvariable=std_k_var, width=4).pack(side='left', padx=(0, 4))
ttk.Checkbutton(stats_frame, text="IQR",      variable=stat_iqr    ).pack(side='left', padx=4)
ttk.Label(stats_frame, text="k=").pack(side='left')
ttk.Entry(stats_frame, textvariable=iqr_k_var, width=4).pack(side='left', padx=(0, 4))
ttk.Checkbutton(stats_frame, text="Min/Max",  variable=stat_minmax ).pack(side='left', padx=4)
ttk.Separator(stats_frame, orient='vertical').pack(side='left', fill='y', padx=8)
ttk.Checkbutton(stats_frame, text="Interp to common grid", variable=interp_var).pack(side='left', padx=4)
ttk.Combobox(stats_frame, textvariable=interp_method_var, values=['pchip', 'linear'],
             width=8, state='readonly').pack(side='left', padx=2)
ttk.Button(stats_frame, text="Export stats…", command=lambda: _export_stats()).pack(side='right', padx=4)

# Compare
cmp_frame = ttk.LabelFrame(main, text="Compare", padding=4)
cmp_frame.grid(row=8, column=0, columnspan=4, sticky="we", pady=(4, 0))
ttk.Radiobutton(cmp_frame, text="% difference (B vs A)", variable=compare_var, value='diff').pack(side='left', padx=4)
ttk.Radiobutton(cmp_frame, text="Overlay only",          variable=compare_var, value='overlay').pack(side='left', padx=4)
ttk.Label(cmp_frame, text="Diff range ±%:").pack(side='left', padx=(12, 2))
ttk.Entry(cmp_frame, textvariable=diff_lim_var, width=5).pack(side='left')
ttk.Button(cmp_frame, text="Plot", command=_plot).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear all", command=_clear_all_filters).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear B",   command=lambda: _clear_group('B')).pack(side='right', padx=4)
ttk.Button(cmp_frame, text="Clear A",   command=lambda: _clear_group('A')).pack(side='right', padx=4)

# Plots open in their own matplotlib window — see _plot() for figure creation.

main.grid_rowconfigure(8, weight=1)
main.grid_columnconfigure(1, weight=1)


root.mainloop()
