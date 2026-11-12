# Script to make AVG PDDs from multiple data sets and output the avg and plot that plus a confidence interval
import os
import webbrowser
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Patch as _MPatch
import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.stats import norm as _norm, t as _t

# Function to handle file selection and populate sheet names

# Function to retrieve clinical criteria from user input (with defaults)
def get_clinical_criteria():
    try:
        dose_tolerance = float(dose_entry.get())
    except ValueError:
        dose_tolerance = 1.0  # Default to 1% if input is invalid

    try:
        pos_tolerance = float(pos_entry.get())
    except ValueError:
        pos_tolerance = 1.0  # Default to 1 mm if input is invalid

    return dose_tolerance, pos_tolerance

def select_file():
    # Open file dialog to select an Excel file
    file_path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)

    if file_path:
        try:
            # Read the entire Excel file once, store in a global or local variable
            global excel_data
            excel_data = pd.read_excel(file_path, sheet_name=None)  # Read all sheets at once
            
            # Extract sheet names and populate the listbox
            sheets = sorted(excel_data.keys())
            sheet_listbox.delete(0, tk.END)
            for sheet in sheets:
                sheet_listbox.insert(tk.END, sheet)
            sheet_listbox.selection_set(0, tk.END)  # Select all by default
            
            # Populate field sizes, scan types, and depths
            populate_fsl()
            populate_scl()
            populate_dl()
        except Exception as e:
            print(f"Error reading sheets: {e}")



def populate_scl():
    try:
        selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
        if not selected_sheets:
            print("No sheets selected!")
            return

        common_scl = None
        for sheet in selected_sheets:
            df = excel_data[sheet]  # Use preloaded data from `select_file()`
            scl_current = set(df['Axis'].unique())

            if common_scl is None:
                common_scl = scl_current
            else:
                common_scl = common_scl.intersection(scl_current)

        # Update the listbox
        scl_listbox.delete(0, tk.END)
        for scl in sorted(common_scl):
            scl_listbox.insert(tk.END, scl)
        scl_listbox.selection_set(0, tk.END)  # Select all by default

    except Exception as e:
        print(f"Error loading scan types: {e}")




def populate_fsl():
    try:
        selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
        if not selected_sheets:
            print("No sheets selected!")
            return

        common_fsl = None
        for sheet in selected_sheets:
            df = excel_data[sheet]  # Use preloaded data from `select_file()`
            fsl_current = set(df['FS'].unique())

            if common_fsl is None:
                common_fsl = fsl_current
            else:
                common_fsl = common_fsl.intersection(fsl_current)

        # Update the listbox
        fsl_listbox.delete(0, tk.END)
        for fsl in sorted(common_fsl):
            fsl_listbox.insert(tk.END, fsl)
        fsl_listbox.selection_set(0, tk.END)

    except Exception as e:
        print(f"Error populating field sizes: {e}")

def populate_dl():
    try:
        selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
        if not selected_sheets:
            print("No sheets selected!")
            return

        common_depths = None
        for sheet in selected_sheets:
            df = excel_data[sheet]  # Use preloaded data from `select_file()`
            if 'Depth' in df.columns:
                depths = set(df['Depth'].unique())

                if common_depths is None:
                    common_depths = depths
                else:
                    common_depths = common_depths.intersection(depths)

        # If no common depths found (e.g. Z-only scans), clear the list silently
        if not common_depths:
            depth_listbox.delete(0, tk.END)
            return

        # Update the listbox
        depth_listbox.delete(0, tk.END)
        for depth in sorted(common_depths):
            depth_listbox.insert(tk.END, depth)
        depth_listbox.selection_set(0, tk.END)

    except Exception as e:
        print(f"Error loading depths: {e}")

def save_average():
    global avg_dose_global, fixed_pos_global, fsl, scl, dl

    if 'avg_dose_global' not in globals() or not avg_dose_global:
        print("Average data does not exist. Run make_avg first to calculate averages.")
        return

    selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]

    rows = []
    for combination_key, avg_dose in avg_dose_global.items():
        fs, axis, depth = combination_key
        fixed_pos = fixed_pos_global[combination_key]

        # Pull Energy and SSD from the first source sheet that has data for this combination
        energy, ssd = '', None
        for sheet in selected_sheets:
            df_src = excel_data.get(sheet)
            if df_src is None:
                continue
            sub = df_src[(df_src['FS'] == fs) & (df_src['Axis'] == axis)]
            if axis != 'Z' and depth is not None:
                sub = sub[sub['Depth'] == depth]
            if sub.empty:
                continue
            if 'Energy' in sub.columns:
                vals = sub['Energy'].dropna()
                energy = str(vals.iloc[0]) if len(vals) else ''
            if 'SSD' in sub.columns:
                vals = sub['SSD'].dropna()
                ssd = float(vals.mean()) if len(vals) else None
            break

        detector_label = 'AVG: ' + ', '.join(selected_sheets)
        row_depth = float(depth) if (depth is not None and axis != 'Z') else 0.0
        for pos_val, dose_val in zip(fixed_pos, avg_dose):
            rows.append({
                'Depth':    row_depth,
                'Pos':      round(float(pos_val),  4),
                'Dose':     round(float(dose_val), 4),
                'FS':       fs,
                'Axis':     axis,
                'Energy':   energy,
                'Detector': detector_label,
                'SSD':      ssd,
            })

    if not rows:
        print("No average data to save.")
        return

    COLS = ['Depth', 'Pos', 'Dose', 'FS', 'Axis', 'Energy', 'Detector', 'SSD']
    df_all      = pd.DataFrame(rows, columns=COLS)
    df_depth    = df_all[df_all['Axis'] == 'Z'].copy()
    df_profiles = df_all[df_all['Axis'] != 'Z'].copy()

    metadata = {
        'Date Created': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'Sheets Used':  ', '.join(selected_sheets),
        'Min Datasets': min_datasets_var.get(),
        'Field Sizes':  ', '.join(map(str, fsl)) if 'fsl' in globals() else '',
        'Scan Types':   ', '.join(scl)           if 'scl' in globals() else '',
        'Depths':       ', '.join(map(str, dl))  if ('dl' in globals() and dl) else 'All Depths',
    }
    df_metadata = pd.DataFrame(list(metadata.items()), columns=['Description', 'Value'])

    output_file = filedialog.asksaveasfilename(
        defaultextension=".xlsx",
        filetypes=[("Excel Files", "*.xlsx")],
        title="Save Average Data As")
    if not output_file:
        return

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_depth.to_excel(   writer, sheet_name='Depth Scans', index=False)
        df_profiles.to_excel(writer, sheet_name='Profiles',    index=False)
        df_metadata.to_excel(writer, sheet_name='Metadata',    index=False)

    print(f"Average data saved to {output_file}")

_manual_shifts      = {}   # {sheet_name: x_shift_cm}
avg_dose_global     = {}   # {combination: avg_array}
fixed_pos_global    = {}   # {combination: position_array}
std_dose_global     = {}   # {combination: std_array}
_raw_stacks_global  = {}   # {combination: {sheet_name: resampled_array}}


def _compute_dta(sheet_dose, avg_dose, fixed_pos, pos_tol_cm):
    """DTA (cm) at each position: smallest shift of avg that matches sheet_dose.
    Search is capped at ±(pos_tol_cm + 0.1 cm), matching the comp.py convention."""
    valid = ~np.isnan(avg_dose) & ~np.isnan(sheet_dose)
    if not valid.any():
        return np.full(len(fixed_pos), np.nan)
    avg_int = PchipInterpolator(fixed_pos[valid], avg_dose[valid], extrapolate=False)
    search   = np.arange(-(pos_tol_cm + 0.1), pos_tol_cm + 0.11, 0.01)
    avg_grid     = avg_int(fixed_pos[:, None] + search[None, :])   # (N, S)
    diffs        = np.abs(sheet_dose[:, None] - avg_grid)          # (N, S)
    diffs_filled = np.where(np.isnan(diffs), np.inf, diffs)        # replace NaN so argmin works
    best_ix      = np.argmin(diffs_filled, axis=1)
    dta_out      = np.abs(search[best_ix])
    # positions where the entire search row was NaN → cap at search boundary
    dta_out[np.all(np.isnan(diffs), axis=1)] = pos_tol_cm + 0.1
    return dta_out   # cm


def _compute_envelope(avg_dose, fixed_pos, dose_tolerance, pos_tolerance, axis):
    """Return (upper_envelope, lower_envelope) for the tolerance band."""
    tol_cm = pos_tolerance / 10.0
    if axis in ('X', 'Y', 'XY', 'YX'):
        sh_l = np.interp(fixed_pos, fixed_pos - tol_cm, avg_dose, left=np.nan, right=np.nan)
        sh_r = np.interp(fixed_pos, fixed_pos + tol_cm, avg_dose, left=np.nan, right=np.nan)
        upper = np.fmax(avg_dose + dose_tolerance, np.fmax(sh_l, sh_r))
        lower = np.fmin(avg_dose - dose_tolerance, np.fmin(sh_l, sh_r))
        upper = np.nan_to_num(upper, nan=avg_dose[-1])
        lower = np.nan_to_num(lower, nan=avg_dose[0])
    else:
        up_shift = np.interp(fixed_pos + tol_cm, fixed_pos, avg_dose, left=np.nan, right=np.nan)
        lo_shift = np.interp(fixed_pos - tol_cm, fixed_pos, avg_dose, left=np.nan, right=np.nan)
        upper = np.fmax(avg_dose + dose_tolerance, np.nan_to_num(up_shift, nan=avg_dose[-1]))
        lower = np.fmin(avg_dose - dose_tolerance, np.nan_to_num(lo_shift, nan=avg_dose[0]))
    return upper, lower


def _preprocess(pos, dose, axis, sheet):
    """Apply pre-processing steps to a single curve before plotting/averaging."""
    pos  = pos.copy().astype(float)
    dose = dose.copy().astype(float)
    is_profile = axis in ('X', 'Y', 'XY', 'YX')

    # Manual shift applies to all curve types
    shift = _manual_shifts.get(sheet, 0.0)
    if shift:
        pos = pos + shift

    if is_profile:
        # Auto-center: shift so FWHM midpoint = 0
        if auto_center_var.get():
            try:
                fine = np.linspace(float(pos.min()), float(pos.max()), 10000)
                fine_dose = PchipInterpolator(pos, dose)(fine)
                half_max = float(fine_dose.max()) / 2.0
                above = np.where(fine_dose >= half_max)[0]
                if len(above) >= 2:
                    pos = pos - (fine[above[0]] + fine[above[-1]]) / 2.0
            except Exception:
                pass

        # Y normalization: CAX (x=0) or max
        if cax_norm_var.get():
            try:
                d_cax = float(PchipInterpolator(pos, dose)(0.0))
                dose = dose / d_cax * 100.0 if d_cax > 0 else dose / dose.max() * 100.0
            except Exception:
                dose = dose / dose.max() * 100.0
        else:
            dose = dose / dose.max() * 100.0

    else:  # PDD (Z axis)
        mode = pdd_norm_var.get()
        if mode == "depth":
            try:
                d_ref = float(pdd_norm_depth_entry.get())
                d_at_ref = float(PchipInterpolator(pos, dose)(d_ref))
                dose = dose / d_at_ref * 100.0 if d_at_ref > 0 else dose / dose.max() * 100.0
            except Exception:
                dose = dose / dose.max() * 100.0
        elif mode == "TPS d10":
            try:
                d_at_10 = float(PchipInterpolator(pos, dose)(10.0))
                tps_val = float(tps_d10_entry.get())
                dose = dose / d_at_10 * tps_val if d_at_10 > 0 else dose / dose.max() * 100.0
            except Exception:
                dose = dose / dose.max() * 100.0
        else:
            dose = dose / dose.max() * 100.0

    return pos, dose


def _open_manual_shifts():
    selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
    if not selected_sheets:
        from tkinter import messagebox
        messagebox.showwarning("No sheets selected", "Select sheets in the list before opening manual shifts.")
        return

    dlg = tk.Toplevel(root)
    dlg.title("Manual Profile X Shifts")
    dlg.resizable(False, False)
    dlg.grab_set()

    ttk.Label(dlg, text="Sheet",         width=32).grid(row=0, column=0, padx=12, pady=4, sticky='w')
    ttk.Label(dlg, text="Shift [cm]",   width=12).grid(row=0, column=1, padx=12, pady=4, sticky='w')

    # Store Entry widgets directly (more reliable than StringVar.get() in Toplevel)
    entry_widgets = {}
    for i, sheet in enumerate(selected_sheets):
        ttk.Label(dlg, text=sheet, width=32).grid(row=i + 1, column=0, padx=12, pady=3, sticky='w')
        cur = _manual_shifts.get(sheet, 0.0)
        ent = ttk.Entry(dlg, width=10)
        if cur != 0.0:
            ent.insert(0, str(cur))
        ent.grid(row=i + 1, column=1, padx=12, pady=3)
        entry_widgets[sheet] = ent

    def _save(*_):
        for sheet, ent in entry_widgets.items():
            txt = ent.get().strip()
            try:
                _manual_shifts[sheet] = float(txt) if txt else 0.0
            except ValueError:
                _manual_shifts[sheet] = 0.0
        print("Manual shifts saved:", {k: v for k, v in _manual_shifts.items() if v != 0.0})
        dlg.destroy()

    def _clear():
        for ent in entry_widgets.values():
            ent.delete(0, tk.END)

    # Bind Enter key and window-close to save
    dlg.bind("<Return>", _save)
    dlg.protocol("WM_DELETE_WINDOW", _save)

    n = len(selected_sheets)
    btn_frame = ttk.Frame(dlg)
    btn_frame.grid(row=n + 1, column=0, columnspan=2, pady=8)
    ttk.Button(btn_frame, text="OK",       command=_save        ).pack(side='left', padx=6)
    ttk.Button(btn_frame, text="Clear All", command=_clear      ).pack(side='left', padx=6)
    ttk.Button(btn_frame, text="Cancel",   command=dlg.destroy  ).pack(side='left', padx=6)


def plot_data():
    global dl, fsl, scl, avg_dose_global, fixed_pos_global

    # Get selected field sizes, scan types, and depths
    selected_fsl_indices = fsl_listbox.curselection()
    selected_scl_indices = scl_listbox.curselection()
    selected_depth_indices = depth_listbox.curselection()

    fsl = [float(fsl_listbox.get(i)) for i in selected_fsl_indices]
    scl = [scl_listbox.get(i) for i in selected_scl_indices]
    dl = [float(depth_listbox.get(i)) for i in selected_depth_indices] if selected_depth_indices else []

    try:
        _ms = float(marker_size_entry.get())
    except ValueError:
        _ms = 3.0

    try:
        _fs = float(font_size_entry.get())
    except ValueError:
        _fs = 16.0
    plt.rcParams.update({'font.size': _fs})

    # Define marker styles to iterate over
    markers = ['o', 's', 'D', '+', 'x', '^', '*', 'p', 'h', 'H', '|', '_']
    
    # Get selected sheets from the listbox
    selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
    
    if not selected_sheets:
        print("No sheets selected!")
        return
    
    # Set up plotting for profiles (X, Y, XY, YX)
    fig_profiles, ax_profiles = None, None
    if any(axis in ['X', 'Y', 'XY', 'YX'] for axis in scl):
        fig_profiles, ax_profiles = plt.subplots(figsize=(10, 5), constrained_layout=True)
        ax_profiles.set_ylabel('Dose [%]')
        ax_profiles.set_xlabel('Off Axis Position [cm]')
        max_field_size = max(fsl) if fsl else 0
        try:
            _pxmin = float(profile_xmin_entry.get())
        except ValueError:
            _pxmin = -max_field_size / 2 - 2
            profile_xmin_entry.delete(0, tk.END)
            profile_xmin_entry.insert(0, f"{_pxmin:.2f}")
        try:
            _pxmax = float(profile_xmax_entry.get())
        except ValueError:
            _pxmax = max_field_size / 2 + 2
            profile_xmax_entry.delete(0, tk.END)
            profile_xmax_entry.insert(0, f"{_pxmax:.2f}")
        ax_profiles.set_xlim([_pxmin, _pxmax])
        ax_profiles.set_ylim([0, 110])

    # Set up plotting for Z profiles
    fig_z, ax_z = None, None
    if 'Z' in scl:
        fig_z = plt.figure(figsize=(10, 5), constrained_layout=True)
        ax_z = fig_z.add_subplot(111)
        ax_z.set_ylabel('Percentage Depth Dose [%]')
        ax_z.set_xlabel('Position [cm]')
        try:
            _zxmin = float(z_xmin_entry.get())
        except ValueError:
            _zxmin = 0
            z_xmin_entry.delete(0, tk.END)
            z_xmin_entry.insert(0, "0")
        try:
            _zxmax = float(z_xmax_entry.get())
        except ValueError:
            _zxmax = 30
            z_xmax_entry.delete(0, tk.END)
            z_xmax_entry.insert(0, "30")
        ax_z.set_xlim([_zxmin, _zxmax])
        ax_z.set_ylim([0, 110])

    # Define a color map for sheets
    color_map = plt.colormaps["tab20"].resampled(len(selected_sheets))

    highlight_sheet = highlight_var.get() if highlight_var.get() != "None" else None

    # Plot actual data from sheets
    for sheet_idx, sheet in enumerate(selected_sheets):
        df = excel_data[sheet].sort_values(by=['FS', 'Axis', 'Pos'])
        color  = color_map(sheet_idx)
        marker = markers[sheet_idx % len(markers)]
        is_highlight = (sheet == highlight_sheet)
        ls    = 'None'
        color = 'red'    if is_highlight else color_map(sheet_idx)
        alpha = 1.0      if is_highlight else 0.5
        ms    = _ms * 2  if is_highlight else _ms
        zord  = 5        if is_highlight else 2
        sheet_labeled = False
        for fs in fsl:
            for axis in scl:
                if axis in ['X', 'Y', 'XY', 'YX']:
                    for depth in dl:
                        if not df[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth)].empty:
                            pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Pos'].values
                            dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Dose'].values
                            pos, dose = _preprocess(pos, dose, axis, sheet)
                            if not sheet_labeled:
                                ax_profiles.plot(pos, dose, linestyle=ls, marker=marker, color=color, ms=ms, alpha=alpha, zorder=zord, label=sheet)
                                sheet_labeled = True
                            else:
                                ax_profiles.plot(pos, dose, linestyle=ls, marker=marker, color=color, ms=ms, alpha=alpha, zorder=zord, label=None)

                elif axis == 'Z':
                    if not df[(df['FS'] == fs) & (df['Axis'] == axis)].empty:
                        pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Pos'].values
                        dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Dose'].values
                        pos, dose = _preprocess(pos, dose, axis, sheet)
                        if not sheet_labeled:
                            ax_z.plot(pos, dose, linestyle=ls, marker=marker, color=color, ms=ms, alpha=alpha, zorder=zord, label=sheet)
                            sheet_labeled = True
                        else:
                            ax_z.plot(pos, dose, linestyle=ls, marker=marker, color=color, ms=ms, alpha=alpha, zorder=zord, label=None)

    # Now plot the averages if they exist
    try:
        if avg_dose_global:  # Only check if average exists
            for fs in fsl:
                for axis in scl:
                    if axis in ['X', 'Y', 'XY', 'YX']:
                        for depth in dl:
                            combination_key = (fs, axis, depth)
                            if combination_key in avg_dose_global:
                                avg_dose = avg_dose_global[combination_key]
                                fixed_pos = fixed_pos_global[combination_key]
                                # Only label the first average plot once
                                label = "Average" if fs == fsl[0] and axis == scl[0] and depth == dl[0] else None
                                ax_profiles.plot(fixed_pos, avg_dose, '-k', label=label, linewidth=2)

                                dose_tolerance, pos_tolerance = get_clinical_criteria()
                                upper_envelope, lower_envelope = _compute_envelope(
                                    avg_dose, fixed_pos, dose_tolerance, pos_tolerance, axis)
                                if show_band_var.get():
                                    ax_profiles.fill_between(fixed_pos, lower_envelope, upper_envelope, facecolor='0.85', edgecolor='none', alpha=1.0)
                                    ax_profiles.plot(fixed_pos, upper_envelope, '--k', linewidth=1)
                                    ax_profiles.plot(fixed_pos, lower_envelope, '--k', linewidth=1)
                                if show_sd_var.get() and combination_key in std_dose_global:
                                    try:
                                        _k = float(sd_mult_entry.get())
                                    except ValueError:
                                        _k = 2.326
                                    _std = std_dose_global[combination_key]
                                    ax_profiles.fill_between(fixed_pos, avg_dose - _k * _std, avg_dose + _k * _std, facecolor='0.75', edgecolor='none', alpha=1.0)
                                    ax_profiles.plot(fixed_pos, avg_dose + _k * _std, '--k', linewidth=1)
                                    ax_profiles.plot(fixed_pos, avg_dose - _k * _std, '--k', linewidth=1)
                                if show_ci_mean_var.get() and combination_key in std_dose_global:
                                    _n_ci = len(_raw_stacks_global.get(combination_key, {}))
                                    if _n_ci >= 2:
                                        try:
                                            _ci_pct_m = float(ci_entry.get())
                                        except ValueError:
                                            _ci_pct_m = 98.0
                                        _t_val = float(_t.ppf((1 + _ci_pct_m / 100) / 2, df=_n_ci - 1))
                                        _sem = std_dose_global[combination_key] / np.sqrt(_n_ci)
                                        ax_profiles.fill_between(fixed_pos, avg_dose - _t_val * _sem, avg_dose + _t_val * _sem, facecolor='0.60', edgecolor='none', alpha=1.0)
                                        ax_profiles.plot(fixed_pos, avg_dose + _t_val * _sem, '--k', linewidth=1)
                                        ax_profiles.plot(fixed_pos, avg_dose - _t_val * _sem, '--k', linewidth=1)

                    elif axis == 'Z':
                        combination_key = (fs, axis, None)
                        if combination_key in avg_dose_global:
                            avg_dose = avg_dose_global[combination_key]
                            fixed_pos = fixed_pos_global[combination_key]
                            # Only label the first Z average plot once
                            label = "Average" if fs == fsl[0] and axis == 'Z' else None
                            ax_z.plot(fixed_pos, avg_dose, '-k', label=label, linewidth=2)

                            dose_tolerance, pos_tolerance = get_clinical_criteria()
                            upper_envelope, lower_envelope = _compute_envelope(
                                avg_dose, fixed_pos, dose_tolerance, pos_tolerance, axis)
                            if show_band_var.get():
                                ax_z.fill_between(fixed_pos, lower_envelope, upper_envelope, facecolor='0.85', edgecolor='none', alpha=1.0)
                                ax_z.plot(fixed_pos, upper_envelope, '--k', linewidth=1)
                                ax_z.plot(fixed_pos, lower_envelope, '--k', linewidth=1)
                            if show_sd_var.get() and combination_key in std_dose_global:
                                try:
                                    _k = float(sd_mult_entry.get())
                                except ValueError:
                                    _k = 2.326
                                _std = std_dose_global[combination_key]
                                ax_z.fill_between(fixed_pos, avg_dose - _k * _std, avg_dose + _k * _std, facecolor='0.75', edgecolor='none', alpha=1.0)
                                ax_z.plot(fixed_pos, avg_dose + _k * _std, '--k', linewidth=1)
                                ax_z.plot(fixed_pos, avg_dose - _k * _std, '--k', linewidth=1)
                            if show_ci_mean_var.get() and combination_key in std_dose_global:
                                _n_ci = len(_raw_stacks_global.get(combination_key, {}))
                                if _n_ci >= 2:
                                    try:
                                        _ci_pct_m = float(ci_entry.get())
                                    except ValueError:
                                        _ci_pct_m = 98.0
                                    _t_val = float(_t.ppf((1 + _ci_pct_m / 100) / 2, df=_n_ci - 1))
                                    _sem = std_dose_global[combination_key] / np.sqrt(_n_ci)
                                    ax_z.fill_between(fixed_pos, avg_dose - _t_val * _sem, avg_dose + _t_val * _sem, facecolor='0.60', edgecolor='none', alpha=1.0)
                                    ax_z.plot(fixed_pos, avg_dose + _t_val * _sem, '--k', linewidth=1)
                                    ax_z.plot(fixed_pos, avg_dose - _t_val * _sem, '--k', linewidth=1)
    except NameError:
        print("Average data does not exist. Run make_avg first to calculate averages.")

    # --- Zoom insets ---
    def _populate_inset(ins, axis_filter):
        """Re-plot all data onto an inset axis (no labels/legend)."""
        for _si, _sh in enumerate(selected_sheets):
            _df  = excel_data[_sh].sort_values(by=['FS', 'Axis', 'Pos'])
            _col = 'red' if _sh == highlight_sheet else color_map(_si)
            _msi = _ms * 2 if _sh == highlight_sheet else _ms
            _alp = 1.0 if _sh == highlight_sheet else 0.5
            _mkr = markers[_si % len(markers)]
            for _fs in fsl:
                for _ax in [a for a in scl if axis_filter(a)]:
                    _dl_loop = dl if _ax in ('X', 'Y', 'XY', 'YX') else [None]
                    for _dep in _dl_loop:
                        if _dep is not None:
                            _sub = _df[(_df['FS'] == _fs) & (_df['Axis'] == _ax) & (_df['Depth'] == _dep)]
                        else:
                            _sub = _df[(_df['FS'] == _fs) & (_df['Axis'] == _ax)]
                        if _sub.empty:
                            continue
                        _p, _d = _preprocess(_sub['Pos'].values, _sub['Dose'].values, _ax, _sh)
                        ins.plot(_p, _d, linestyle='None', marker=_mkr, color=_col, ms=_msi, alpha=_alp)
        if avg_dose_global:
            _dt, _pt = get_clinical_criteria()
            try: _ki = float(sd_mult_entry.get())
            except ValueError: _ki = 2.326
            for _fs in fsl:
                for _ax in [a for a in scl if axis_filter(a)]:
                    _dep = None if _ax == 'Z' else None
                    _dl_loop = dl if _ax in ('X', 'Y', 'XY', 'YX') else [None]
                    for _dep in _dl_loop:
                        _key = (_fs, _ax, _dep)
                        if _key not in avg_dose_global:
                            continue
                        _avg, _fp = avg_dose_global[_key], fixed_pos_global[_key]
                        ins.plot(_fp, _avg, '-k', linewidth=2)
                        if show_band_var.get():
                            _ue, _le = _compute_envelope(_avg, _fp, _dt, _pt, _ax)
                            ins.fill_between(_fp, _le, _ue, facecolor='0.85', edgecolor='none')
                            ins.plot(_fp, _ue, '--k', linewidth=1)
                            ins.plot(_fp, _le, '--k', linewidth=1)
                        if show_sd_var.get() and _key in std_dose_global:
                            _st = std_dose_global[_key]
                            ins.fill_between(_fp, _avg - _ki*_st, _avg + _ki*_st, facecolor='0.75', edgecolor='none')
                            ins.plot(_fp, _avg + _ki*_st, '--k', linewidth=1)
                            ins.plot(_fp, _avg - _ki*_st, '--k', linewidth=1)
                        if show_ci_mean_var.get() and _key in std_dose_global:
                            _n_ci = len(_raw_stacks_global.get(_key, {}))
                            if _n_ci >= 2:
                                try:
                                    _ci_pct_m = float(ci_entry.get())
                                except ValueError:
                                    _ci_pct_m = 98.0
                                _t_val = float(_t.ppf((1 + _ci_pct_m / 100) / 2, df=_n_ci - 1))
                                _sem = std_dose_global[_key] / np.sqrt(_n_ci)
                                ins.fill_between(_fp, _avg - _t_val*_sem, _avg + _t_val*_sem, facecolor='0.60', edgecolor='none')
                                ins.plot(_fp, _avg + _t_val*_sem, '--k', linewidth=1)
                                ins.plot(_fp, _avg - _t_val*_sem, '--k', linewidth=1)

    try:
        _pix1, _pix2 = float(profile_inset_xmin.get()), float(profile_inset_xmax.get())
    except ValueError:
        _pix1 = _pix2 = None

    try:
        _zix1, _zix2 = float(pdd_inset_xmin.get()), float(pdd_inset_xmax.get())
    except ValueError:
        _zix1 = _zix2 = None

    if fig_profiles and show_profile_inset_var.get() and _pix1 is not None and _pix1 < _pix2:
        _p_bounds = list(_INSET_POSITIONS.get(profile_inset_pos_var.get(), [0.53, 0.37, 0.45, 0.58]))
        try: _p_bounds[0] = float(profile_inset_x0.get())
        except ValueError: pass
        try: _p_bounds[1] = float(profile_inset_y0.get())
        except ValueError: pass
        try: _p_bounds[2] = float(profile_inset_w.get())
        except ValueError: pass
        try: _p_bounds[3] = float(profile_inset_h.get())
        except ValueError: pass
        ins_p = ax_profiles.inset_axes(_p_bounds)
        _populate_inset(ins_p, lambda a: a in ('X', 'Y', 'XY', 'YX'))
        ins_p.set_xlim([_pix1, _pix2])
        try:
            ins_p.set_ylim([float(profile_inset_ymin.get()), float(profile_inset_ymax.get())])
        except ValueError:
            pass
        ins_p.tick_params(labelsize=max(8, _fs - 4))
        _, _p_connectors = ax_profiles.indicate_inset_zoom(ins_p, edgecolor='black')
        for _c in _p_connectors:
            _c.set_visible(True)

    if fig_z and show_pdd_inset_var.get() and _zix1 is not None and _zix1 < _zix2:
        _z_bounds = list(_INSET_POSITIONS.get(pdd_inset_pos_var.get(), [0.53, 0.37, 0.45, 0.58]))
        try: _z_bounds[0] = float(pdd_inset_x0.get())
        except ValueError: pass
        try: _z_bounds[1] = float(pdd_inset_y0.get())
        except ValueError: pass
        try: _z_bounds[2] = float(pdd_inset_w.get())
        except ValueError: pass
        try: _z_bounds[3] = float(pdd_inset_h.get())
        except ValueError: pass
        ins_z = ax_z.inset_axes(_z_bounds)
        _populate_inset(ins_z, lambda a: a == 'Z')
        ins_z.set_xlim([_zix1, _zix2])
        try:
            ins_z.set_ylim([float(pdd_inset_ymin.get()), float(pdd_inset_ymax.get())])
        except ValueError:
            pass
        ins_z.tick_params(labelsize=max(8, _fs - 4))
        _, _z_connectors = ax_z.indicate_inset_zoom(ins_z, edgecolor='black')
        for _c in _z_connectors:
            _c.set_visible(True)

    try:
        _pix1b, _pix2b = float(profile_inset2_xmin.get()), float(profile_inset2_xmax.get())
    except ValueError:
        _pix1b = _pix2b = None

    try:
        _zix1b, _zix2b = float(pdd_inset2_xmin.get()), float(pdd_inset2_xmax.get())
    except ValueError:
        _zix1b = _zix2b = None

    if fig_profiles and show_profile_inset2_var.get() and _pix1b is not None and _pix1b < _pix2b:
        _p2_bounds = list(_INSET_POSITIONS.get(profile_inset2_pos_var.get(), [0.02, 0.37, 0.45, 0.58]))
        try: _p2_bounds[0] = float(profile_inset2_x0.get())
        except ValueError: pass
        try: _p2_bounds[1] = float(profile_inset2_y0.get())
        except ValueError: pass
        try: _p2_bounds[2] = float(profile_inset2_w.get())
        except ValueError: pass
        try: _p2_bounds[3] = float(profile_inset2_h.get())
        except ValueError: pass
        ins_p2 = ax_profiles.inset_axes(_p2_bounds)
        _populate_inset(ins_p2, lambda a: a in ('X', 'Y', 'XY', 'YX'))
        ins_p2.set_xlim([_pix1b, _pix2b])
        try:
            ins_p2.set_ylim([float(profile_inset2_ymin.get()), float(profile_inset2_ymax.get())])
        except ValueError:
            pass
        ins_p2.tick_params(labelsize=max(8, _fs - 4))
        _, _p2_connectors = ax_profiles.indicate_inset_zoom(ins_p2, edgecolor='black')
        for _c in _p2_connectors:
            _c.set_visible(True)

    if fig_z and show_pdd_inset2_var.get() and _zix1b is not None and _zix1b < _zix2b:
        _z2_bounds = list(_INSET_POSITIONS.get(pdd_inset2_pos_var.get(), [0.53, 0.04, 0.45, 0.55]))
        try: _z2_bounds[0] = float(pdd_inset2_x0.get())
        except ValueError: pass
        try: _z2_bounds[1] = float(pdd_inset2_y0.get())
        except ValueError: pass
        try: _z2_bounds[2] = float(pdd_inset2_w.get())
        except ValueError: pass
        try: _z2_bounds[3] = float(pdd_inset2_h.get())
        except ValueError: pass
        ins_z2 = ax_z.inset_axes(_z2_bounds)
        _populate_inset(ins_z2, lambda a: a == 'Z')
        ins_z2.set_xlim([_zix1b, _zix2b])
        try:
            ins_z2.set_ylim([float(pdd_inset2_ymin.get()), float(pdd_inset2_ymax.get())])
        except ValueError:
            pass
        ins_z2.tick_params(labelsize=max(8, _fs - 4))
        _, _z2_connectors = ax_z.indicate_inset_zoom(ins_z2, edgecolor='black')
        for _c in _z2_connectors:
            _c.set_visible(True)

    try:
        _pix1c, _pix2c = float(profile_inset3_xmin.get()), float(profile_inset3_xmax.get())
    except ValueError:
        _pix1c = _pix2c = None

    try:
        _zix1c, _zix2c = float(pdd_inset3_xmin.get()), float(pdd_inset3_xmax.get())
    except ValueError:
        _zix1c = _zix2c = None

    if fig_profiles and show_profile_inset3_var.get() and _pix1c is not None and _pix1c < _pix2c:
        _p3_bounds = list(_INSET_POSITIONS.get(profile_inset3_pos_var.get(), [0.27, 0.37, 0.45, 0.58]))
        try: _p3_bounds[0] = float(profile_inset3_x0.get())
        except ValueError: pass
        try: _p3_bounds[1] = float(profile_inset3_y0.get())
        except ValueError: pass
        try: _p3_bounds[2] = float(profile_inset3_w.get())
        except ValueError: pass
        try: _p3_bounds[3] = float(profile_inset3_h.get())
        except ValueError: pass
        ins_p3 = ax_profiles.inset_axes(_p3_bounds)
        _populate_inset(ins_p3, lambda a: a in ('X', 'Y', 'XY', 'YX'))
        ins_p3.set_xlim([_pix1c, _pix2c])
        try:
            ins_p3.set_ylim([float(profile_inset3_ymin.get()), float(profile_inset3_ymax.get())])
        except ValueError:
            pass
        ins_p3.tick_params(labelsize=max(8, _fs - 4))
        _, _p3_connectors = ax_profiles.indicate_inset_zoom(ins_p3, edgecolor='black')
        for _c in _p3_connectors:
            _c.set_visible(True)

    if fig_z and show_pdd_inset3_var.get() and _zix1c is not None and _zix1c < _zix2c:
        _z3_bounds = list(_INSET_POSITIONS.get(pdd_inset3_pos_var.get(), [0.27, 0.04, 0.45, 0.55]))
        try: _z3_bounds[0] = float(pdd_inset3_x0.get())
        except ValueError: pass
        try: _z3_bounds[1] = float(pdd_inset3_y0.get())
        except ValueError: pass
        try: _z3_bounds[2] = float(pdd_inset3_w.get())
        except ValueError: pass
        try: _z3_bounds[3] = float(pdd_inset3_h.get())
        except ValueError: pass
        ins_z3 = ax_z.inset_axes(_z3_bounds)
        _populate_inset(ins_z3, lambda a: a == 'Z')
        ins_z3.set_xlim([_zix1c, _zix2c])
        try:
            ins_z3.set_ylim([float(pdd_inset3_ymin.get()), float(pdd_inset3_ymax.get())])
        except ValueError:
            pass
        ins_z3.tick_params(labelsize=max(8, _fs - 4))
        _, _z3_connectors = ax_z.indicate_inset_zoom(ins_z3, edgecolor='black')
        for _c in _z3_connectors:
            _c.set_visible(True)

    # Build proxy patch handles for bands
    _band_patch    = _MPatch(facecolor='0.85', edgecolor='black', linestyle='--', linewidth=1)
    _sd_patch      = _MPatch(facecolor='0.75', edgecolor='black', linestyle='--', linewidth=1)
    _ci_mean_patch = _MPatch(facecolor='0.60', edgecolor='black', linestyle='--', linewidth=1)
    try:
        _k_leg = float(sd_mult_entry.get())
    except ValueError:
        _k_leg = 2.326
    try:
        _ci_pct_leg = float(ci_entry.get())
    except ValueError:
        _ci_pct_leg = 98.0
    _cov_pct_leg = float((2 * _norm.cdf(_k_leg) - 1) * 100)
    dose_tolerance_leg, pos_tolerance_leg = get_clinical_criteria()

    # Finalize the plots
    if fig_profiles:
        h, l = ax_profiles.get_legend_handles_labels()
        if show_band_var.get() and avg_dose_global:
            h.append(_band_patch); l.append(f'({dose_tolerance_leg}%/{pos_tolerance_leg}mm) Envelope')
        if show_sd_var.get() and std_dose_global:
            h.append(_sd_patch);      l.append(f'±{_k_leg:.3f}·SD (Prediction, {_cov_pct_leg:.0f}%)')
        if show_ci_mean_var.get() and std_dose_global:
            h.append(_ci_mean_patch); l.append(f'{_ci_pct_leg:.0f}% CI on mean (t-dist)')
        ax_profiles.legend(h, l, loc='upper right', fontsize=_fs, markerscale=1.5)
        fig_profiles.show()

    if fig_z:
        h, l = ax_z.get_legend_handles_labels()
        if show_band_var.get() and avg_dose_global:
            h.append(_band_patch); l.append(f'({dose_tolerance_leg}%/{pos_tolerance_leg}mm) Envelope')
        if show_sd_var.get() and std_dose_global:
            h.append(_sd_patch);      l.append(f'±{_k_leg:.3f}·SD (Prediction, {_cov_pct_leg:.0f}%)')
        if show_ci_mean_var.get() and std_dose_global:
            h.append(_ci_mean_patch); l.append(f'{_ci_pct_leg:.0f}% CI on mean (t-dist)')
        ax_z.legend(h, l, loc='upper right', fontsize=_fs, markerscale=1.5)
        fig_z.show()
    print("Data plotting completed!")



    
def _modified_zscore(values):
    """Iglewicz-Hoaglin modified Z-score. Returns array of M scores; |M| > 3.5 flags an outlier."""
    vals = np.array(values, dtype=float)
    med  = np.median(vals)
    mad  = np.median(np.abs(vals - med))
    if mad == 0:
        std = np.std(vals)
        return np.zeros_like(vals) if std == 0 else (vals - np.mean(vals)) / std
    return 0.6745 * (vals - med) / mad


def _build_report_lines(dose_tolerance, pos_tolerance, sd_mult):
    ci_pct      = float((2 * _norm.cdf(sd_mult) - 1) * 100)
    pos_tol_cm  = pos_tolerance / 10.0
    tol_hdr     = f"Tol ({dose_tolerance}%/{pos_tolerance}mm)"
    ci_hdr      = f"Pred ({ci_pct:.0f}%)"
    tol_w       = max(14, len(tol_hdr))
    ci_w        = max(12, len(ci_hdr))
    _total_w    = 2 + 32 + 2 + 6 + 2 + 9 + 2 + 9 + 2 + 9 + 2 + 9 + 2 + 10 + 2 + tol_w + 2 + ci_w
    sep         = f"  {'-'*32}  {'-'*6}  {'-'*9}  {'-'*9}  {'-'*9}  {'-'*9}  {'-'*10}  {'-'*tol_w}  {'-'*ci_w}"

    lines = []
    lines.append("=" * _total_w)
    lines.append("BEAM DATA SPREAD REPORT")
    lines.append(f"Tolerance: {dose_tolerance}% / {pos_tolerance} mm   |   "
                 f"Prediction band: ±{sd_mult:.4f}·SD ({ci_pct:.1f}% coverage)   |   "
                 f"CI on mean: {ci_pct:.1f}% (t-dist, n varies per combination)")
    lines.append(f"† Max dev† = max deviation at points where DTA > {pos_tolerance} mm "
                 f"(search cap ±{pos_tolerance/10+0.1:.2f} cm); "
                 f"if all points are within {pos_tolerance} mm DTA, reports deviation at the point of maximum DTA")
    lines.append("=" * _total_w)

    # accumulators for summary: {sheet: [sq_devs, max_devs, dta_max_devs, n_in_tol, n_in_ci, n_tot]}
    summary = {}

    for combination in sorted(avg_dose_global.keys(), key=lambda c: (c[0], c[1], c[2] or 0)):
        fs, axis, depth = combination
        avg_dose  = avg_dose_global[combination]
        fixed_pos = fixed_pos_global[combination]
        std       = std_dose_global.get(combination, np.zeros_like(avg_dose))
        sheets    = _raw_stacks_global.get(combination, {})

        if axis == 'Z':
            lines.append(f"\nFS={fs}, Axis={axis}")
        else:
            lines.append(f"\nFS={fs}, Axis={axis}, Depth={depth} cm")

        valid_mask = ~np.isnan(avg_dose)
        mean_sd = float(np.nanmean(std[valid_mask])) if valid_mask.any() else float('nan')
        max_sd  = float(np.nanmax(std[valid_mask]))  if valid_mask.any() else float('nan')
        lines.append(f"  Spread across datasets:  mean SD = {mean_sd:.3f}%  |  max SD = {max_sd:.3f}%")
        n_combo = len(sheets)
        if n_combo >= 2 and valid_mask.any():
            _t_ci      = float(_t.ppf((1 + ci_pct / 100) / 2, df=n_combo - 1))
            _sem_ci    = std / np.sqrt(n_combo)
            ci_hw_mean = float(np.nanmean(_t_ci * _sem_ci[valid_mask]))
            ci_hw_max  = float(np.nanmax(_t_ci * _sem_ci[valid_mask]))
            lines.append(f"  CI on mean ({ci_pct:.0f}%, t={_t_ci:.3f}, n={n_combo}):  "
                         f"mean ±{ci_hw_mean:.3f}%  |  max ±{ci_hw_max:.3f}%")

        if sheets:
            upper_env, lower_env = _compute_envelope(avg_dose, fixed_pos, dose_tolerance, pos_tolerance, axis)
            hdr = (f"  {'Sheet':<32}  {'n':>6}  {'Mean dev':>9}  {'RMS dev':>9}  "
                   f"{'P95 dev':>9}  {'Max dev':>9}  {'Max dev†':>10}  {tol_hdr:>{tol_w}}  {ci_hdr:>{ci_w}}")
            lines.append(hdr)
            lines.append(sep)
            comb_sq, comb_mx, comb_dta_mx, comb_tol, comb_ci, comb_n, comb_diffs = [], [], [], 0, 0, 0, []
            for sn in sorted(sheets.keys()):
                sheet_dose = sheets[sn]
                valid = ~np.isnan(avg_dose) & ~np.isnan(sheet_dose)
                if not valid.any():
                    continue
                diff   = sheet_dose[valid] - avg_dose[valid]
                mean_d = float(np.mean(diff))
                rms    = float(np.sqrt(np.mean(diff ** 2)))
                p95_d  = float(np.percentile(np.abs(diff), 95))
                mxd    = float(np.max(np.abs(diff)))
                n_tot  = int(valid.sum())

                # DTA-filtered max dev
                dta_val    = _compute_dta(sheet_dose, avg_dose, fixed_pos, pos_tol_cm)
                dta_filter = dta_val > pos_tol_cm
                combo_filt = valid & dta_filter
                if combo_filt.any():
                    mxd_dta = float(np.max(np.abs((sheet_dose - avg_dose)[combo_filt])))
                else:
                    # All points pass DTA — report deviation at the point of highest DTA
                    masked_dta = np.where(valid, dta_val, -np.inf)
                    worst_idx  = int(np.argmax(masked_dta))
                    mxd_dta    = float(np.abs(sheet_dose[worst_idx] - avg_dose[worst_idx]))

                n_in_tol = int(np.sum(
                    (sheet_dose[valid] >= lower_env[valid]) &
                    (sheet_dose[valid] <= upper_env[valid])))
                pct_tol = 100.0 * n_in_tol / n_tot if n_tot > 0 else float('nan')

                ci_upper = avg_dose[valid] + sd_mult * std[valid]
                ci_lower = avg_dose[valid] - sd_mult * std[valid]
                n_in_ci  = int(np.sum((sheet_dose[valid] >= ci_lower) & (sheet_dose[valid] <= ci_upper)))
                pct_ci   = 100.0 * n_in_ci / n_tot if n_tot > 0 else float('nan')

                lines.append(f"  {sn:<32}  {n_tot:>6}  {mean_d:>+8.3f}%  {rms:>8.3f}%  "
                             f"{p95_d:>8.3f}%  {mxd:>8.3f}%  "
                             f"{mxd_dta:>9.3f}%  {pct_tol:>{tol_w-1}.1f}%  {pct_ci:>{ci_w-1}.1f}%")

                # accumulate combination ALL
                comb_sq.extend(diff ** 2)
                comb_diffs.extend(diff)
                comb_mx.append(mxd)
                comb_dta_mx.append(mxd_dta)
                comb_tol += n_in_tol
                comb_ci  += n_in_ci
                comb_n   += n_tot

                # accumulate summary
                acc = summary.setdefault(sn, [[], [], [], 0, 0, 0, []])
                acc[0].extend(diff ** 2)
                acc[1].append(mxd)
                acc[2].append(mxd_dta)
                acc[3] += n_in_tol
                acc[4] += n_in_ci
                acc[5] += n_tot
                acc[6].extend(diff)

            # ALL row for this combination
            if comb_n > 0:
                mean_d_c  = float(np.mean(comb_diffs))
                rms_c     = float(np.sqrt(np.mean(comb_sq)))
                p95_d_c   = float(np.percentile(np.abs(comb_diffs), 95))
                mxd_c     = float(np.max(comb_mx))
                mxd_dta_c = float(np.max(comb_dta_mx)) if comb_dta_mx else float('nan')
                pct_tol_c = 100.0 * comb_tol / comb_n
                pct_ci_c  = 100.0 * comb_ci  / comb_n
                lines.append(sep)
                lines.append(f"  {'ALL':<32}  {comb_n:>6}  {mean_d_c:>+8.3f}%  {rms_c:>8.3f}%  "
                              f"{p95_d_c:>8.3f}%  {mxd_c:>8.3f}%  "
                              f"{mxd_dta_c:>9.3f}%  {pct_tol_c:>{tol_w-1}.1f}%  {pct_ci_c:>{ci_w-1}.1f}%")

        lines.append("")

    # --- Overall summary ---
    if summary:
        lines.append("=" * _total_w)
        lines.append("SUMMARY (all combinations)")
        hdr = (f"  {'Sheet':<32}  {'n':>6}  {'Mean dev':>9}  {'RMS dev':>9}  "
               f"{'P95 dev':>9}  {'Max dev':>9}  {'Max dev†':>10}  {tol_hdr:>{tol_w}}  {ci_hdr:>{ci_w}}")
        lines.append(hdr)
        lines.append(sep)
        all_sq, all_mx, all_dta_mx, all_tol, all_ci, all_n, all_diffs = [], [], [], 0, 0, 0, []
        fleet_metrics = {}   # {sheet_name: (mean_d, rms, p95)}
        for sn in sorted(summary.keys()):
            sq_devs, max_devs, dta_max_devs, n_in_tol, n_in_ci, n_tot, diffs_all = summary[sn]
            mean_d_all = float(np.mean(diffs_all))                              if diffs_all else float('nan')
            rms_all    = float(np.sqrt(np.mean(sq_devs)))                       if sq_devs   else float('nan')
            p95_d_all  = float(np.percentile(np.abs(diffs_all), 95))            if diffs_all else float('nan')
            mxd_all    = float(np.max(max_devs))                                if max_devs  else float('nan')
            mxd_dta_a  = float(np.max(dta_max_devs))                            if dta_max_devs else float('nan')
            pct_tol    = 100.0 * n_in_tol / n_tot if n_tot > 0 else float('nan')
            pct_ci     = 100.0 * n_in_ci  / n_tot if n_tot > 0 else float('nan')
            lines.append(f"  {sn:<32}  {n_tot:>6}  {mean_d_all:>+8.3f}%  {rms_all:>8.3f}%  "
                         f"{p95_d_all:>8.3f}%  {mxd_all:>8.3f}%  "
                         f"{mxd_dta_a:>9.3f}%  {pct_tol:>{tol_w-1}.1f}%  {pct_ci:>{ci_w-1}.1f}%")
            fleet_metrics[sn] = (mean_d_all, rms_all, p95_d_all)
            all_sq.extend(sq_devs)
            all_diffs.extend(diffs_all)
            all_mx.append(mxd_all)
            all_dta_mx.append(mxd_dta_a)
            all_tol += n_in_tol
            all_ci  += n_in_ci
            all_n   += n_tot

        mean_d_tot  = float(np.mean(all_diffs))                   if all_diffs   else float('nan')
        rms_tot     = float(np.sqrt(np.mean(all_sq)))             if all_sq      else float('nan')
        p95_d_tot   = float(np.percentile(np.abs(all_diffs), 95)) if all_diffs   else float('nan')
        mxd_tot     = float(np.max(all_mx))                       if all_mx      else float('nan')
        mxd_dta_tot = float(np.max(all_dta_mx))                   if all_dta_mx  else float('nan')
        pct_tol_tot = 100.0 * all_tol / all_n if all_n > 0 else float('nan')
        pct_ci_tot  = 100.0 * all_ci  / all_n if all_n > 0 else float('nan')
        lines.append(sep)
        lines.append(f"  {'ALL':<32}  {all_n:>6}  {mean_d_tot:>+8.3f}%  {rms_tot:>8.3f}%  "
                     f"{p95_d_tot:>8.3f}%  {mxd_tot:>8.3f}%  "
                     f"{mxd_dta_tot:>9.3f}%  {pct_tol_tot:>{tol_w-1}.1f}%  {pct_ci_tot:>{ci_w-1}.1f}%")
        lines.append("")

        # --- Statistical outlier FLAGS (modified Z-score, Iglewicz-Hoaglin threshold 3.5) ---
        if len(fleet_metrics) >= 3:
            sheet_names = sorted(fleet_metrics.keys())
            means = [fleet_metrics[s][0] for s in sheet_names]
            rmss  = [fleet_metrics[s][1] for s in sheet_names]
            p95s  = [fleet_metrics[s][2] for s in sheet_names]
            z_mean = _modified_zscore(means)
            z_rms  = _modified_zscore(rmss)
            z_p95  = _modified_zscore(p95s)
            _THRESH = 3.5
            flags = {}
            for i, sn in enumerate(sheet_names):
                triggers = []
                if abs(z_mean[i]) > _THRESH:
                    triggers.append(f"Mean dev (M={z_mean[i]:+.1f})")
                if abs(z_rms[i]) > _THRESH:
                    triggers.append(f"RMS dev (M={z_rms[i]:+.1f})")
                if abs(z_p95[i]) > _THRESH:
                    triggers.append(f"P95 dev (M={z_p95[i]:+.1f})")
                if triggers:
                    flags[sn] = triggers
            lines.append("=" * _total_w)
            lines.append("STATISTICAL OUTLIER FLAGS  (modified Z-score threshold: |M| > 3.5)")
            lines.append("-" * _total_w)
            if flags:
                for sn, triggers in sorted(flags.items()):
                    lines.append(f"  *** {sn:<32}  {', '.join(triggers)}")
            else:
                lines.append("  No statistical outliers detected.")
            lines.append("")

    lines.append("=" * _total_w)
    return lines


def _get_sd_mult():
    try:
        return float(sd_mult_entry.get())
    except (ValueError, tk.TclError):
        return 2.326


def print_report():
    if not avg_dose_global:
        print("No average data. Run Make Avg first.")
        return
    dose_tolerance, pos_tolerance = get_clinical_criteria()
    for line in _build_report_lines(dose_tolerance, pos_tolerance, _get_sd_mult()):
        print(line)


def save_report():
    if not avg_dose_global:
        print("No average data. Run Make Avg first.")
        return
    dose_tolerance, pos_tolerance = get_clinical_criteria()
    out_path = filedialog.asksaveasfilename(
        defaultextension=".txt",
        filetypes=[("Text Files", "*.txt")],
        title="Save Report As")
    if not out_path:
        return
    with open(out_path, 'w') as f:
        f.write('\n'.join(_build_report_lines(dose_tolerance, pos_tolerance, _get_sd_mult())))
    print(f"Report saved to {out_path}")


def make_avg():
    dose_tolerance, pos_tolerance = get_clinical_criteria()
    global avg_dose_global, fixed_pos_global, std_dose_global, _raw_stacks_global, fsl, scl, dl

    selected_fsl_indices = fsl_listbox.curselection()
    selected_scl_indices = scl_listbox.curselection()
    selected_depth_indices = depth_listbox.curselection()

    fsl = [float(fsl_listbox.get(i)) for i in selected_fsl_indices]
    scl = [scl_listbox.get(i) for i in selected_scl_indices]
    dl = [float(depth_listbox.get(i)) for i in selected_depth_indices]

    if not fsl or not scl:
        print("Field size or scan types not selected!")
        return

    selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
    if not selected_sheets:
        print("No sheets selected!")
        return

    avg_sheets = (
        [s for s in selected_sheets if "TPS" not in s.upper()]
        if exclude_tps_var.get() else selected_sheets
    )
    if not avg_sheets:
        print("No non-TPS sheets selected for averaging!")
        return

    # Minimum number of datasets required at each position to include in average
    min_n_str = min_datasets_var.get()
    min_n = len(avg_sheets) if min_n_str == 'All' else max(1, int(min_n_str))

    fixed_pos_z = np.arange(0, 30.05, 0.05)
    max_field_size = max([int(fs) for fs in fsl]) if fsl else 0
    fixed_pos_profiles = np.arange(-max_field_size / 2 - 20, max_field_size / 2 + 20, 0.05)

    resampled_doses_by_combination = {}   # {combination: {sheet: array}}

    for sheet in avg_sheets:
        df = excel_data[sheet].sort_values(by=['FS', 'Axis', 'Pos'])
        for fs in fsl:
            for axis in scl:
                if axis in ['X', 'Y', 'XY', 'YX']:
                    for depth in dl:
                        subset = df[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth)]
                        if not subset.empty:
                            pos  = subset['Pos'].values
                            dose = subset['Dose'].values
                            pos, dose = _preprocess(pos, dose, axis, sheet)
                            interpolator = interp1d(pos, dose, bounds_error=False, fill_value=np.nan)
                            resampled_doses_by_combination.setdefault((fs, axis, depth), {})[sheet] = \
                                interpolator(fixed_pos_profiles)
                elif axis == 'Z':
                    subset = df[(df['FS'] == fs) & (df['Axis'] == axis)]
                    if not subset.empty:
                        pos  = subset['Pos'].values
                        dose = subset['Dose'].values
                        pos, dose = _preprocess(pos, dose, axis, sheet)
                        interpolator = interp1d(pos, dose, bounds_error=False, fill_value=np.nan)
                        resampled_doses_by_combination.setdefault((fs, axis, None), {})[sheet] = \
                            interpolator(fixed_pos_z)

    avg_dose_global  = {}
    fixed_pos_global = {}
    std_dose_global  = {}
    _raw_stacks_global = {}

    for combination, sheet_doses in resampled_doses_by_combination.items():
        fs, axis, depth = combination
        fixed_pos = fixed_pos_profiles if axis in ['X', 'Y', 'XY', 'YX'] else fixed_pos_z

        sheet_names = list(sheet_doses.keys())
        stack  = np.array(list(sheet_doses.values()))     # (n_sheets, n_pos)
        counts = np.sum(~np.isnan(stack), axis=0)
        with np.errstate(all='ignore'):
            avg = np.nanmean(stack, axis=0)
            std = np.nanstd(stack, axis=0, ddof=1) if stack.shape[0] > 1 else np.zeros_like(avg)
        avg[counts < min_n] = np.nan

        valid = ~np.isnan(avg)
        if valid.any():
            start = int(np.argmax(valid))
            end   = int(len(valid) - np.argmax(valid[::-1]))
            fixed_pos_global[combination]  = fixed_pos[start:end]
            avg_dose_global[combination]   = avg[start:end]
            std_dose_global[combination]   = std[start:end]
            _raw_stacks_global[combination] = {sn: stack[i, start:end]
                                               for i, sn in enumerate(sheet_names)}

    print(f"Average calculation completed! (min datasets = {min_n})")



# Initialize main window
root = tk.Tk()
root.title("Beam Data Averaging ToolV1.0")

# File selection section
file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

# Sheet selection section (Updated)
sheet_frame = ttk.Frame(root, padding="10")
sheet_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))  # Padding to separate from the file frame

# Listbox for selecting multiple sheets
ttk.Label(sheet_frame, text="Select Sheets:").grid(row=0, column=0, sticky="w")
sheet_listbox = tk.Listbox(sheet_frame, selectmode="multiple", width=40, height=6, exportselection=0)
sheet_listbox.grid(row=1, column=0, padx=5, pady=5)
sheet_listbox.bind("<<ListboxSelect>>", lambda _: (populate_fsl(), populate_scl(), populate_dl(), _update_highlight_combo()))

def toggle_sheets():
    if sheet_listbox.curselection() == tuple(range(sheet_listbox.size())):
        sheet_listbox.selection_clear(0, tk.END)
        sheet_toggle_btn.config(text="Select All")
    else:
        sheet_listbox.selection_set(0, tk.END)
        sheet_toggle_btn.config(text="Deselect All")

sheet_toggle_btn = ttk.Button(sheet_frame, text="Select All", command=toggle_sheets)
sheet_toggle_btn.grid(row=2, column=0, pady=2)

# Populating the listbox with sheet names from the Excel file (Done in `select_file()` function)


# Data Selection Frame
data_selection_frame = ttk.LabelFrame(root, text="Data Selection", padding="10")
data_selection_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))  # Adjust row number to fit below file and sheet sections

# Field Size Selection
ttk.Label(data_selection_frame, text="Select Field Sizes").grid(row=0, column=0, sticky="w")
fsl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=6, exportselection=0)
fsl_listbox.grid(row=1, column=0, padx=5, pady=5)

def toggle_fsl():
    if fsl_listbox.curselection() == tuple(range(fsl_listbox.size())):
        fsl_listbox.selection_clear(0, tk.END)
        fsl_toggle_btn.config(text="Select All")
    else:
        fsl_listbox.selection_set(0, tk.END)
        fsl_toggle_btn.config(text="Deselect All")

fsl_toggle_btn = ttk.Button(data_selection_frame, text="Select All", command=toggle_fsl)
fsl_toggle_btn.grid(row=2, column=0, pady=2)

# Scan Types Selection
ttk.Label(data_selection_frame, text="Select Scan Types").grid(row=0, column=1, sticky="w", padx=20)
scl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=6, exportselection=0)
scl_listbox.grid(row=1, column=1, padx=5, pady=5)

def toggle_scl():
    if scl_listbox.curselection() == tuple(range(scl_listbox.size())):
        scl_listbox.selection_clear(0, tk.END)
        scl_toggle_btn.config(text="Select All")
    else:
        scl_listbox.selection_set(0, tk.END)
        scl_toggle_btn.config(text="Deselect All")

scl_toggle_btn = ttk.Button(data_selection_frame, text="Select All", command=toggle_scl)
scl_toggle_btn.grid(row=2, column=1, pady=2)

# Depth Selection
ttk.Label(data_selection_frame, text="Select Depths").grid(row=0, column=2, sticky="w", padx=20)
depth_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=6, exportselection=0)
depth_listbox.grid(row=1, column=2, padx=5, pady=5)

def toggle_depth():
    if depth_listbox.curselection() == tuple(range(depth_listbox.size())):
        depth_listbox.selection_clear(0, tk.END)
        depth_toggle_btn.config(text="Select All")
    else:
        depth_listbox.selection_set(0, tk.END)
        depth_toggle_btn.config(text="Deselect All")

depth_toggle_btn = ttk.Button(data_selection_frame, text="Select All", command=toggle_depth)
depth_toggle_btn.grid(row=2, column=2, pady=2)


# Averaging threshold + tolerance controls in one frame
controls_frame = ttk.Frame(root, padding="5")
controls_frame.grid(row=5, column=0, sticky="w", padx=10)

ttk.Label(controls_frame, text="Min. datasets for avg:").grid(row=0, column=0, sticky="w", padx=(0, 5))
min_datasets_var = tk.StringVar(value="All")
min_ds_combo = ttk.Combobox(controls_frame, textvariable=min_datasets_var, width=6, state="readonly")
min_ds_combo['values'] = ['All'] + [str(i) for i in range(1, 21)]
min_ds_combo.grid(row=0, column=1, sticky="w", padx=(0, 20))

ttk.Label(controls_frame, text="Dose Tolerance [%]:").grid(row=0, column=2, sticky="w", padx=(0, 5))
dose_entry = ttk.Entry(controls_frame, width=6)
dose_entry.insert(0, "1")
dose_entry.grid(row=0, column=3, sticky="w", padx=(0, 20))

ttk.Label(controls_frame, text="Position Tolerance [mm]:").grid(row=0, column=4, sticky="w", padx=(0, 5))
pos_entry = ttk.Entry(controls_frame, width=6)
pos_entry.insert(0, "1")
pos_entry.grid(row=0, column=5, sticky="w", padx=(0, 20))

ttk.Label(controls_frame, text="Marker size:").grid(row=0, column=6, sticky="w", padx=(0, 5))
marker_size_entry = ttk.Entry(controls_frame, width=5)
marker_size_entry.insert(0, "3")
marker_size_entry.grid(row=0, column=7, sticky="w", padx=(0, 20))

ttk.Label(controls_frame, text="Font size:").grid(row=0, column=8, sticky="w", padx=(0, 5))
font_size_entry = ttk.Entry(controls_frame, width=5)
font_size_entry.insert(0, "16")
font_size_entry.grid(row=0, column=9, sticky="w", padx=(0, 20))

show_band_var = tk.BooleanVar(value=True)
ttk.Checkbutton(controls_frame, text="Show tolerance band", variable=show_band_var).grid(
    row=0, column=10, sticky="w")

ttk.Label(controls_frame, text="Profile x [cm]:").grid(row=1, column=0, sticky="w", padx=(0, 5))
profile_xmin_entry = ttk.Entry(controls_frame, width=6)
profile_xmin_entry.grid(row=1, column=1, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="to").grid(row=1, column=2, sticky="w")
profile_xmax_entry = ttk.Entry(controls_frame, width=6)
profile_xmax_entry.grid(row=1, column=3, sticky="w", padx=(0, 20))

ttk.Label(controls_frame, text="PDD x [cm]:").grid(row=1, column=4, sticky="w", padx=(0, 5))
z_xmin_entry = ttk.Entry(controls_frame, width=6)
z_xmin_entry.insert(0, "0")
z_xmin_entry.grid(row=1, column=5, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="to").grid(row=1, column=6, sticky="w")
z_xmax_entry = ttk.Entry(controls_frame, width=6)
z_xmax_entry.insert(0, "30")
z_xmax_entry.grid(row=1, column=7, sticky="w", padx=(0, 20))

show_sd_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="Show SD band", variable=show_sd_var).grid(
    row=1, column=8, sticky="w", padx=(0, 8))

ttk.Label(controls_frame, text="CI [%]:").grid(row=1, column=9, sticky="w", padx=(0, 3))
ci_entry = ttk.Entry(controls_frame, width=6)
ci_entry.insert(0, "98")
ci_entry.grid(row=1, column=10, sticky="w", padx=(0, 8))

ttk.Label(controls_frame, text="SD mult:").grid(row=1, column=11, sticky="w", padx=(0, 3))
sd_mult_entry = ttk.Entry(controls_frame, width=6)
sd_mult_entry.insert(0, "2.326")
sd_mult_entry.grid(row=1, column=12, sticky="w", padx=(0, 8))

def _ci_changed(*_):
    try:
        ci = float(ci_entry.get())
        z = float(_norm.ppf((1 + ci / 100) / 2))
        sd_mult_entry.delete(0, tk.END)
        sd_mult_entry.insert(0, f"{z:.4f}")
    except Exception:
        pass

def _sdmult_changed(*_):
    try:
        z = float(sd_mult_entry.get())
        ci = float((2 * _norm.cdf(z) - 1) * 100)
        ci_entry.delete(0, tk.END)
        ci_entry.insert(0, f"{ci:.2f}")
    except Exception:
        pass

ci_entry.bind("<FocusOut>", _ci_changed)
ci_entry.bind("<Return>", _ci_changed)
sd_mult_entry.bind("<FocusOut>", _sdmult_changed)
sd_mult_entry.bind("<Return>", _sdmult_changed)

show_ci_mean_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="Show CI on mean (t-dist)", variable=show_ci_mean_var).grid(
    row=1, column=13, sticky="w", padx=(8, 0))

_INSET_POSITIONS = {
    "Upper Right":  [0.53, 0.37, 0.45, 0.58],
    "Upper Left":   [0.02, 0.37, 0.45, 0.58],
    "Upper Center": [0.27, 0.37, 0.45, 0.58],
    "Lower Right":  [0.53, 0.04, 0.45, 0.55],
    "Lower Left":   [0.02, 0.04, 0.45, 0.55],
    "Lower Center": [0.27, 0.04, 0.45, 0.55],
}

_CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BeamDataAvgTool_config.json")

def _inset_dict(show_var, xmin, xmax, ymin, ymax, pos_var, x0, y0, w, h):
    return {"show": show_var.get(), "xmin": xmin.get(), "xmax": xmax.get(),
            "ymin": ymin.get(), "ymax": ymax.get(), "pos": pos_var.get(),
            "x0": x0.get(), "y0": y0.get(), "w": w.get(), "h": h.get()}

def _load_inset_dict(d, show_var, xmin, xmax, ymin, ymax, pos_var, x0, y0, w, h):
    def _set(entry, val):
        entry.delete(0, tk.END)
        if val: entry.insert(0, val)
    show_var.set(d.get("show", False))
    _set(xmin, d.get("xmin", "")); _set(xmax, d.get("xmax", ""))
    _set(ymin, d.get("ymin", "")); _set(ymax, d.get("ymax", ""))
    if d.get("pos") in _INSET_POSITIONS: pos_var.set(d["pos"])
    _set(x0, d.get("x0", "")); _set(y0, d.get("y0", ""))
    _set(w,  d.get("w",  "")); _set(h,  d.get("h",  ""))

def _save_inset_config():
    import json
    cfg = {
        "profile_inset":  _inset_dict(show_profile_inset_var,  profile_inset_xmin,  profile_inset_xmax,  profile_inset_ymin,  profile_inset_ymax,  profile_inset_pos_var,  profile_inset_x0,  profile_inset_y0,  profile_inset_w,  profile_inset_h),
        "pdd_inset":      _inset_dict(show_pdd_inset_var,      pdd_inset_xmin,      pdd_inset_xmax,      pdd_inset_ymin,      pdd_inset_ymax,      pdd_inset_pos_var,      pdd_inset_x0,      pdd_inset_y0,      pdd_inset_w,      pdd_inset_h),
        "profile_inset2": _inset_dict(show_profile_inset2_var, profile_inset2_xmin, profile_inset2_xmax, profile_inset2_ymin, profile_inset2_ymax, profile_inset2_pos_var, profile_inset2_x0, profile_inset2_y0, profile_inset2_w, profile_inset2_h),
        "pdd_inset2":     _inset_dict(show_pdd_inset2_var,     pdd_inset2_xmin,     pdd_inset2_xmax,     pdd_inset2_ymin,     pdd_inset2_ymax,     pdd_inset2_pos_var,     pdd_inset2_x0,     pdd_inset2_y0,     pdd_inset2_w,     pdd_inset2_h),
        "profile_inset3": _inset_dict(show_profile_inset3_var, profile_inset3_xmin, profile_inset3_xmax, profile_inset3_ymin, profile_inset3_ymax, profile_inset3_pos_var, profile_inset3_x0, profile_inset3_y0, profile_inset3_w, profile_inset3_h),
        "pdd_inset3":     _inset_dict(show_pdd_inset3_var,     pdd_inset3_xmin,     pdd_inset3_xmax,     pdd_inset3_ymin,     pdd_inset3_ymax,     pdd_inset3_pos_var,     pdd_inset3_x0,     pdd_inset3_y0,     pdd_inset3_w,     pdd_inset3_h),
    }
    try:
        with open(_CONFIG_PATH, "w") as f:
            json.dump(cfg, f, indent=2)
        print(f"Inset settings saved to {_CONFIG_PATH}")
    except Exception as e:
        print(f"Could not save inset settings: {e}")

def _load_inset_config():
    import json
    try:
        with open(_CONFIG_PATH) as f:
            cfg = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return
    _load_inset_dict(cfg.get("profile_inset",  {}), show_profile_inset_var,  profile_inset_xmin,  profile_inset_xmax,  profile_inset_ymin,  profile_inset_ymax,  profile_inset_pos_var,  profile_inset_x0,  profile_inset_y0,  profile_inset_w,  profile_inset_h)
    _load_inset_dict(cfg.get("pdd_inset",      {}), show_pdd_inset_var,      pdd_inset_xmin,      pdd_inset_xmax,      pdd_inset_ymin,      pdd_inset_ymax,      pdd_inset_pos_var,      pdd_inset_x0,      pdd_inset_y0,      pdd_inset_w,      pdd_inset_h)
    _load_inset_dict(cfg.get("profile_inset2", {}), show_profile_inset2_var, profile_inset2_xmin, profile_inset2_xmax, profile_inset2_ymin, profile_inset2_ymax, profile_inset2_pos_var, profile_inset2_x0, profile_inset2_y0, profile_inset2_w, profile_inset2_h)
    _load_inset_dict(cfg.get("pdd_inset2",     {}), show_pdd_inset2_var,     pdd_inset2_xmin,     pdd_inset2_xmax,     pdd_inset2_ymin,     pdd_inset2_ymax,     pdd_inset2_pos_var,     pdd_inset2_x0,     pdd_inset2_y0,     pdd_inset2_w,     pdd_inset2_h)
    _load_inset_dict(cfg.get("profile_inset3", {}), show_profile_inset3_var, profile_inset3_xmin, profile_inset3_xmax, profile_inset3_ymin, profile_inset3_ymax, profile_inset3_pos_var, profile_inset3_x0, profile_inset3_y0, profile_inset3_w, profile_inset3_h)
    _load_inset_dict(cfg.get("pdd_inset3",     {}), show_pdd_inset3_var,     pdd_inset3_xmin,     pdd_inset3_xmax,     pdd_inset3_ymin,     pdd_inset3_ymax,     pdd_inset3_pos_var,     pdd_inset3_x0,     pdd_inset3_y0,     pdd_inset3_w,     pdd_inset3_h)

show_profile_inset_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="Profile inset x [cm]:", variable=show_profile_inset_var).grid(row=3, column=0, sticky="w", padx=(0, 4))
profile_inset_xmin = ttk.Entry(controls_frame, width=6)
profile_inset_xmin.grid(row=3, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=3, column=2, sticky="w")
profile_inset_xmax = ttk.Entry(controls_frame, width=6)
profile_inset_xmax.grid(row=3, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=3, column=4, sticky="w", padx=(0, 4))
profile_inset_ymin = ttk.Entry(controls_frame, width=6)
profile_inset_ymin.grid(row=3, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=3, column=6, sticky="w")
profile_inset_ymax = ttk.Entry(controls_frame, width=6)
profile_inset_ymax.grid(row=3, column=7, sticky="w", padx=(0, 10))
profile_inset_pos_var = tk.StringVar(value="Upper Right")
ttk.Combobox(controls_frame, textvariable=profile_inset_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=3, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=3, column=9, sticky="w", padx=(0, 2))
profile_inset_x0 = ttk.Entry(controls_frame, width=5)
profile_inset_x0.grid(row=3, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=3, column=11, sticky="w", padx=(0, 2))
profile_inset_y0 = ttk.Entry(controls_frame, width=5)
profile_inset_y0.grid(row=3, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=3, column=13, sticky="w", padx=(0, 2))
profile_inset_w = ttk.Entry(controls_frame, width=5)
profile_inset_w.grid(row=3, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=3, column=15, sticky="w", padx=(0, 2))
profile_inset_h = ttk.Entry(controls_frame, width=5)
profile_inset_h.grid(row=3, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=3, column=17, padx=(4, 0))

show_pdd_inset_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="PDD inset x [cm]:", variable=show_pdd_inset_var).grid(row=4, column=0, sticky="w", padx=(0, 4))
pdd_inset_xmin = ttk.Entry(controls_frame, width=6)
pdd_inset_xmin.grid(row=4, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=4, column=2, sticky="w")
pdd_inset_xmax = ttk.Entry(controls_frame, width=6)
pdd_inset_xmax.grid(row=4, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=4, column=4, sticky="w", padx=(0, 4))
pdd_inset_ymin = ttk.Entry(controls_frame, width=6)
pdd_inset_ymin.grid(row=4, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=4, column=6, sticky="w")
pdd_inset_ymax = ttk.Entry(controls_frame, width=6)
pdd_inset_ymax.grid(row=4, column=7, sticky="w", padx=(0, 10))
pdd_inset_pos_var = tk.StringVar(value="Upper Right")
ttk.Combobox(controls_frame, textvariable=pdd_inset_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=4, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=4, column=9, sticky="w", padx=(0, 2))
pdd_inset_x0 = ttk.Entry(controls_frame, width=5)
pdd_inset_x0.grid(row=4, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=4, column=11, sticky="w", padx=(0, 2))
pdd_inset_y0 = ttk.Entry(controls_frame, width=5)
pdd_inset_y0.grid(row=4, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=4, column=13, sticky="w", padx=(0, 2))
pdd_inset_w = ttk.Entry(controls_frame, width=5)
pdd_inset_w.grid(row=4, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=4, column=15, sticky="w", padx=(0, 2))
pdd_inset_h = ttk.Entry(controls_frame, width=5)
pdd_inset_h.grid(row=4, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=4, column=17, padx=(4, 0))

show_profile_inset2_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="Profile inset 2 x [cm]:", variable=show_profile_inset2_var).grid(row=5, column=0, sticky="w", padx=(0, 4))
profile_inset2_xmin = ttk.Entry(controls_frame, width=6)
profile_inset2_xmin.grid(row=5, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=5, column=2, sticky="w")
profile_inset2_xmax = ttk.Entry(controls_frame, width=6)
profile_inset2_xmax.grid(row=5, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=5, column=4, sticky="w", padx=(0, 4))
profile_inset2_ymin = ttk.Entry(controls_frame, width=6)
profile_inset2_ymin.grid(row=5, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=5, column=6, sticky="w")
profile_inset2_ymax = ttk.Entry(controls_frame, width=6)
profile_inset2_ymax.grid(row=5, column=7, sticky="w", padx=(0, 10))
profile_inset2_pos_var = tk.StringVar(value="Upper Left")
ttk.Combobox(controls_frame, textvariable=profile_inset2_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=5, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=5, column=9, sticky="w", padx=(0, 2))
profile_inset2_x0 = ttk.Entry(controls_frame, width=5)
profile_inset2_x0.grid(row=5, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=5, column=11, sticky="w", padx=(0, 2))
profile_inset2_y0 = ttk.Entry(controls_frame, width=5)
profile_inset2_y0.grid(row=5, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=5, column=13, sticky="w", padx=(0, 2))
profile_inset2_w = ttk.Entry(controls_frame, width=5)
profile_inset2_w.grid(row=5, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=5, column=15, sticky="w", padx=(0, 2))
profile_inset2_h = ttk.Entry(controls_frame, width=5)
profile_inset2_h.grid(row=5, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=5, column=17, padx=(4, 0))

show_pdd_inset2_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="PDD inset 2 x [cm]:", variable=show_pdd_inset2_var).grid(row=6, column=0, sticky="w", padx=(0, 4))
pdd_inset2_xmin = ttk.Entry(controls_frame, width=6)
pdd_inset2_xmin.grid(row=6, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=6, column=2, sticky="w")
pdd_inset2_xmax = ttk.Entry(controls_frame, width=6)
pdd_inset2_xmax.grid(row=6, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=6, column=4, sticky="w", padx=(0, 4))
pdd_inset2_ymin = ttk.Entry(controls_frame, width=6)
pdd_inset2_ymin.grid(row=6, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=6, column=6, sticky="w")
pdd_inset2_ymax = ttk.Entry(controls_frame, width=6)
pdd_inset2_ymax.grid(row=6, column=7, sticky="w", padx=(0, 10))
pdd_inset2_pos_var = tk.StringVar(value="Lower Right")
ttk.Combobox(controls_frame, textvariable=pdd_inset2_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=6, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=6, column=9, sticky="w", padx=(0, 2))
pdd_inset2_x0 = ttk.Entry(controls_frame, width=5)
pdd_inset2_x0.grid(row=6, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=6, column=11, sticky="w", padx=(0, 2))
pdd_inset2_y0 = ttk.Entry(controls_frame, width=5)
pdd_inset2_y0.grid(row=6, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=6, column=13, sticky="w", padx=(0, 2))
pdd_inset2_w = ttk.Entry(controls_frame, width=5)
pdd_inset2_w.grid(row=6, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=6, column=15, sticky="w", padx=(0, 2))
pdd_inset2_h = ttk.Entry(controls_frame, width=5)
pdd_inset2_h.grid(row=6, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=6, column=17, padx=(4, 0))

show_profile_inset3_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="Profile inset 3 x [cm]:", variable=show_profile_inset3_var).grid(row=7, column=0, sticky="w", padx=(0, 4))
profile_inset3_xmin = ttk.Entry(controls_frame, width=6)
profile_inset3_xmin.grid(row=7, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=7, column=2, sticky="w")
profile_inset3_xmax = ttk.Entry(controls_frame, width=6)
profile_inset3_xmax.grid(row=7, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=7, column=4, sticky="w", padx=(0, 4))
profile_inset3_ymin = ttk.Entry(controls_frame, width=6)
profile_inset3_ymin.grid(row=7, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=7, column=6, sticky="w")
profile_inset3_ymax = ttk.Entry(controls_frame, width=6)
profile_inset3_ymax.grid(row=7, column=7, sticky="w", padx=(0, 10))
profile_inset3_pos_var = tk.StringVar(value="Upper Center")
ttk.Combobox(controls_frame, textvariable=profile_inset3_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=7, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=7, column=9, sticky="w", padx=(0, 2))
profile_inset3_x0 = ttk.Entry(controls_frame, width=5)
profile_inset3_x0.grid(row=7, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=7, column=11, sticky="w", padx=(0, 2))
profile_inset3_y0 = ttk.Entry(controls_frame, width=5)
profile_inset3_y0.grid(row=7, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=7, column=13, sticky="w", padx=(0, 2))
profile_inset3_w = ttk.Entry(controls_frame, width=5)
profile_inset3_w.grid(row=7, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=7, column=15, sticky="w", padx=(0, 2))
profile_inset3_h = ttk.Entry(controls_frame, width=5)
profile_inset3_h.grid(row=7, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=7, column=17, padx=(4, 0))

show_pdd_inset3_var = tk.BooleanVar(value=False)
ttk.Checkbutton(controls_frame, text="PDD inset 3 x [cm]:", variable=show_pdd_inset3_var).grid(row=8, column=0, sticky="w", padx=(0, 4))
pdd_inset3_xmin = ttk.Entry(controls_frame, width=6)
pdd_inset3_xmin.grid(row=8, column=1, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=8, column=2, sticky="w")
pdd_inset3_xmax = ttk.Entry(controls_frame, width=6)
pdd_inset3_xmax.grid(row=8, column=3, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="y [%]:").grid(row=8, column=4, sticky="w", padx=(0, 4))
pdd_inset3_ymin = ttk.Entry(controls_frame, width=6)
pdd_inset3_ymin.grid(row=8, column=5, sticky="w", padx=(0, 3))
ttk.Label(controls_frame, text="to").grid(row=8, column=6, sticky="w")
pdd_inset3_ymax = ttk.Entry(controls_frame, width=6)
pdd_inset3_ymax.grid(row=8, column=7, sticky="w", padx=(0, 10))
pdd_inset3_pos_var = tk.StringVar(value="Lower Center")
ttk.Combobox(controls_frame, textvariable=pdd_inset3_pos_var,
             values=list(_INSET_POSITIONS.keys()), width=12, state="readonly").grid(
    row=8, column=8, sticky="w", padx=(0, 10))
ttk.Label(controls_frame, text="pos x:").grid(row=8, column=9, sticky="w", padx=(0, 2))
pdd_inset3_x0 = ttk.Entry(controls_frame, width=5)
pdd_inset3_x0.grid(row=8, column=10, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="y:").grid(row=8, column=11, sticky="w", padx=(0, 2))
pdd_inset3_y0 = ttk.Entry(controls_frame, width=5)
pdd_inset3_y0.grid(row=8, column=12, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="w:").grid(row=8, column=13, sticky="w", padx=(0, 2))
pdd_inset3_w = ttk.Entry(controls_frame, width=5)
pdd_inset3_w.grid(row=8, column=14, sticky="w", padx=(0, 4))
ttk.Label(controls_frame, text="h:").grid(row=8, column=15, sticky="w", padx=(0, 2))
pdd_inset3_h = ttk.Entry(controls_frame, width=5)
pdd_inset3_h.grid(row=8, column=16, sticky="w", padx=(0, 4))
ttk.Button(controls_frame, text="Save", command=_save_inset_config).grid(row=8, column=17, padx=(4, 0))

ttk.Label(controls_frame, text="Highlight machine:").grid(row=2, column=0, sticky="w", padx=(0, 5))
highlight_var = tk.StringVar(value="None")
highlight_combo = ttk.Combobox(controls_frame, textvariable=highlight_var, width=28, state="readonly")
highlight_combo['values'] = ["None"]
highlight_combo.grid(row=2, column=1, columnspan=4, sticky="w", padx=(0, 20))

def _update_highlight_combo(*_):
    selected = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
    cur = highlight_var.get()
    highlight_combo['values'] = ["None"] + selected
    highlight_var.set(cur if cur in selected else "None")

# Pre-processing frame
preproc_frame = ttk.LabelFrame(root, text="Pre-processing", padding="5")
preproc_frame.grid(row=6, column=0, sticky="ew", padx=10, pady=(0, 6))

auto_center_var = tk.BooleanVar(value=True)
ttk.Checkbutton(preproc_frame, text="Auto-center profiles (FWHM)", variable=auto_center_var).grid(
    row=0, column=0, sticky="w", padx=(0, 16))

cax_norm_var = tk.BooleanVar(value=True)
ttk.Checkbutton(preproc_frame, text="Normalize profiles to CAX", variable=cax_norm_var).grid(
    row=0, column=1, sticky="w", padx=(0, 16))

_TPS_D10 = {"6X": 66.5, "10X": 73.6, "15X": 76.5, "6FFF": 63.4, "8FFF": 68.2, "10FFF": 70.1}

ttk.Label(preproc_frame, text="PDD norm:").grid(row=0, column=2, sticky="w", padx=(0, 4))
pdd_norm_var = tk.StringVar(value="dmax")
pdd_norm_combo = ttk.Combobox(preproc_frame, textvariable=pdd_norm_var, width=8, state="readonly")
pdd_norm_combo['values'] = ["dmax", "depth", "TPS d10"]
pdd_norm_combo.grid(row=0, column=3, sticky="w", padx=(0, 4))

ttk.Label(preproc_frame, text="Depth [cm]:").grid(row=0, column=4, sticky="w", padx=(0, 2))
pdd_norm_depth_entry = ttk.Entry(preproc_frame, width=6)
pdd_norm_depth_entry.insert(0, "10.0")
pdd_norm_depth_entry.grid(row=0, column=5, sticky="w", padx=(0, 16))

# TPS d10 widgets (shown only when "TPS d10" selected)
_tps_energy_lbl  = ttk.Label(preproc_frame, text="Energy:")
_tps_energy_var  = tk.StringVar(value="6X")
_tps_energy_combo = ttk.Combobox(preproc_frame, textvariable=_tps_energy_var,
                                  width=7, state="readonly",
                                  values=list(_TPS_D10.keys()))
_tps_val_lbl  = ttk.Label(preproc_frame, text="%dd(10):")
tps_d10_entry = ttk.Entry(preproc_frame, width=6)
tps_d10_entry.insert(0, str(_TPS_D10["6X"]))

def _on_tps_energy(*_):
    tps_d10_entry.delete(0, tk.END)
    tps_d10_entry.insert(0, str(_TPS_D10.get(_tps_energy_var.get(), "")))
_tps_energy_var.trace_add("write", _on_tps_energy)

def _on_pdd_norm_change(*_):
    mode = pdd_norm_var.get()
    if mode == "depth":
        _tps_energy_lbl.grid_remove(); _tps_energy_combo.grid_remove()
        _tps_val_lbl.grid_remove();    tps_d10_entry.grid_remove()
        pdd_norm_depth_entry.grid(row=0, column=5, sticky="w", padx=(0, 16))
    elif mode == "TPS d10":
        pdd_norm_depth_entry.grid_remove()
        _tps_energy_lbl.grid( row=0, column=5, sticky="w", padx=(0, 2))
        _tps_energy_combo.grid(row=0, column=6, sticky="w", padx=(0, 4))
        _tps_val_lbl.grid(    row=0, column=7, sticky="w", padx=(0, 2))
        tps_d10_entry.grid(   row=0, column=8, sticky="w", padx=(0, 16))
    else:  # dmax
        pdd_norm_depth_entry.grid_remove()
        _tps_energy_lbl.grid_remove(); _tps_energy_combo.grid_remove()
        _tps_val_lbl.grid_remove();    tps_d10_entry.grid_remove()

pdd_norm_combo.bind("<<ComboboxSelected>>", _on_pdd_norm_change)
_on_pdd_norm_change()  # set initial state

ttk.Button(preproc_frame, text="Manual Profile Shifts…", command=_open_manual_shifts).grid(
    row=0, column=9, sticky="w", padx=(0, 4))

exclude_tps_var = tk.BooleanVar(value=True)
ttk.Checkbutton(preproc_frame, text="Exclude TPS tabs from avg", variable=exclude_tps_var).grid(
    row=0, column=10, sticky="w", padx=(16, 0))

def clear_avg():
    global avg_dose_global, fixed_pos_global, std_dose_global, _raw_stacks_global
    avg_dose_global    = {}
    fixed_pos_global   = {}
    std_dose_global    = {}
    _raw_stacks_global = {}
    print("Average cleared.")

def open_help():
    help_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BeamDataAvgTool_Help.html")
    webbrowser.open(f"file:///{help_path.replace(os.sep, '/')}")

buttons_frame = ttk.Frame(root, padding="5")
buttons_frame.grid(row=8, column=0, sticky="w", padx=10, pady=5)
ttk.Button(buttons_frame, text="Plot Data",    command=plot_data).grid(   row=0, column=0, padx=4)
ttk.Button(buttons_frame, text="Make Avg",     command=make_avg).grid(    row=0, column=1, padx=4)
ttk.Button(buttons_frame, text="Save Avg",     command=save_average).grid(row=0, column=2, padx=4)
ttk.Button(buttons_frame, text="Clear Avg",    command=clear_avg).grid(   row=0, column=3, padx=4)
ttk.Button(buttons_frame, text="Print Report", command=print_report).grid(row=0, column=4, padx=4)
ttk.Button(buttons_frame, text="Save Report",  command=save_report).grid( row=0, column=5, padx=4)
ttk.Button(buttons_frame, text="Help",         command=open_help).grid(   row=0, column=6, padx=4)

# Restore saved inset settings (if any)
_load_inset_config()

# Start the GUI loop
root.mainloop()