# Script to make AVG PDDs from multiple data sets and output the avg and plot that plus a confidence interval
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
from scipy.interpolate import interp1d

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
            sheets = list(excel_data.keys())
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

def plot_data():
    global dl, fsl, scl, avg_dose_global, fixed_pos_global

    # Get selected field sizes, scan types, and depths
    selected_fsl_indices = fsl_listbox.curselection()
    selected_scl_indices = scl_listbox.curselection()
    selected_depth_indices = depth_listbox.curselection()

    fsl = [float(fsl_listbox.get(i)) for i in selected_fsl_indices]
    scl = [scl_listbox.get(i) for i in selected_scl_indices]
    dl = [float(depth_listbox.get(i)) for i in selected_depth_indices] if selected_depth_indices else []

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
        fig_profiles, ax_profiles = plt.subplots(figsize=(10, 5))
        ax_profiles.set_title("Profile Scans (X, Y, XY, YX)")
        ax_profiles.set_ylabel('Dose [%]')
        ax_profiles.set_xlabel('Off Axis Position [cm]')
        max_field_size = max(fsl) if fsl else 0
        ax_profiles.set_xlim([-max_field_size/2-20, max_field_size/2+20])
        ax_profiles.set_ylim([0, 110])

    # Set up plotting for Z profiles
    fig_z, ax_z = None, None
    if 'Z' in scl:
        fig_z = plt.figure(figsize=(10, 5))
        ax_z = fig_z.add_subplot(111)
        ax_z.set_title("Z Profiles")
        ax_z.set_ylabel('Percentage Depth Dose [%]')
        ax_z.set_xlabel('Position [cm]')
        ax_z.set_xlim([0, 30])
        ax_z.set_ylim([0, 110])

    # Define a color map for sheets
    color_map = plt.cm.get_cmap("tab20", len(selected_sheets))

    # Plot actual data from sheets
    for sheet_idx, sheet in enumerate(selected_sheets):
        df = excel_data[sheet].sort_values(by=['FS', 'Axis', 'Pos'])
        color = color_map(sheet_idx)
        marker = markers[sheet_idx % len(markers)]
        sheet_labeled = False  # To label the sheet only once in the legend
        for fs in fsl:
            for axis in scl:
                if axis in ['X', 'Y', 'XY', 'YX']:
                    for depth in dl:
                        if not df[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth)].empty:
                            pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Pos'].values
                            dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Dose'].values
                            dose = dose / dose.max() * 100  # Normalize dose
                            if not sheet_labeled:
                                ax_profiles.plot(pos, dose, linestyle='None', marker=marker, color=color, ms=5, alpha=0.5, label=sheet)
                                sheet_labeled = True  # Label only once
                            else:
                                ax_profiles.plot(pos, dose, linestyle='None', marker=marker, color=color, ms=5, alpha=0.5, label=None)  # No label after the first

                elif axis == 'Z':
                    if not df[(df['FS'] == fs) & (df['Axis'] == axis)].empty:
                        pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Pos'].values
                        dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Dose'].values
                        dose = dose / dose.max() * 100  # Normalize dose
                        if not sheet_labeled:  # Use the same `sheet_labeled` logic here
                            ax_z.plot(pos, dose, linestyle='None', marker=marker, color=color, ms=5, alpha=0.5, label=sheet)
                            sheet_labeled = True  # Mark as labeled after the first
                        else:
                            ax_z.plot(pos, dose, linestyle='None', marker=marker, color=color, ms=5, alpha=0.5, label=None)  # No label after the first

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
                                label = f"Average FS: {fs}, Axis: {axis}, Depth: {depth}" if fs == fsl[0] and axis == scl[0] and depth == dl[0] else None
                                ax_profiles.plot(fixed_pos, avg_dose, '-k', label=label, linewidth=2)

                                # Plot clinical tolerance envelope (combined dose and positional tolerance)
                                dose_tolerance, pos_tolerance = get_clinical_criteria()  # Assuming clinical criteria function

                                # Convert positional tolerance from mm to cm
                                pos_tolerance_cm = pos_tolerance / 10.0

                                # Shift the entire profile to the left and right by positional tolerance
                                pos_shift_left = fixed_pos - pos_tolerance_cm
                                pos_shift_right = fixed_pos + pos_tolerance_cm

                                # Interpolate to create the shifted curves
                                avg_shift_left = np.interp(fixed_pos, pos_shift_left, avg_dose, left=np.nan, right=np.nan)
                                avg_shift_right = np.interp(fixed_pos, pos_shift_right, avg_dose, left=np.nan, right=np.nan)

                                # Combine the dose tolerance and positional shifts to create the final envelope
                                upper_envelope = np.fmax(avg_dose + dose_tolerance, np.fmax(avg_shift_left, avg_shift_right))
                                lower_envelope = np.fmin(avg_dose - dose_tolerance, np.fmin(avg_shift_left, avg_shift_right))

                                # Ensure that NaN values at edges do not cause asymmetry in the envelope
                                upper_envelope = np.nan_to_num(upper_envelope, nan=avg_dose[-1])
                                lower_envelope = np.nan_to_num(lower_envelope, nan=avg_dose[0])

                                # Plot the final envelope
                                ax_profiles.fill_between(fixed_pos, lower_envelope, upper_envelope, color='gray', alpha=0.2, label=f'Clinical Tolerance Envelope ({dose_tolerance}%/{pos_tolerance}mm)' if label else None)

                    elif axis == 'Z':
                        combination_key = (fs, axis, None)
                        if combination_key in avg_dose_global:
                            avg_dose = avg_dose_global[combination_key]
                            fixed_pos = fixed_pos_global[combination_key]
                            # Only label the first Z average plot once
                            label = "Average" if fs == fsl[0] and axis == 'Z' else None
                            ax_z.plot(fixed_pos, avg_dose, '-k', label=label, linewidth=2)

                            # Plot clinical tolerance envelope (combined dose and positional tolerance)
                            dose_tolerance, pos_tolerance = get_clinical_criteria()  # Assuming clinical criteria function

                            # Convert positional tolerance from mm to cm
                            pos_tolerance_cm = pos_tolerance / 10.0

                            # Shift the average curve by dose tolerance (up and down)
                            upper_bound = avg_dose + dose_tolerance
                            lower_bound = avg_dose - dose_tolerance

                            # Shift the average curve by position tolerance (left and right)
                            upper_pos_shift_right = np.interp(fixed_pos + pos_tolerance_cm, fixed_pos, avg_dose, left=np.nan, right=np.nan)
                            lower_pos_shift_left = np.interp(fixed_pos - pos_tolerance_cm, fixed_pos, avg_dose, left=np.nan, right=np.nan)

                            # Ensure that NaN values at edges do not cause asymmetry in the envelope
                            upper_pos_shift_right = np.nan_to_num(upper_pos_shift_right, nan=avg_dose[-1])
                            lower_pos_shift_left = np.nan_to_num(lower_pos_shift_left, nan=avg_dose[0])

                            # Combine all bounds to create the final envelope
                            upper_envelope = np.fmax(upper_bound, upper_pos_shift_right)
                            lower_envelope = np.fmin(lower_bound, lower_pos_shift_left)

                            # Plot the final envelope
                            ax_z.fill_between(fixed_pos, lower_envelope, upper_envelope, color='gray', alpha=0.2, label=f'Clinical Tolerance Envelope ({dose_tolerance}%/{pos_tolerance}mm)' if label else None)
    except NameError:
        print("Average data does not exist. Run make_avg first to calculate averages.")

    # Finalize the plots
    if fig_profiles:
        ax_profiles.legend(loc='upper right', fontsize=10, markerscale=1.5)
        fig_profiles.show()
    
    if fig_z:
        ax_z.legend(loc='upper right', fontsize=10, markerscale=1.5)
        fig_z.show()

    plt.tight_layout()
    print("Data plotting completed!")



    
def make_avg():
    dose_tolerance, pos_tolerance = get_clinical_criteria()
    global avg_dose_global, fixed_pos_global, fsl, scl, dl

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

    # Minimum number of datasets required at each position to include in average
    min_n_str = min_datasets_var.get()
    min_n = len(selected_sheets) if min_n_str == 'All' else max(1, int(min_n_str))

    fixed_pos_z = np.arange(0, 30.05, 0.05)
    max_field_size = max([int(fs) for fs in fsl]) if fsl else 0
    fixed_pos_profiles = np.arange(-max_field_size / 2 - 20, max_field_size / 2 + 20, 0.05)

    resampled_doses_by_combination = {}

    for sheet in selected_sheets:
        df = excel_data[sheet].sort_values(by=['FS', 'Axis', 'Pos'])
        for fs in fsl:
            for axis in scl:
                if axis in ['X', 'Y', 'XY', 'YX']:
                    for depth in dl:
                        subset = df[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth)]
                        if not subset.empty:
                            pos  = subset['Pos'].values
                            dose = subset['Dose'].values / subset['Dose'].max() * 100
                            interpolator = interp1d(pos, dose, bounds_error=False, fill_value=np.nan)
                            resampled_doses_by_combination.setdefault((fs, axis, depth), []).append(
                                interpolator(fixed_pos_profiles))
                elif axis == 'Z':
                    subset = df[(df['FS'] == fs) & (df['Axis'] == axis)]
                    if not subset.empty:
                        pos  = subset['Pos'].values
                        dose = subset['Dose'].values / subset['Dose'].max() * 100
                        interpolator = interp1d(pos, dose, bounds_error=False, fill_value=np.nan)
                        resampled_doses_by_combination.setdefault((fs, axis, None), []).append(
                            interpolator(fixed_pos_z))

    avg_dose_global  = {}
    fixed_pos_global = {}

    for combination, doses in resampled_doses_by_combination.items():
        fs, axis, depth = combination
        fixed_pos = fixed_pos_profiles if axis in ['X', 'Y', 'XY', 'YX'] else fixed_pos_z

        stack  = np.array(doses)                          # (n_sheets, n_pos)
        counts = np.sum(~np.isnan(stack), axis=0)         # how many datasets at each point
        avg    = np.nanmean(stack, axis=0)
        avg[counts < min_n] = np.nan                      # blank positions below threshold

        # Trim to the contiguous valid extent
        valid = ~np.isnan(avg)
        if valid.any():
            start = int(np.argmax(valid))
            end   = int(len(valid) - np.argmax(valid[::-1]))
            fixed_pos_global[combination] = fixed_pos[start:end]
            avg_dose_global[combination]  = avg[start:end]

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
sheet_listbox = tk.Listbox(sheet_frame, selectmode="multiple", width=40, height=10, exportselection=0)
sheet_listbox.grid(row=1, column=0, padx=5, pady=5)
sheet_listbox.bind("<<ListboxSelect>>", lambda _: (populate_fsl(), populate_scl(), populate_dl()))

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
fsl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
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
scl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
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
depth_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
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
pos_entry.grid(row=0, column=5, sticky="w")

# Adding the buttons to the GUI

# Button to plot data (without affecting the average)
plot_button = ttk.Button(root, text="Plot Data", command=plot_data)
plot_button.grid(row=8, column=0, pady=10)

# Button to calculate the average and plot it
make_avg_button = ttk.Button(root, text="Make Avg", command=make_avg)
make_avg_button.grid(row=9, column=0, pady=10)

# Button to save the average data
save_avg_button = ttk.Button(root, text="Save Avg", command=save_average)
save_avg_button.grid(row=10, column=0, pady=10)

# Start the GUI loop
root.mainloop()