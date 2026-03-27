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

        # If no common depths found, return
        if not common_depths:
            print("No common depths found!")
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

    # Ensure required data is available
    if 'avg_dose_global' not in globals() or not avg_dose_global:
        print("Average data does not exist. Run make_avg first to calculate averages.")
        return
    if 'fsl' not in globals() or not fsl:
        print("Field sizes are not defined. Run make_avg first to define them.")
        return
    if 'scl' not in globals() or not scl:
        print("Scan types are not defined. Run make_avg first to define them.")
        return
    if 'dl' not in globals():
        print("Depths are not defined. Run make_avg first to define them.")
        return

    # Create a list to store all average data
    avg_data_list = []

    # Store metadata information
    metadata = {
        'Date Created': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'Sheets Used': ', '.join([sheet_listbox.get(i) for i in sheet_listbox.curselection()]),
        'Field Sizes': ', '.join(map(str, fsl)),
        'Scan Types': ', '.join(scl),
        'Depths': ', '.join(map(str, dl)) if dl else 'All Depths'
    }

    # Loop through each field size, axis, and depth combination in the average data
    for combination_key, avg_dose in avg_dose_global.items():
        fs, axis, depth = combination_key
        fixed_pos = fixed_pos_global[combination_key]

        # Create a dictionary to store the average data in the same format as the input data
        avg_data = {
            'FS': fs,
            'Axis': axis,
            'Depth': depth if depth is not None and axis != 'Z' else '',
            'Pos': fixed_pos,
            'Dose': avg_dose
        }
        avg_data_list.append(pd.DataFrame(avg_data))

    # Concatenate all average data into a single DataFrame
    df_avg_all = pd.concat(avg_data_list, ignore_index=True)

    # Create a new DataFrame for metadata
    df_metadata = pd.DataFrame(list(metadata.items()), columns=['Description', 'Value'])

    # Prompt user for save location
    output_file = filedialog.asksaveasfilename(
        defaultextension=".xlsx",
        filetypes=[("Excel Files", "*.xlsx")],
        title="Save Average Data As")
    if not output_file:
        return  # User cancelled

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Write the average data to a single sheet named 'Average Data'
        df_avg_all.to_excel(writer, sheet_name='Average Data', index=False)

        # Write the metadata to a separate sheet named 'Metadata'
        df_metadata.to_excel(writer, sheet_name='Metadata', index=False)

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
    color_map = plt.colormaps["tab20"].resampled(len(selected_sheets))

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
    # Get clinical criteria from user input (or use defaults)
    dose_tolerance, pos_tolerance = get_clinical_criteria()
    global avg_dose_global, fixed_pos_global, fsl, scl, dl  # Store the global average and position data


    # Retrieve selected field sizes, scan types, and depths from the GUI
    selected_fsl_indices = fsl_listbox.curselection()
    selected_scl_indices = scl_listbox.curselection()
    selected_depth_indices = depth_listbox.curselection()

    # Ensure field size, scan type, and depth are selected
    fsl = [float(fsl_listbox.get(i)) for i in selected_fsl_indices]  
    scl = [scl_listbox.get(i) for i in selected_scl_indices]  
    dl = [float(depth_listbox.get(i)) for i in selected_depth_indices]  

    if not fsl or not scl:
        print("Field size or scan types not selected!")
        return

    # Set up fixed position grids for Z profiles and other profiles
    fixed_pos_z = np.arange(0, 30.05, 0.05)  # Z profiles (0 to 30 cm)
    max_field_size = max([int(fs) for fs in fsl]) if fsl else 0
    fixed_pos_profiles = np.arange(-max_field_size / 2 - 20, max_field_size / 2 + 20, 0.05)  # Profile grids

    resampled_doses_by_combination = {}  # Store resampled doses by FS, Axis, and Depth
    profile_ranges_by_combination = {}  # Store the minimum profile range for each combination

    # Get selected sheets from the GUI
    selected_sheets = [sheet_listbox.get(i) for i in sheet_listbox.curselection()]
    if not selected_sheets:
        print("No sheets selected!")
        return

    # Loop over selected sheets and resample data
    for sheet in selected_sheets:
        df = excel_data[sheet].sort_values(by=['FS', 'Axis', 'Pos'])

        # Resample doses for each FS, Axis, and Depth combination
        for fs in fsl:
            for axis in scl:
                if axis in ['X', 'Y', 'XY', 'YX']:
                    for depth in dl:
                        # Resample for profiles
                        if not df[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth)].empty:
                            pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Pos'].values
                            dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis) & (df['Depth'] == depth), 'Dose'].values
                            dose = dose / dose.max() * 100  # Normalize dose

                            interpolator = interp1d(pos, dose, bounds_error=False, fill_value="extrapolate")
                            resampled_dose = interpolator(fixed_pos_profiles)

                            # Store resampled data
                            combination_key = (fs, axis, depth)
                            if combination_key not in resampled_doses_by_combination:
                                resampled_doses_by_combination[combination_key] = []
                            resampled_doses_by_combination[combination_key].append(resampled_dose)

                            # Update the minimum profile range for this specific combination
                            profile_range = pos[-1] - pos[0]
                            if combination_key not in profile_ranges_by_combination:
                                profile_ranges_by_combination[combination_key] = profile_range
                            else:
                                profile_ranges_by_combination[combination_key] = min(profile_ranges_by_combination[combination_key], profile_range)

                elif axis == 'Z':
                    # Resample for Z profiles
                    if not df[(df['FS'] == fs) & (df['Axis'] == axis)].empty:
                        pos = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Pos'].values
                        dose = df.loc[(df['FS'] == fs) & (df['Axis'] == axis), 'Dose'].values
                        dose = dose / dose.max() * 100  # Normalize dose

                        interpolator = interp1d(pos, dose, bounds_error=False, fill_value="extrapolate")
                        resampled_dose = interpolator(fixed_pos_z)

                        # Store resampled data
                        combination_key = (fs, axis, None)  # No depth for Z profiles
                        if combination_key not in resampled_doses_by_combination:
                            resampled_doses_by_combination[combination_key] = []
                        resampled_doses_by_combination[combination_key].append(resampled_dose)

    # Calculate averages and store them globally
    avg_dose_global = {}
    fixed_pos_global = {}  # Store corresponding fixed positions

    for combination, doses in resampled_doses_by_combination.items():
        avg_dose = np.mean(doses, axis=0)  # Average resampled doses
        avg_dose_global[combination] = avg_dose

        # Store corresponding fixed positions for profiles and Z
        fs, axis, depth = combination
        if axis in ['X', 'Y', 'XY', 'YX']:
            fixed_pos_global[combination] = fixed_pos_profiles
        elif axis == 'Z':
            fixed_pos_global[combination] = fixed_pos_z

    # Post-process to remove data beyond the shortest profile range for each combination
    for combination in avg_dose_global.keys():
        fs, axis, depth = combination
        if axis in ['X', 'Y', 'XY', 'YX']:
            # Get the current positions and average dose values
            positions = fixed_pos_global[combination]
            avg_dose = avg_dose_global[combination]

            # Determine the valid range based on the shortest profile range for this combination
            min_profile_range = profile_ranges_by_combination[combination]
            valid_indices = np.where((positions >= -min_profile_range / 2) & (positions <= min_profile_range / 2))[0]

            if len(valid_indices) > 0:
                start_index = valid_indices[0]
                end_index = valid_indices[-1]

                # Keep only the valid portion of the profile
                fixed_pos_global[combination] = positions[start_index:end_index + 1]
                avg_dose_global[combination] = avg_dose[start_index:end_index + 1]

    print("Average calculation and refined post-processing completed!")



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

# Populating the listbox with sheet names from the Excel file (Done in `select_file()` function)


# Data Selection Frame
data_selection_frame = ttk.LabelFrame(root, text="Data Selection", padding="10")
data_selection_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))  # Adjust row number to fit below file and sheet sections

# Field Size Selection
fsl_label = ttk.Label(data_selection_frame, text="Select Field Sizes")
fsl_label.grid(row=0, column=0, sticky="w")
fsl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
fsl_listbox.grid(row=1, column=0, padx=5, pady=5)

# Scan Types Selection
scl_label = ttk.Label(data_selection_frame, text="Select Scan Types")
scl_label.grid(row=0, column=1, sticky="w", padx=20)
scl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
scl_listbox.grid(row=1, column=1, padx=5, pady=5)

# Depth Selection 
depth_label = ttk.Label(data_selection_frame, text="Select Depths")
depth_label.grid(row=0, column=2, sticky="w", padx=20)
depth_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
depth_listbox.grid(row=1, column=2, padx=5, pady=5)


# Add input fields for dose and position tolerance
ttk.Label(root, text="Dose Tolerance [%]:").grid(row=6, column=0, sticky="w")
dose_entry = ttk.Entry(root)
dose_entry.insert(0, "1")  # Default value of 1%
dose_entry.grid(row=6, column=0)

ttk.Label(root, text="Position Tolerance [mm]:").grid(row=7, column=0, sticky="w")
pos_entry = ttk.Entry(root)
pos_entry.insert(0, "1")  # Default value of 1 mm
pos_entry.grid(row=7, column=0)

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