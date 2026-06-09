import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
from scipy import interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from center import center as center_function
from comp import dosedif as dosedif
from comp import dta as dtafunc
from scipy.optimize import minimize_scalar
from gamma import gamma as g  # your existing gamma function

# GUI Initialization
root = tk.Tk()
root.title("Profile Compare Tool (Phantom, Site, Profile Type & Analysis with Gamma Criteria)")

# File Selection Frame
def select_file():
    file_path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)
    if file_path:
        try:
            sheets = pd.ExcelFile(file_path).sheet_names
            sheet1_combo['values'] = sheets
            sheet2_combo['values'] = sheets
            sheet1_combo.set(sheets[0])
            sheet2_combo.set(sheets[1] if len(sheets) > 1 else sheets[0])
            populate_common_combinations()

        except Exception as e:
            print(f"Error reading sheets: {e}")

file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

# Sheet Selection Frame
sheet_frame = ttk.Frame(root, padding="10")
sheet_frame.grid(row=1, column=0, sticky="ew")
ttk.Label(sheet_frame, text="Select Sheet 1:").grid(row=0, column=0, sticky="w")
sheet1_combo = ttk.Combobox(sheet_frame, width=30)
sheet1_combo.grid(row=0, column=1, padx=5)
ttk.Label(sheet_frame, text="Select Sheet 2:").grid(row=1, column=0, sticky="w")
sheet2_combo = ttk.Combobox(sheet_frame, width=30)
sheet2_combo.grid(row=1, column=1, padx=5)

# Data Selection Frame
data_selection_frame = ttk.LabelFrame(root, text="Data Selection", padding="10")
data_selection_frame.grid(row=2, column=0, sticky="ew")

# Phantom Selection
phantom_label = ttk.Label(data_selection_frame, text="Select Phantom Types")
phantom_label.grid(row=0, column=0, sticky="w")
phantom_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
phantom_listbox.grid(row=1, column=0, padx=5, pady=5)

# Site Selection
site_label = ttk.Label(data_selection_frame, text="Select Sites")
site_label.grid(row=0, column=1, sticky="w", padx=20)
site_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
site_listbox.grid(row=1, column=1, padx=5, pady=5)

# Profile Type Selection
profile_label = ttk.Label(data_selection_frame, text="Select Profile Types (Axial, Sagittal, Coronal)")
profile_label.grid(row=0, column=2, sticky="w", padx=20)
profile_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
profile_listbox.grid(row=1, column=2, padx=5, pady=5)

# Analysis Type Selection
analysis_frame = ttk.LabelFrame(root, text="Analysis Type", padding="10")
analysis_frame.grid(row=3, column=0, sticky="ew")
analysis_var = tk.StringVar(master=root,value="None")
ttk.Radiobutton(analysis_frame, text="None", variable=analysis_var, value="None").grid(row=4, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Dose Difference", variable=analysis_var, value="Dose Difference").grid(row=0, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="DTA", variable=analysis_var, value="DTA").grid(row=1, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Composite", variable=analysis_var, value="Composite").grid(row=2, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Gamma", variable=analysis_var, value="Gamma").grid(row=3, column=0, sticky="w")

# Gamma Criteria Inputs
gamma_criteria_frame = ttk.LabelFrame(root, text="Gamma Criteria", padding="10")
gamma_criteria_frame.grid(row=4, column=0, sticky="ew")
ttk.Label(gamma_criteria_frame, text="Dose Difference [%]:").grid(row=0, column=0, sticky="w")
gamma_dd_entry = ttk.Entry(gamma_criteria_frame, width=10)
gamma_dd_entry.grid(row=0, column=1, padx=5)
gamma_dd_entry.insert(0, "7")
ttk.Label(gamma_criteria_frame, text="Distance to Agreement [mm]:").grid(row=1, column=0, sticky="w")
gamma_dta_entry = ttk.Entry(gamma_criteria_frame, width=10)
gamma_dta_entry.grid(row=1, column=1, padx=5)
gamma_dta_entry.insert(0, "5")

# Processing Options Frame
processing_frame = ttk.LabelFrame(root, text="Processing Options", padding="10")
processing_frame.grid(row=5, column=0, sticky="ew")  # Move Run button to row=6

ttk.Label(processing_frame, text="Shift Sheet 1 [cm]:").grid(row=0, column=0, sticky="w")
shift1_entry = ttk.Entry(processing_frame, width=10)
shift1_entry.grid(row=0, column=1, padx=5)
shift1_entry.insert(0, "0")

ttk.Label(processing_frame, text="Shift Sheet 2 [cm]:").grid(row=0, column=2, sticky="w", padx=20)
shift2_entry = ttk.Entry(processing_frame, width=10)
shift2_entry.grid(row=0, column=3, padx=5)
shift2_entry.insert(0, "0")
auto_align_var = tk.IntVar()
auto_align_check = tk.Checkbutton(processing_frame, text="Auto-align Sheet 2 to Sheet 1", variable=auto_align_var)
auto_align_check.grid(row=2, column=0, sticky="w", padx=5, pady=2)





def populate_common_combinations():
    try:
        df1 = pd.read_excel(file_entry.get(), sheet_name=sheet1_combo.get())
        df2 = pd.read_excel(file_entry.get(), sheet_name=sheet2_combo.get())

        # Extract real combinations from each sheet
        comb1 = df1[['Phantom', 'Site', 'ProfileType']].drop_duplicates()
        comb2 = df2[['Phantom', 'Site', 'ProfileType']].drop_duplicates()

        # Only keep exact matches from both sheets
        valid_combos = pd.merge(comb1, comb2, how='inner', on=['Phantom', 'Site', 'ProfileType'])

        # Pull valid values from the filtered combinations
        common_phantoms = sorted(valid_combos['Phantom'].unique())
        common_sites = sorted(valid_combos['Site'].unique())
        common_profiles = sorted(valid_combos['ProfileType'].unique())

        # Populate listboxes
        phantom_listbox.delete(0, tk.END)
        site_listbox.delete(0, tk.END)
        profile_listbox.delete(0, tk.END)

        for p in common_phantoms:
            phantom_listbox.insert(tk.END, p)
            #phantom_listbox.selection_set(0, tk.END)

        for s in common_sites:
            site_listbox.insert(tk.END, s)
            #site_listbox.selection_set(0, tk.END)

        for pt in common_profiles:
            profile_listbox.insert(tk.END, pt)
            #profile_listbox.selection_set(0, tk.END)

    except Exception as e:
        print(f"Error populating common combinations: {e}")




def compute_optimal_shift(Pos1, Dose1, Pos2, Dose2):
    dd = 0.01          # 1% dose diff
    dta = 0.1          # 1 mm DTA
    norm = 2           # normalize to Dmax
    thres = 0.0
    interpinterval = 0.1
    shift_range = (-0.5,0.5)

    def objective(shift):
        shifted_pos2 = Pos2 + shift

        # Compute overlapping region between Pos1 and shifted Pos2
        common_min = max(min(Pos1), min(shifted_pos2))
        common_max = min(max(Pos1), max(shifted_pos2))

        if common_max <= common_min:
            return 1.0  # no overlap, max penalty

        # Mask for overlapping part of Pos1
        mask = (Pos1 >= common_min) & (Pos1 <= common_max)
        x_ref = Pos1[mask]
        y_ref = Dose1[mask]

        # Interpolate Dose2 to match x_ref
        try:
            y_test = np.interp(x_ref, shifted_pos2, Dose2, left=np.nan, right=np.nan)
        except Exception:
            return 1.0

        if np.any(np.isnan(y_test)) or len(y_test) == 0:
            return 1.0

        try:
            _, gamma_vals = g(x_ref, y_ref, x_ref, y_test, dd, dta, norm, thres, interpinterval)
        except Exception as e:
            print(f"Gamma failed at shift {shift:.2f}: {e}")
            return 1.0

        valid = ~np.isnan(gamma_vals)
        if np.sum(valid) == 0:
            return 1.0

        pass_rate = np.sum(gamma_vals[valid] <= 1.0) / np.sum(valid)
        #print(f"Shift {shift:.2f} → Gamma pass rate: {pass_rate:.3f}")
        return 1.0 - pass_rate

    result = minimize_scalar(objective, bounds=shift_range, method='bounded', options={'xatol': 0.01})
    return result.x


# Run Comparison
def run_comparison():

    # Pull basic info
    dd = float(gamma_dd_entry.get()) / 100
    dta = float(gamma_dta_entry.get()) / 10
    analysis = analysis_var.get()

    fn = file_entry.get()
    sn1 = sheet1_combo.get()
    sn2 = sheet2_combo.get()

    df1 = pd.read_excel(fn, sheet_name=sn1)
    df2 = pd.read_excel(fn, sheet_name=sn2)
    df1 = df1.sort_values(by=['Phantom', 'Site', 'ProfileType', 'Pos'])
    df2 = df2.sort_values(by=['Phantom', 'Site', 'ProfileType', 'Pos'])

    selected_phantoms = [phantom_listbox.get(i) for i in phantom_listbox.curselection()]
    selected_sites = [site_listbox.get(i) for i in site_listbox.curselection()]
    selected_profiles = [profile_listbox.get(i) for i in profile_listbox.curselection()]

    print(f"Starting comparison using analysis: {analysis}")
    print(f"DD: {dd}, DTA: {dta}")

    for phantom in selected_phantoms:
        for site in selected_sites:
            for profile in selected_profiles:
                #print(f"Processing: {phantom} | {site} | {profile}")

                sub_df1 = df1[(df1['Phantom'] == phantom) & (df1['Site'] == site) & (df1['ProfileType'] == profile)]
                sub_df2 = df2[(df2['Phantom'] == phantom) & (df2['Site'] == site) & (df2['ProfileType'] == profile)]

                if sub_df1.empty or sub_df2.empty:
                    print(f"Skipping — data missing in one of the sheets: {phantom} | {site} | {profile}")

                    continue

                x1 = sub_df1['Pos'].values
                y1 = sub_df1['Dose'].values
                x2 = sub_df2['Pos'].values
                y2 = sub_df2['Dose'].values
                # Apply user-provided shift (convert to float, default to 0 if blank or invalid)
                try:
                    shift1 = float(shift1_entry.get())
                except ValueError:
                    shift1 = 0
                if auto_align_var.get():
                    print(f"Auto-aligning: {phantom} | {site} | {profile}")
                    best_shift = compute_optimal_shift(x1, y1, x2, y2)  # using g(...) inside
                    shift2 = best_shift
                
                    # Optional: update GUI to reflect calculated shift
                    shift2_entry.delete(0, tk.END)
                    shift2_entry.insert(0, f"{best_shift:.2f}")
                else:
                    try:
                        shift2 = float(shift2_entry.get())
                    except ValueError:
                        shift2 = 0
                
                x1 = x1 + shift1
                x2 = x2 + shift2

                               
                if analysis == "Composite" or analysis == "Gamma":
                    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [1.5, 1, 1]})
                elif analysis in ["Dose Difference", "DTA"]:
                    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 6), gridspec_kw={'height_ratios': [1.5, 1]})
                else:
                    fig, ax0 = plt.subplots(1, 1, figsize=(12, 5))  # 'None'
    
    
                # Plot profiles
                ax0.plot(x1, y1, '+r', label=f'{sn1}')
                ax0.plot(x2, y2, '.k', label=f'{sn2}')
                ax0.set_title(f'{phantom} | {site} | {profile} | {analysis}')
                ax0.set_ylabel('Dose [Gy]')
                ax0.set_xlabel('Position [cm]')
                ax0.legend()
                ax0.grid(False)
                shared_xlim = ax0.get_xlim()
                
                if analysis == "None":
                    plt.tight_layout()
                    plt.show()
                    continue
    
                elif analysis == "Dose Difference":
                    difx, difv = dosedif(x1, y1, x2, y2, 3)  # norm=3
                
                    pass_mask = np.abs(difv) <= dd          # dd is fraction (e.g., 0.01)
                    ax1.plot(difx[pass_mask], difv[pass_mask] * 100, '.g', )
                    ax1.plot(difx[~pass_mask], difv[~pass_mask] * 100, '.r')
                
                    ax1.set_ylabel('Dose Difference [% of Max]')
                    ax1.set_xlabel(f"Position [cm]  (Tol = ±{dd*100:.1f}%)")
                    ax1.set_xlim(shared_xlim)
                   
                    ax1.grid(False)
                    plt.tight_layout()
                    plt.show()
                
                elif analysis == "DTA":
                    dtax, dtav = dtafunc(x1, y1, x2, y2,dta)
                
                    pass_mask = dtav <= dta                 # dta is in cm
                    ax1.plot(dtax[pass_mask], dtav[pass_mask], '.g' )
                    ax1.plot(dtax[~pass_mask], dtav[~pass_mask], '.r')
                
                    
                    ax1.set_xlabel(f"Position [cm]  (Tol = {dta:.2f} cm)")
                    ax1.set_xlim(shared_xlim)
                    ax1.grid(False)
                    plt.tight_layout()
                    plt.show()

    
                elif analysis == "Composite":
                    difx, difv = dosedif(x1, y1, x2, y2, 3)
                    dtax, dtav = dtafunc(x1, y1, x2, y2,dta)
                    
                    # Align onto the same x for fair, pointwise comparison
                    x0 = max(difx[0], dtax[0])
                    x1_ = min(difx[-1], dtax[-1])
                    mask = (difx >= x0) & (difx <= x1_)
                    difx_aligned = difx[mask]
                    difv_aligned = difv[mask]
                
                    # Interpolate DTA onto difx grid
                    dtav_on_difx = np.interp(difx_aligned, dtax, dtav)
                
                    # Fail if BOTH exceed tolerances (dd as fraction, dta in cm)
                    fail_mask = (np.abs(difv_aligned) > dd) & (dtav_on_difx > dta)
                    total_pts = len(difx_aligned)
                    num_fails = int(fail_mask.sum())
                    pass_rate = 100 * (1 - num_fails / total_pts) if total_pts else 0
                
                    # Plots (unchanged style)
                    ax1.plot(dtax, dtav * 10, '.g')  # DTA in mm
                    ax1.set_ylabel('DTA [mm]')
                    ax1.set_xlim(shared_xlim)
                    ax1.grid(False)
                
                    ax2.plot(difx, difv * 100, '.g')  # Dose diff in %
                    # mark composite failures in red at the aligned x’s
                    ax1.plot(difx_aligned[fail_mask], (dtav_on_difx[fail_mask] * 10), '.r')
                    ax2.plot(difx_aligned[fail_mask], (difv_aligned[fail_mask] * 100), '.r')
                
                    ax2.set_ylabel('Dose Difference [%]')
                    ax2.set_xlabel(f'Points outside of {dd*100:.1f}% & {dta*10:.1f}mm '
                                   f'{num_fails}/{total_pts}   Pass Rate: {pass_rate:.2f}%')
                
                    print(f"{phantom} | {site} | {profile} → Composite Failures: "
                          f"{num_fails}/{total_pts}   Pass Rate: {pass_rate:.1f}%")
                
                    ax2.set_xlim(shared_xlim)
                    plt.tight_layout(pad=0)
                    plt.show(block=False)
                    plt.pause(0.1)



    
    
                elif analysis == "Gamma":
                    gx, gv = g(x1, y1, x2, y2, dd, dta, 3, 0, 0.01)
                    
                    gv = np.nan_to_num(gv, nan=2.0)
                    gx1 = gx[gv > 1]
                    gx2 = gx[gv <= 1]
                    gv1 = gv[gv > 1]
                    gv2 = gv[gv <= 1]
                    ax1.plot(gx1, gv1, '.r')
                    ax1.plot(gx2, gv2, '.g')
                    # Format axis label dynamically
                    dd_label = f"{dd * 100:.1f}%"
                    dta_label = f"{dta * 10:.1f}mm"
                    gamma_label = rf"$\Gamma$({dd_label}/{dta_label})"
                 
                    ax1.set_ylabel(gamma_label)
                    num_total = len(gv)
                    num_failed = np.sum(gv > 1)
                    pass_rate = 100 * (1 - num_failed / num_total)
                    ax1.set_xlabel(f'Points with Γ > 1 : {num_failed}/{num_total} Pass Rate : {pass_rate:.1f}%')
                    print(f"{phantom} | {site} | {profile} → Points with Γ > 1: {num_failed}/{num_total}   Pass Rate: {pass_rate:.1f}%")
    
    
                    ax1.set_xlim(shared_xlim)
                    ax1.grid(False)
                 
                    # Histogram
                    bins = 0.1
                    ax2.hist(gv, bins=np.arange(0, np.max(gv)+bins, bins),
                             weights=np.ones_like(gv)/len(gv), color='gray', edgecolor='black')
                    ax2.set_xlabel(gamma_label)
                    ax2.set_ylabel('Normalized Incidence')
                    ax2.set_xlim(left=0)
                    ax2.grid(False)
                   
    
                    plt.tight_layout()
                    plt.show(block=False)
                    plt.pause(0.1)  # Allow rendering



# Run Button
run_button = ttk.Button(root, text="Run Comparison", command=run_comparison)
run_button.grid(row=6, column=0, pady=10)

root.mainloop()
