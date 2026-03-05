# -*- coding: utf-8 -*-
"""
Created on [7/9/2025]

@author: nknutson Nels Knutson

Code to read IBA .opax XML files and convert to standard format.
"""

import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
import os
import xml.etree.ElementTree as ET
import datetime


def select_file_or_folder():
    path = filedialog.askopenfilename(filetypes=[("IBA OPAX Files", "*.opax")])
    if not path:
        path = filedialog.askdirectory()
    file_entry.delete(0, tk.END)
    file_entry.insert(0, path)
    print("Selected:", path)

def run_conversion():
    path = file_entry.get()
    cl = ['Depth', 'Pos', 'Dose', 'FS', 'Axis', 'SSD','Detector','Energy']


    df = pd.DataFrame(columns=cl)
    processed_files = 0

    files_to_process = []

    if os.path.isfile(path):
        files_to_process.append(path)
    elif os.path.isdir(path):
        files_to_process = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.opax')]
    else:
        status_label.config(text="Invalid file or folder selected.")
        return

    progress_bar['maximum'] = len(files_to_process)

    for fn in files_to_process:
        try:
            tree = ET.parse(fn)
            xml_root = tree.getroot()

            measurement_params = xml_root.find('MeasurementParams')
            field_params = xml_root.find('FieldParams')
            equipment_params = xml_root.find('EquipmentParams')
            # --- Energy parsing (MV + X/FFF) ---
            energy_key = ""
            energy_mv = None
            
            # Energy is typically in FieldParams/Energy
            try:
                e_txt = field_params.find('Energy').text if field_params is not None else None
                if e_txt:
                    energy_mv = float(str(e_txt).split()[0])   # accepts "6", "6.0", "6 MV"
            except Exception:
                energy_mv = None
            
            # FFF flag is in FieldParams/IsFlatteningFilterFree
            fff_flag = False
            try:
                fff_txt = field_params.find('IsFlatteningFilterFree').text if field_params is not None else None
                if fff_txt is not None:
                    fff_flag = str(fff_txt).strip().lower() in ("true", "1", "yes")
            except Exception:
                pass
            
            if energy_mv is not None:
                mv = int(round(energy_mv))
                suffix = "FFF" if fff_flag else "X"
                energy_key = f"{mv}{suffix}"

            # Detector comes from EquipmentParams
            detector_elem = equipment_params.find('FieldDetectorName') if equipment_params is not None else None
            detector_val = detector_elem.text.strip() if detector_elem is not None else "Unknown"
            
            # SSD comes from FieldParams
            ssd_elem = field_params.find('SSD') if field_params is not None else None
            ssd_val = float(ssd_elem.text) / 10.0 if ssd_elem is not None else None  # mm → cm
            
            # Keep ScanType and StartDepth from MeasurementParams
            scan_type = measurement_params.find('ScanType').text if measurement_params is not None else "N/A"
            start_depth = measurement_params.find('StartPositionDepth').text if measurement_params is not None else "N/A"

            
            # 🔧 Updated axis assignment
            scan_type = (measurement_params.find('ScanType').text or "").lower()

            st = (measurement_params.find('ScanType').text or "").strip().lower()

            if "crossline" in st:
                axis = "X"
            elif "inline" in st:
                axis = "Y"
            elif "diagonal" in st:
                if "mmpp" in st or "ppmm" in st:
                    axis = "XY"   # one diagonal
                elif "mppm" in st or "pmmp" in st:
                    axis = "YX"   # the other diagonal
                else:
                    axis = "XY"   # safe fallback
            else:
                axis = "Z"


            field_params = xml_root.find('FieldParams')
            fs_crossline = field_params.find('FieldSizeCrossline').text if field_params is not None else "N/A"

            positions_block = xml_root.find('Positions')

            if positions_block is not None:
                for pos in positions_block.findall('Position'):
                    crossline = pos.find('Crossline')
                    inline = pos.find('Inline')
                    depth_elem = pos.find('Depth')
                    dose = pos.find('Dose')

                    dose_val = dose.text.strip() if dose is not None and dose.text is not None else None

                    if axis == 'Z':  # PDD scan
                        depth_val = depth_elem.text.strip() if depth_elem is not None and depth_elem.text is not None else start_depth
                        pos_val = depth_val
                    
                    elif axis == 'X':
                        depth_val = start_depth
                        pos_val = crossline.text.strip() if crossline is not None and crossline.text is not None else None
                    
                    elif axis == 'Y':
                        depth_val = start_depth
                        pos_val = inline.text.strip() if inline is not None and inline.text is not None else None
                    
                    elif axis in ['XY', 'YX']:
                        depth_val = start_depth
                        crossline_val = crossline.text.strip() if crossline is not None and crossline.text is not None else None
                        inline_val = inline.text.strip() if inline is not None and inline.text is not None else None
                    
                        if crossline_val and inline_val:
                            try:
                                x = float(crossline_val)
                                y = float(inline_val)
                                pos_val = np.sqrt(x**2 + y**2)*np.sign(x)
                            except ValueError:
                                pos_val = None
                        else:
                            pos_val = None
                    if pos_val and dose_val and depth_val:
                        try:
                            df.loc[df.shape[0]] = [
                                float(depth_val) / 10.0,
                                float(pos_val) / 10.0,
                                float(dose_val),
                                float(fs_crossline) / 10.0,
                                axis,
                                ssd_val,
                                detector_val,
                                energy_key,                                
                            ]
                        except ValueError:
                            continue

            processed_files += 1

        except Exception as e:
            print(f"Error processing {fn}: {e}")

        # Update progress bar
        progress_bar['value'] = processed_files
        root.update_idletasks()

    # Save to Excel
    # Save to Excel with two sheets
    # wherever you set your output filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M")
    output_path = os.path.join(os.path.dirname(path), f"IBAOutput_{timestamp}.xlsx")
    #output_path = os.path.join(os.path.dirname(path), 'IBAOutput.xlsx')
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df[df['Axis'] == 'Z'].to_excel(writer, sheet_name='Depth Scans', index=False)
        df[df['Axis'] != 'Z'].to_excel(writer, sheet_name='Profile Scans', index=False)


    # Build simplified PTW-style summary with total scans
    summary = f"Processed Files: {processed_files}\n"

    if not df.empty:
        unique_axes = df['Axis'].unique()
        unique_fs = df['FS'].unique()
        unique_detectors = df['Detector'].unique()
        unique_ssd = df['SSD'].dropna().unique()
        unique_energy = df['Energy'].dropna().unique()
        summary += "Energies: " + ", ".join(sorted(map(str, unique_energy))) + "\n"

        summary += "Detectors Used: " + ", ".join(unique_detectors) + "\n"
        summary += "SSDs Detected: " + ", ".join([str(np.round(s,2)) for s in sorted(unique_ssd)]) + " cm\n"

        summary += f"Scan Types: {', '.join(unique_axes)}\n"

        # Separate Z (depth) scans from lateral scans
        lateral_depths = df[df['Axis'] != 'Z']['Depth'].unique()
        depth_scan_present = 'Z' in unique_axes

        if depth_scan_present:
            summary += "Depths Found: PDD Scan (variable depths)"
            if len(lateral_depths) > 0:
                summary += " + " + ", ".join([str(np.round(d,2)) for d in sorted(lateral_depths)]) + " cm\n"
            else:
                summary += "\n"
        else:
            summary += "Depths Found: " + ", ".join([str(np.round(d,2)) for d in sorted(lateral_depths)]) + " cm\n"

        summary += "Field Sizes Detected: " + ", ".join([str(np.round(fs,2)) for fs in sorted(unique_fs)]) + " cm\n"
        summary += f"Total Scans Processed: {processed_files}\n"

    summary += f"\nOutput saved as {os.path.basename(output_path)}\n"


    # Display in summary_text box
    summary_text.config(state=tk.NORMAL)
    summary_text.delete(1.0, tk.END)
    summary_text.insert(tk.END, summary)
    summary_text.config(state=tk.DISABLED)

    # Update final status
    status_label.config(text="✅ Processing complete and saved.")

# Initialize main window
root = tk.Tk()
root.title("IBA data conversion tool V1.0")

# File selection section
file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file_or_folder).grid(row=0, column=2, padx=5)

# Convert Data button
run_button = ttk.Button(root, text="Convert Data", command=run_conversion)
run_button.grid(row=6, column=0, pady=10)

# Progress bar
progress_bar = ttk.Progressbar(root, orient='horizontal', mode='determinate', length=400)
progress_bar.grid(row=7, column=0, padx=10, pady=5)

status_label = ttk.Label(root, text="Status: Ready")
status_label.grid(row=8, column=0, columnspan=2, pady=10)

summary_text = tk.Text(root, width=60, height=10, wrap='word')
summary_text.grid(row=9, column=0, columnspan=2, padx=10, pady=10)
def _on_close():
    # release Spyder/IPython console cleanly
    try:
        root.quit()      # stop mainloop
    finally:
        root.destroy()   # destroy the window

root.protocol("WM_DELETE_WINDOW", _on_close)
root.mainloop()

# Start the GUI loop
root.mainloop()
