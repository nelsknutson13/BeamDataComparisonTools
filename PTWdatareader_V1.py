# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 08:46:36 2019
@author: nknutson

PTW CSV/MCC reader -> lets user pick which scans to write out.
- Detects each BEGIN_DATA..END_DATA block as one scan (scan_id).
- Shows selection dialog listing all scans found.
- Writes only selected scans to Excel (Depth Scans / Profiles), preserving columns.
"""
import os, datetime, hashlib
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import pandas as pd
import numpy as np

# ---------------- UI: file picker & main window ----------------
def select_file():
    file_path = filedialog.askopenfilename(filetypes=[("PTW MCC files", "*.mcc"), ("All files","*.*")])
    if not file_path:
        return
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)
    status_label.config(text="Status: Ready")
    summary_text.config(state=tk.NORMAL); summary_text.delete(1.0, tk.END); summary_text.config(state=tk.DISABLED)

def run_comparison():
    fn = file_entry.get()
    if not fn or not os.path.isfile(fn):
        status_label.config(text="Status: Invalid file")
        return

    # Count lines for progress
    try:
        with open(fn, 'r', errors='ignore') as _f:
            total_lines = sum(1 for _ in _f)
    except Exception:
        total_lines = 0
    processed_lines = 0

    cols = ['Depth', 'Pos', 'Dose', 'FS', 'Axis', 'Energy', 'scan_id']  # +Energy

    df = pd.DataFrame(columns=cols)

    # Header/meta extracted
    ssd_mm = None
    detector_name = None

    # live state
    data = False
    fs = 0.0
    axis = 'X'
    depth = 0.0
    scan_id = -1
    scan_type = None  # 'TPR', 'PDD', 'CROSSPLANE', 'INPLANE', etc.
    energy_mv = None        # e.g., 6.00
    filter_tag = None       # e.g., 'FFF' or 'FF'
    energy_key = ""      # e.g., '6FFF' or '6X'
    scan_diagonal = None  # 'FIRST_DIAGONAL', 'SECOND_DIAGONAL', or None
    scan_angle = None   # float degrees or None

    # tallies (for initial overview)
    total_scans = 0
    scan_counts = {'TPR': 0, 'PDD': 0, 'CROSSPLANE': 0, 'INPLANE': 0}
    axes_seen = set()

    # we’ll collect per-scan info as we go to show in the selector
    scans_info = []   # list of dict

    # temp accumulators per scan
    cur_points = []

    def _finalize_scan():
        """Called on END_DATA to record current scan’s metadata + points."""
        nonlocal cur_points, scans_info, scan_id, scan_type, axis, depth, fs
        if not cur_points:
            return
        pts = np.array(cur_points, dtype=float)  # [[pos_mm, dose], ...]
        # stats
        npts = pts.shape[0]
        pos_min = float(np.nanmin(pts[:,0])) if npts else np.nan
        pos_max = float(np.nanmax(pts[:,0])) if npts else np.nan
        dose_min = float(np.nanmin(pts[:,1])) if npts else np.nan
        dose_max = float(np.nanmax(pts[:,1])) if npts else np.nan
        scans_info.append({
            "scan_id": scan_id,
            "Type": scan_type if scan_type else "UNKNOWN",
            "Axis": axis,
            "Energy": energy_key,
            "Depth_mm": float(depth),
            "FS_mm": float(fs),
            "npts": int(npts),
            "pos_min_mm": pos_min,
            "pos_max_mm": pos_max,
            "dose_min": dose_min,
            "dose_max": dose_max,
        })
        cur_points = []

    # pass 1: parse
    try:
        rows =[]
        with open(fn, 'r', errors='ignore') as f:
            for line in f:
                processed_lines += 1
                if total_lines and (processed_lines == 1 or processed_lines == total_lines or processed_lines % 200 == 0):
                    progress_percent = (processed_lines / total_lines) * 100
                    status_label.config(text=f"Processing... {progress_percent:.1f}%")
                    root.update_idletasks()

                s = line.strip()

                # capture header/meta BEFORE data
                if not data:
                    if s.startswith('SSD='):
                        try: ssd_mm = float(s.split('=', 1)[1])
                        except Exception: pass
                    elif s.startswith('DETECTOR_NAME='):
                        detector_name = s.split('=', 1)[1].strip()
                    elif s.upper().startswith('ENERGY='):
                        val = s.split('=', 1)[1].strip()
                        try:
                            energy_mv = float(val.split()[0])          # accepts "6.00" or "6.00 MV"
                        except Exception:
                            energy_mv = None
                        # compute energy_key immediately
                        if energy_mv is not None:
                            mv = int(round(energy_mv))
                            suffix = "FFF" if (filter_tag or "").strip().upper() == "FFF" else "X"
                            energy_key = f"{mv}{suffix}"
                        else:
                            energy_key = ""
                        continue
                
                    elif s.upper().startswith('FILTER='):
                        # grab the first token (e.g., "FFF" or "FF"), not the first character
                        filter_tag = s.split('=', 1)[1].split()[0].strip().upper()
                        # recompute only if ENERGY already known
                        if energy_mv is not None:
                            mv = int(round(energy_mv))
                            suffix = "FFF" if filter_tag == "FFF" else "X"
                            energy_key = f"{mv}{suffix}"
                        continue
           
                # state switches
                if 'END_DATA' in s:
                    data = False
                    _finalize_scan()
                    continue

                if 'BEGIN_DATA' in s:
                    data = True
                    scan_id += 1
                    total_scans += 1
                    cur_points = []
                    scan_diagonal = None
                    scan_angle = None   # float degrees or None

                    if energy_mv is not None:
                        mv = int(round(energy_mv))
                        suffix = "FFF" if (filter_tag or "").strip().upper() == "FFF" else "X"
                        energy_key = f"{mv}{suffix}"
                    
                    continue

                # inside a data block: parse points
                if data:
                    parts = s.replace('\t\t', '\t').split('\t')
                    if len(parts) >= 2:
                        try:
                            pos_mm = float(parts[0])
                            dose   = float(parts[1])
                            cur_points.append([pos_mm, dose])
                            #df.loc[len(df)] = [depth, pos_mm, dose, fs, axis, energy_key, scan_id]
                            rows.append([depth, pos_mm, dose, fs, axis, energy_key, scan_id])
                        except Exception:
                            pass
                    continue

                # outside data: scan settings / meta
                if 'CURVETYPE' in s:
                    if 'TPR' in s:
                        axis = 'Z'; scan_type = 'TPR'; scan_counts['TPR'] += 1
                    elif 'PDD' in s:
                        axis = 'Z'; scan_type = 'PDD'; scan_counts['PDD'] += 1
                    elif 'CROSSPLANE' in s:
                        axis = 'X'; scan_type = 'CROSSPLANE'; scan_counts['CROSSPLANE'] += 1
                    elif 'INPLANE' in s:
                        axis = 'Y'; scan_type = 'INPLANE'; scan_counts['INPLANE'] += 1
                    axes_seen.add(axis)

                
 
                if s.upper().startswith('SCAN_DIAGONAL'):
                    raw = s.split('=', 1)[1].strip().upper()
                    if raw.startswith('FIRST'):
                        scan_diagonal = 'FIRST_DIAGONAL'
                    elif raw.startswith('SECOND'):
                        scan_diagonal = 'SECOND_DIAGONAL'
                    elif raw in ('', 'NONE', 'NO', 'FALSE'):
                        scan_diagonal = None
                    else:
                        scan_diagonal = raw
                
                    # Fix: set axis immediately if scan_type is known
                    if scan_type in ('INPLANE','CROSSPLANE') and scan_diagonal in ('FIRST_DIAGONAL','SECOND_DIAGONAL'):
                        if scan_type == 'INPLANE':
                            axis = 'YX' if scan_diagonal == 'FIRST_DIAGONAL' else 'XY'
                        else:
                            axis = 'XY' if scan_diagonal == 'FIRST_DIAGONAL' else 'YX'
                        axes_seen.add(axis)

                if s.upper().startswith('SCAN_ANGLE'):
                    raw = s.split('=', 1)[1].strip()
                    raw = raw.replace('°', '').replace(',', '.').split()[0]
                    try:
                        scan_angle = float(raw)
                    except Exception:
                        scan_angle = None
                
                    # Fix: assign XY/YX if no diagonal but we know the type
                    if scan_type in ('INPLANE','CROSSPLANE') and scan_diagonal not in ('FIRST_DIAGONAL','SECOND_DIAGONAL') and scan_angle is not None:
                        ang = scan_angle % 180.0
                        if scan_type == 'INPLANE':
                            if ang == 45.0: axis = 'YX'
                            elif ang == 135.0: axis = 'XY'
                        else:  # CROSSPLANE
                            if ang == 45.0: axis = 'XY'
                            elif ang == 135.0: axis = 'YX'
                        axes_seen.add(axis)

                if 'SCAN_DEPTH' in s and 'SCAN_DEPTH_EXTENDED' not in s:
                    try: depth = float(s.split('=', 1)[1])
                    except Exception: pass

                if s.split('=')[0] in ('FIELD_INPLANE', 'FIELD_CROSSPLANE'):
                    try: fs = float(s.split('=', 1)[1])
                    except Exception: pass

        # in case file didn’t end with END_DATA
        if data:
            _finalize_scan()
        # BUILD DF ONCE (ADD THIS BLOCK)
        if rows:
            df = pd.DataFrame(rows, columns=cols)
        else:
            df = pd.DataFrame(columns=cols)

    except Exception as e:
        messagebox.showerror("Read error", f"Failed to parse file:\n{e}")
        return

    if df.empty:
        status_label.config(text="Status: No data found")
        return

    # ---- Unit conversion (mm -> cm) ----
    df['Pos']   = df['Pos'] / 10.0
    df['Depth'] = df['Depth'] / 10.0
    df['FS']    = df['FS'] / 10.0

    # SSD & Detector
    ssd_cm = (ssd_mm / 10.0) if ssd_mm is not None else None
    df['SSD'] = ssd_cm
    df['Detector'] = detector_name if detector_name else "Detector (unspecified)"

    # ---- Show selection dialog so user chooses which scans to keep ----
    # Convert scans_info mm → cm for display
    for si in scans_info:
        si["Depth_cm"] = si.pop("Depth_mm") / 10.0 if si["Depth_mm"] is not None else np.nan
        si["FS_cm"]    = si.pop("FS_mm") / 10.0 if si["FS_mm"] is not None else np.nan
        si["pos_min_cm"] = si.pop("pos_min_mm") / 10.0 if not np.isnan(si["pos_min_mm"]) else np.nan
        si["pos_max_cm"] = si.pop("pos_max_mm") / 10.0 if not np.isnan(si["pos_max_mm"]) else np.nan

    selected_ids = scan_selection_dialog(root, scans_info)
    if selected_ids is None:
        # user cancelled → keep all
        selected_ids = {s["scan_id"] for s in scans_info}

    # filter to selected
    keep_mask = df['scan_id'].isin(selected_ids)
    df_sel = df.loc[keep_mask].copy()

    # Split outputs
    pdd_tmr_df  = df_sel[df_sel['Axis'] == 'Z'].copy()
    profiles_df = df_sel[df_sel['Axis'].isin(['X', 'Y', 'XY', 'YX'])].copy()


    # ---- UI summary ----
    status_label.config(text="Processing complete (ready to save).")

    summary = []
    summary.append(f"Total Scans Found: {len(scans_info)}")
    summary.append(f"Selected Scans: {len(selected_ids)}")
    summary.append(f"Scan Types: {scan_counts}")
    axes_present = sorted({s['Axis'] for s in scans_info})
    summary.append(f"Axes Present: {', '.join(axes_present) if axes_present else '—'}")
    if ssd_cm is not None: summary.append(f"SSD: {ssd_cm:.1f} cm")
    summary.append(f"Detector: {df['Detector'].iloc[0]}")
    summary.append("Will write: 'Depth Scans' (Axis=Z) and 'Profiles' (Axis X/Y/XY/YX)")
    summary_text.config(state=tk.NORMAL)
    summary_text.delete(1.0, tk.END)
    summary_text.insert(tk.END, "\n".join(summary))
    summary_text.config(state=tk.DISABLED)

    # ---- Save output ----
    try:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_filename = f"PTWOutput_{timestamp}.xlsx"
        output_path = os.path.join(os.path.dirname(fn), output_filename)
        with pd.ExcelWriter(output_path) as writer:
            # write chosen scans only
            pdd_tmr_df.drop(columns=["scan_id"]).to_excel(writer, sheet_name='Depth Scans', index=False)
            profiles_df.drop(columns=["scan_id"]).to_excel(writer, sheet_name='Profiles', index=False)

            # write a "ScanSelection" tab for provenance
            sel_rows = [s for s in scans_info if s["scan_id"] in selected_ids]
            all_rows = [s for s in scans_info]
            pd.DataFrame(sel_rows).to_excel(writer, sheet_name="ScanSelection", index=False)
            pd.DataFrame(all_rows).to_excel(writer, sheet_name="ScanCatalog", index=False)

        status_label.config(text=f"Saved: {output_path}")
        print(f"Data saved successfully to {output_path}")
    except Exception as e:
        messagebox.showerror("Write error", f"An error occurred:\n{e}")


# ---------------- selection dialog ----------------
def scan_selection_dialog(parent, scans_info):
    """
    Modal dialog to choose which scans to keep.
    Returns a set of selected scan_ids, or None if cancelled.
    """
    dlg = tk.Toplevel(parent)
    dlg.title("Select scans to export")
    dlg.transient(parent)
    dlg.grab_set()

    # columns to show
    cols = ("Keep", "ID", "Type", "Axis","Energy", "Depth(cm)", "FS(cm)", "npts", "Pos range (cm)", "Dose range")
    tree = ttk.Treeview(dlg, columns=cols, show="headings", height=16, selectmode="extended")
    for c in cols:
        tree.heading(c, text=c)
    widths = [60, 50, 90, 60,80, 90, 80, 70, 150, 150]
    for c, w in zip(cols, widths):
        tree.column(c, width=w, anchor="w")

    # we’ll manage selection with a set of ids
    selected_ids = {s["scan_id"] for s in scans_info}  # default: select all

    def fmt_rng(a, b, prec):
        if np.isnan(a) or np.isnan(b): return "—"
        return f"{a:.{prec}f} … {b:.{prec}f}"

    # insert rows
    for s in scans_info:
        row = (
            "✔", s["scan_id"], s["Type"], s["Axis"],s.get("Energy",""),
            f"{s['Depth_cm']:.2f}", f"{s['FS_cm']:.2f}",
            s["npts"],
            fmt_rng(s["pos_min_cm"], s["pos_max_cm"], 3),
            fmt_rng(s["dose_min"], s["dose_max"], 4),
        )
        tree.insert("", "end", iid=str(s["scan_id"]), values=row)

    # toggle selection on double-click
    def on_toggle(event=None):
        item = tree.identify_row(event.y) if event else tree.focus()
        if not item:
            return
        sid = int(item)
        if sid in selected_ids:
            selected_ids.remove(sid)
            vals = list(tree.item(item, "values")); vals[0] = " "
            tree.item(item, values=vals)
        else:
            selected_ids.add(sid)
            vals = list(tree.item(item, "values")); vals[0] = "✔"
            tree.item(item, values=vals)

    tree.bind("<Double-1>", on_toggle)

    # buttons
    btn_frame = ttk.Frame(dlg)
    ttk.Button(btn_frame, text="Select All", command=lambda: _select_all(tree, selected_ids, True)).pack(side="left", padx=4)
    ttk.Button(btn_frame, text="Select None", command=lambda: _select_all(tree, selected_ids, False)).pack(side="left", padx=4)
    ttk.Button(btn_frame, text="OK", command=dlg.destroy).pack(side="right", padx=4)
    def on_cancel():
        nonlocal selected_ids
        selected_ids = None
        dlg.destroy()
    ttk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side="right", padx=4)

    tree.grid(row=0, column=0, sticky="nsew", padx=10, pady=(10,0))
    btn_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=10)
    dlg.grid_rowconfigure(0, weight=1)
    dlg.grid_columnconfigure(0, weight=1)
    dlg.geometry("+%d+%d" % (parent.winfo_rootx()+60, parent.winfo_rooty()+60))

    dlg.wait_window()
    return selected_ids

def _select_all(tree, selected_ids, make_selected: bool):
    # helper: toggle all rows
    selected_ids.clear()
    for iid in tree.get_children(""):
        if make_selected:
            selected_ids.add(int(iid))
            vals = list(tree.item(iid, "values")); vals[0] = "✔"
        else:
            vals = list(tree.item(iid, "values")); vals[0] = " "
        tree.item(iid, values=vals)

# ---------------- main window ----------------
root = tk.Tk()
root.title("PTW data conversion tool (select scans)")

file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=60)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

run_button = ttk.Button(root, text="Convert Data (choose scans)", command=run_comparison)
run_button.grid(row=1, column=0, pady=8)

status_label = ttk.Label(root, text="Status: Ready")
status_label.grid(row=2, column=0, columnspan=2, pady=4)

summary_text = tk.Text(root, width=70, height=10, wrap='word')
summary_text.grid(row=3, column=0, columnspan=2, padx=10, pady=10)
summary_text.config(state=tk.DISABLED)

def _on_close():
    try:
        root.quit()
    finally:
        root.destroy()
root.protocol("WM_DELETE_WINDOW", _on_close)

root.mainloop()
