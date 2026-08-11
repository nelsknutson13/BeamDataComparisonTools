# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 08:46:36 2019
@author: nknutson

PTW CSV/MCC reader -> lets user pick which scans to write out.
- Detects each BEGIN_DATA..END_DATA block as one scan (scan_id).
- Shows selection dialog listing all scans found.
- Writes only selected scans to Excel (Depth Scans / Profiles), preserving columns.

UPDATE (2026-02-10):
- Multi-file selection enabled (Browse selects multiple .mcc files)
- Convert runs through all selected files and writes one Excel per file
"""
import os, datetime, hashlib
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import pandas as pd
import numpy as np

# ---------------- Multi-file selection state ----------------
SELECTED_FILES = []   # list of file paths
SELECTED_ROOT  = None # folder chosen via Browse Folder (None when using Browse Files)

INSTRUCTIONS_TEXT = (
    "PTW Data Conversion Tool\n"
    "------------------------\n"
    "This tool reads PTW MCC files (*.mcc) and converts selected scans to Unified Beam Data Format.\n"
    "\n"
    "Instructions:\n"
    "1) Click 'Browse' and select ONE or MULTIPLE .mcc files.\n"
    "   - Use Ctrl or Shift to multi-select files.\n"
    "2) Click 'Convert Data (choose scans)'.\n"
    "3) Select which scans to export in the dialog window.\n"
    "\n"
    "\nLog:\n"
)


# ---------------- UI: file picker & main window ----------------
def append_summary(msg: str):
    summary_text.config(state=tk.NORMAL)
    summary_text.insert(tk.END, msg + "\n")
    summary_text.see(tk.END)
    summary_text.config(state=tk.DISABLED)

def clear_summary():
    summary_text.config(state=tk.NORMAL)
    summary_text.delete(1.0, tk.END)
    summary_text.config(state=tk.DISABLED)


def _set_selected_files(paths, root_folder=None):
    global SELECTED_FILES, SELECTED_ROOT
    SELECTED_FILES = list(paths)
    SELECTED_ROOT  = root_folder
    if not SELECTED_FILES:
        return
    display = SELECTED_FILES[0]
    if len(SELECTED_FILES) > 1:
        display = f"{SELECTED_FILES[0]}  (+{len(SELECTED_FILES)-1} more)"
    file_entry.delete(0, tk.END)
    file_entry.insert(0, display)
    status_label.config(text=f"Status: Ready ({len(SELECTED_FILES)} file(s) selected)")


def select_file():
    """Select one OR multiple .mcc files."""
    file_paths = filedialog.askopenfilenames(
        filetypes=[("PTW MCC files", "*.mcc"), ("All files", "*.*")]
    )
    if file_paths:
        _set_selected_files(file_paths, root_folder=None)


def select_folder():
    """Select a folder and find all .mcc files in it (including subfolders)."""
    folder = filedialog.askdirectory(title="Select folder containing .mcc files")
    if not folder:
        return
    found = sorted(str(p) for p in Path(folder).rglob("*.mcc"))
    if not found:
        messagebox.showinfo("No files found", f"No .mcc files found in:\n{folder}")
        return
    _set_selected_files(found, root_folder=folder)


def run_comparison():
    global SELECTED_FILES
    if not SELECTED_FILES:
        maybe = file_entry.get()
        if maybe and os.path.isfile(maybe):
            SELECTED_FILES = [maybe]
    if not SELECTED_FILES:
        status_label.config(text="Status: No valid file(s) selected")
        return

    # Parse all files, combining into one dataset
    all_dfs        = []
    all_scans_info = []
    combined_counts = {'TPR': 0, 'PDD': 0, 'CROSSPLANE': 0, 'INPLANE': 0}
    out_dir        = SELECTED_ROOT  # set by Browse Folder; falls back below
    id_offset      = 0
    fail_list      = []

    for fn in SELECTED_FILES:
        try:
            df, scans_info, ssd_cm, scan_counts = _parse_mcc_file(fn, id_offset)
        except Exception as e:
            fail_list.append((fn, str(e)))
            messagebox.showerror("Error", f"{os.path.basename(fn)} failed:\n{e}")
            continue
        for si in scans_info:
            si['Source'] = os.path.basename(fn)
        all_dfs.append(df)
        all_scans_info.extend(scans_info)
        id_offset += len(scans_info)
        for k in combined_counts:
            combined_counts[k] += scan_counts.get(k, 0)
        if out_dir is None:
            out_dir = os.path.dirname(fn)

    if not all_dfs:
        status_label.config(text="Status: No data loaded.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)

    # One selection dialog covering all files
    selected_ids = scan_selection_dialog(root, all_scans_info)
    if selected_ids is None:
        status_label.config(text="Status: Cancelled.")
        return

    df_sel      = combined_df[combined_df['scan_id'].isin(selected_ids)].copy()
    pdd_tmr_df  = df_sel[df_sel['Axis'] == 'Z'].copy()
    profiles_df = df_sel[df_sel['Axis'].isin(['X', 'Y', 'XY', 'YX'])].copy()

    # Summary log
    summary = []
    summary.append(f"Files processed: {len(all_dfs)}" +
                   (f"  ({len(fail_list)} failed)" if fail_list else ""))
    summary.append(f"Total scans found: {len(all_scans_info)}  |  Selected: {len(selected_ids)}")
    summary.append(f"Scan types: {combined_counts}")
    axes_present = sorted({s['Axis'] for s in all_scans_info})
    summary.append(f"Axes present: {', '.join(axes_present) if axes_present else '—'}")
    det_list = sorted(set(df_sel['Detector'].astype(str)))
    summary.append("Detector(s): " + ", ".join(det_list))
    summary.append("Will write: 'Depth Scans' (Axis=Z) and 'Profiles' (Axis X/Y/XY/YX)")
    append_summary("\n".join(summary))
    append_summary("-" * 60)

    # Save one combined output
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if len(SELECTED_FILES) == 1:
        base = os.path.splitext(os.path.basename(SELECTED_FILES[0]))[0]
        out_name = f"{base}_PTWOutput_{timestamp}.xlsx"
    else:
        out_name = f"PTWOutput_combined_{timestamp}.xlsx"
    output_path = os.path.join(out_dir, out_name)

    try:
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            pdd_tmr_df.drop(columns=["scan_id"]).to_excel(writer, sheet_name='Depth Scans', index=False)
            profiles_df.drop(columns=["scan_id"]).to_excel(writer, sheet_name='Profiles', index=False)
            sel_rows = [s for s in all_scans_info if s["scan_id"] in selected_ids]
            pd.DataFrame(sel_rows).to_excel(writer, sheet_name="ScanSelection", index=False)
            pd.DataFrame(all_scans_info).to_excel(writer, sheet_name="ScanCatalog", index=False)
        status_label.config(text=f"Saved: {output_path}")
        try:
            os.startfile(output_path)
        except Exception:
            pass
    except Exception as e:
        messagebox.showerror("Error", f"Write error:\n{e}")


def _parse_mcc_file(fn, id_offset=0):
    """Parse one .mcc file. Returns (df, scans_info, ssd_cm, scan_counts).
    scan_ids start at id_offset so they stay unique across multiple files.
    """
    if not fn or not os.path.isfile(fn):
        raise FileNotFoundError(f"Invalid file: {fn}")
    try:
        with open(fn, 'r', errors='ignore') as _f:
            total_lines = sum(1 for _ in _f)
    except Exception:
        total_lines = 0
    processed_lines = 0

    cols = ['Depth', 'Pos', 'Dose', 'FS', 'Axis', 'Energy', 'Detector', 'scan_id']

    ssd_mm = None
    current_detector = None
    scan_detector = "Detector (unspecified)"

    data = False
    fs = 0.0
    axis = 'X'
    depth = 0.0
    scan_id = id_offset - 1
    scan_type = None
    energy_mv = None
    filter_tag = None
    energy_key = ""
    scan_diagonal = None
    scan_angle = None

    total_scans = 0
    scan_counts = {'TPR': 0, 'PDD': 0, 'CROSSPLANE': 0, 'INPLANE': 0}
    axes_seen = set()
    scans_info = []
    cur_points = []

    def _finalize_scan():
        nonlocal cur_points, scan_id, scan_type, axis, depth, fs
        if not cur_points:
            return
        pts = np.array(cur_points, dtype=float)
        npts = pts.shape[0]
        pos_min = float(np.nanmin(pts[:,0])) if npts else np.nan
        pos_max = float(np.nanmax(pts[:,0])) if npts else np.nan
        dose_min = float(np.nanmin(pts[:,1])) if npts else np.nan
        dose_max = float(np.nanmax(pts[:,1])) if npts else np.nan
        scans_info.append({
            "scan_id":    scan_id,
            "Type":       scan_type if scan_type else "UNKNOWN",
            "Axis":       axis,
            "Energy":     energy_key,
            "Depth_mm":   float(depth),
            "FS_mm":      float(fs),
            "npts":       int(npts),
            "pos_min_mm": pos_min,
            "pos_max_mm": pos_max,
            "dose_min":   dose_min,
            "dose_max":   dose_max,
            "Detector":   scan_detector,
        })
        cur_points.clear()

    try:
        rows = []
        with open(fn, 'r', errors='ignore') as f:
            for line in f:
                processed_lines += 1
                if total_lines and (processed_lines == 1 or processed_lines == total_lines or processed_lines % 200 == 0):
                    status_label.config(text=f"Parsing {os.path.basename(fn)}… {processed_lines/total_lines*100:.0f}%")
                    root.update_idletasks()

                s = line.strip()

                if not data:
                    if s.startswith('SSD='):
                        try: ssd_mm = float(s.split('=', 1)[1])
                        except Exception: pass
                    elif s.startswith('DETECTOR_NAME='):
                        current_detector = s.split('=', 1)[1].strip()
                    elif s.upper().startswith('ENERGY='):
                        val = s.split('=', 1)[1].strip()
                        try:
                            energy_mv = float(val.split()[0])
                        except Exception:
                            energy_mv = None
                        if energy_mv is not None:
                            mv = int(round(energy_mv))
                            suffix = "FFF" if (filter_tag or "").strip().upper() == "FFF" else "X"
                            energy_key = f"{mv}{suffix}"
                        else:
                            energy_key = ""
                        continue
                    elif s.upper().startswith('FILTER='):
                        filter_tag = s.split('=', 1)[1].split()[0].strip().upper()
                        if energy_mv is not None:
                            mv = int(round(energy_mv))
                            suffix = "FFF" if filter_tag == "FFF" else "X"
                            energy_key = f"{mv}{suffix}"
                        continue

                if 'END_DATA' in s:
                    data = False
                    _finalize_scan()
                    continue

                if 'BEGIN_DATA' in s:
                    data = True
                    scan_detector = current_detector if current_detector else "Detector (unspecified)"
                    scan_id += 1
                    total_scans += 1
                    cur_points = []
                    scan_diagonal = None
                    scan_angle = None
                    if energy_mv is not None:
                        mv = int(round(energy_mv))
                        suffix = "FFF" if (filter_tag or "").strip().upper() == "FFF" else "X"
                        energy_key = f"{mv}{suffix}"
                    continue

                if data:
                    parts = s.replace('\t\t', '\t').split('\t')
                    if len(parts) >= 2:
                        try:
                            pos_mm = float(parts[0])
                            dose   = float(parts[1])
                            cur_points.append([pos_mm, dose])
                            rows.append([depth, pos_mm, dose, fs, axis, energy_key, scan_detector, scan_id])
                        except Exception:
                            pass
                    continue

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
                    if scan_type in ('INPLANE', 'CROSSPLANE') and scan_diagonal in ('FIRST_DIAGONAL', 'SECOND_DIAGONAL'):
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
                    if scan_type in ('INPLANE', 'CROSSPLANE') and scan_diagonal not in ('FIRST_DIAGONAL', 'SECOND_DIAGONAL') and scan_angle is not None:
                        ang = scan_angle % 180.0
                        if scan_type == 'INPLANE':
                            if ang == 45.0:  axis = 'YX'
                            elif ang == 135.0: axis = 'XY'
                        else:
                            if ang == 45.0:  axis = 'XY'
                            elif ang == 135.0: axis = 'YX'
                        axes_seen.add(axis)

                if 'SCAN_DEPTH' in s and 'SCAN_DEPTH_EXTENDED' not in s:
                    try: depth = float(s.split('=', 1)[1])
                    except Exception: pass

                if s.split('=')[0] in ('FIELD_INPLANE', 'FIELD_CROSSPLANE'):
                    try: fs = float(s.split('=', 1)[1])
                    except Exception: pass

        if data:
            _finalize_scan()
        df = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)

    except Exception as e:
        raise RuntimeError(f"Failed to parse file:\n{fn}\n\n{e}")

    if df.empty:
        raise RuntimeError(f"No data found in {os.path.basename(fn)}")

    df['Pos']   = df['Pos']   / 10.0
    df['Depth'] = df['Depth'] / 10.0
    df['FS']    = df['FS']    / 10.0
    ssd_cm = (ssd_mm / 10.0) if ssd_mm is not None else None
    df['SSD'] = ssd_cm

    for si in scans_info:
        si["Depth_cm"]   = si.pop("Depth_mm")   / 10.0 if si["Depth_mm"]   is not None else np.nan
        si["FS_cm"]      = si.pop("FS_mm")       / 10.0 if si["FS_mm"]      is not None else np.nan
        si["pos_min_cm"] = si.pop("pos_min_mm")  / 10.0 if not np.isnan(si["pos_min_mm"]) else np.nan
        si["pos_max_cm"] = si.pop("pos_max_mm")  / 10.0 if not np.isnan(si["pos_max_mm"]) else np.nan

    return df, scans_info, ssd_cm, scan_counts


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

    show_source = any("Source" in s for s in scans_info)
    cols = ("Keep", "ID", "Type", "Axis", "Energy", "Detector", "Depth(cm)", "FS(cm)", "npts", "Pos range (cm)", "Dose range")
    if show_source:
        cols = cols + ("Source",)
    widths = [60, 50, 90, 60, 80, 160, 90, 80, 70, 150, 150]
    if show_source:
        widths = widths + [200]
    tree = ttk.Treeview(dlg, columns=cols, show="headings", height=16, selectmode="extended")
    for c in cols:
        tree.heading(c, text=c)
    for c, w in zip(cols, widths):
        tree.column(c, width=w, anchor="w")

    # we'll manage selection with a set of ids
    selected_ids = {s["scan_id"] for s in scans_info}  # default: select all

    def fmt_rng(a, b, prec):
        if np.isnan(a) or np.isnan(b): return "—"
        return f"{a:.{prec}f} … {b:.{prec}f}"

    # insert rows
    for s in scans_info:
        row = (
            "✔", s["scan_id"], s["Type"], s["Axis"], s.get("Energy", ""),
            s.get("Detector", ""),
            f"{s['Depth_cm']:.2f}", f"{s['FS_cm']:.2f}",
            s["npts"],
            fmt_rng(s["pos_min_cm"], s["pos_max_cm"], 3),
            fmt_rng(s["dose_min"], s["dose_max"], 4),
        )
        if show_source:
            row = row + (s.get("Source", ""),)
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
root.title("PTW data conversion tool")
root.grid_rowconfigure(4, weight=1)
root.grid_columnconfigure(0, weight=1)

file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File(s):").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=60)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse Files…", command=select_file).grid(row=0, column=2, padx=5)
ttk.Button(file_frame, text="Browse Folder…", command=select_folder).grid(row=0, column=3, padx=5)

run_button = ttk.Button(root, text="Convert Data (choose scans)", command=run_comparison)
run_button.grid(row=1, column=0, pady=8)

clear_button = ttk.Button(root, text="Clear Log", command=clear_summary)
clear_button.grid(row=2, column=0, pady=2)

status_label = ttk.Label(root, text="Status: Ready")
status_label.grid(row=3, column=0, columnspan=2, pady=4)

summary_text = tk.Text(root, width=70, height=10, wrap='word')
summary_text.grid(row=4, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")
summary_text.insert(tk.END, INSTRUCTIONS_TEXT)
summary_text.config(state=tk.DISABLED)

def _on_close():
    try:
        root.quit()
    finally:
        root.destroy()
root.protocol("WM_DELETE_WINDOW", _on_close)

root.mainloop()
