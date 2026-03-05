import re
import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, subprocess

last_df = None
last_meta = None

# Cache for parsed groups and labels
last_groups = None                # dict: key(tuple)-> {voltage: [arrays]}
last_group_keys = []              # list of keys in stable order
last_label_map = {}               # label -> key

# Facets the UI can group by (only ones we parse today)
available_facets = [
    "Energy", "Detector", "SSD_cm", "FieldSize_cm",
    "NominalDoseRate_MUmin", "Detector_SN", "Linac", "ScanType", "ProfileDepth_cm"
]

group_facets = ["Energy", "Detector", "SSD_cm", "FieldSize_cm"]  # default selection (leave MU/min out)




# Normalize CF at this depth (cm) per energy label
NORM_DEPTH_CM = {
    "6MV": 1.5,
    "6MV FFF": 1.4,
    "8MV FFF": 1.9,
    "10MV": 2.4,
    "10MV FFF": 2.2,
    "15MV": 3.0,
}
# Use the same mapping as dmax defaults; adjust later if you prefer
DMAX_DEPTH_CM = dict(NORM_DEPTH_CM)

def _dmax_for_energy(energy_label):
    """Return dmax (cm) for an energy label like '6MV', '6MV FFF', etc."""
    return DMAX_DEPTH_CM.get(energy_label)


# ------------------------------
# Legend: click to toggle lines
# ------------------------------
def _attach_click_legend(fig, lines, labels):
    leg = fig.axes[0].legend(lines, labels)
    lined = {}
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(True)
        lined[legline] = origline

    def on_pick(event):
        legline = event.artist
        orig = lined.get(legline)
        if orig is None:
            return
        vis = not orig.get_visible()
        orig.set_visible(vis)
        legline.set_alpha(1.0 if vis else 0.3)
        fig.canvas.draw()

    fig.canvas.mpl_connect("pick_event", on_pick)

# ------------------------------
# Parsing helpers
# ------------------------------
def parse_scan(block):
    """Extract voltage, data array, and metadata from one BEGIN_SCAN block.
    Returns: (voltage, arr, det_name, ssd_cm, energy_mv, filt, fs_cm, nominal_mu_min, det_sn, linac)
      - arr columns: depth_cm, charge_C (may be NaN), doseRate_GyPerMin
    """
    # --- Basic headers ---
    m = re.search(r'DETECTOR_HV=([0-9.]+)', block); voltage = float(m.group(1)) if m else None
    m = re.search(r'DETECTOR_NAME=([^\n\r]+)', block); det_name = m.group(1).strip() if m else ""
    m = re.search(r'SSD=([0-9.]+)', block); ssd_cm = (float(m.group(1))/10.0) if m else None
    m = re.search(r'ENERGY=([0-9.]+)', block); energy_mv = float(m.group(1)) if m else None
    m = re.search(r'FILTER=([^\n\r]+)', block); filt = m.group(1).strip().upper() if m else None
    #Looking for profile scans and if it has a depth
    m = re.search(r'SCAN_CURVETYPE=([^\n\r]+)', block)
    scan_type = m.group(1).strip().upper() if m else "PDD"
    #print(f"SCAN_CURVETYPE parsed as: '{scan_type}'")

    m = re.search(r'SCAN_DEPTH=([0-9.]+)', block)
    profile_depth_cm = float(m.group(1)) / 10.0 if m else None
    # Electrometer / timing (optional for charge calc)
    m = re.search(r'DETECTOR_CALIBRATION=([0-9.Ee+-]+)', block); det_cal = float(m.group(1)) if m else None
    m = re.search(r'MEAS_TIME=([0-9.]+)', block); meas_time = float(m.group(1)) if m else None

    # NEW: Nominal dose rate (Gy/min) -> MU/min (×100)
    m = re.search(r'EXPECTED_MAX_DOSE_RATE=([0-9.]+)', block)
    nominal_mu_min = float(m.group(1)) * 100 if m else None

    # Field size (mm) -> cm (equivalent square)
    fs_cm = None
    mfx = re.search(r'FIELD_INPLANE\s*=\s*([0-9.]+)', block)
    mfy = re.search(r'FIELD_CROSSPLANE\s*=\s*([0-9.]+)', block)
    if not (mfx and mfy):
        mfx = mfx or re.search(r'REF_FIELD_INPLANE\s*=\s*([0-9.]+)', block)
        mfy = mfy or re.search(r'REF_FIELD_CROSSPLANE\s*=\s*([0-9.]+)', block)
    if mfx and mfy:
        fx_cm = float(mfx.group(1)) / 10.0
        fy_cm = float(mfy.group(1)) / 10.0
        fs_cm = (2*fx_cm*fy_cm/(fx_cm + fy_cm)) if (fx_cm + fy_cm) else None

    # NEW: Detector serial + Linac
    m = re.search(r'DETECTOR_SN=([^\n\r]+)', block)
    det_sn = m.group(1).strip() if m else None
    m = re.search(r'LINAC=([^\n\r]+)', block)
    linac = m.group(1).strip() if m else None

    # Data block
    rows = []
    dm = re.search(r'BEGIN_DATA(.*?)END_DATA', block, flags=re.S)
    if dm:
        for ln in dm.group(1).strip().splitlines():
            parts = [p for p in re.split(r'\s+', ln.strip()) if p]
            if len(parts) >= 2:
                try:
                    depth_mm = float(parts[0])
                    dose_rate_gpm = float(parts[1])         # Gy/min
                    dose_rate_gps = dose_rate_gpm / 60.0    # Gy/s
                    charge_c = (dose_rate_gps * meas_time / det_cal) if (det_cal and meas_time) else np.nan
                    rows.append((depth_mm/10.0, charge_c, dose_rate_gpm))
                except ValueError:
                    pass
    arr = np.array(rows) if rows else np.empty((0, 3), dtype=float)

    return (voltage, arr, det_name, ssd_cm, energy_mv, filt,fs_cm, 
            nominal_mu_min, det_sn, linac, scan_type, profile_depth_cm)




def _energy_label(energy_mv, filt):
    if energy_mv is None:
        return None
    is_fff = (filt or "").strip().upper() == "FFF"
    return f"{energy_mv:g}MV" + (" FFF" if is_fff else "")

def parse_mcc(filepath):
    with open(filepath, "r", errors="ignore") as f:
        text = f.read()
    blocks = [b.split('END_SCAN')[0] for b in re.split(r'BEGIN_SCAN\s+\d+\s*', text)[1:]]
    scans = [parse_scan(b) for b in blocks if "BEGIN_DATA" in b]
    return scans
def parse_mcc_path(path, recurse=False):
    """Load scans from a single .mcc file or from all .mcc files in a folder."""
    if os.path.isdir(path):
        scans = []
        if recurse:
            walker = os.walk(path)
        else:
            walker = [(path, [], os.listdir(path))]  # only this folder
        file_count = 0
        for root, _, files in walker:
            for fn in sorted(files):
                if fn.lower().endswith(".mcc"):
                    file_count += 1
                    full = os.path.join(root, fn)
                    try:
                        scans.extend(parse_mcc(full))
                    except Exception as e:
                        print(f"[skip] {full}: {e}")
        if file_count == 0:
            print("No .mcc files found in the selected folder.")
        else:
            print(f"Loaded {file_count} .mcc file(s); {len(scans)} scan(s).")
        return scans
    else:
        scans = parse_mcc(path)
        print(f"Loaded 1 .mcc file; {len(scans)} scan(s).")
        return scans

# ------------------------------
# Grouping helpers
# ------------------------------
def _get_facet(key, facet):
    try:
        return key[group_facets.index(facet)]
    except ValueError:
        return None

def _group_label(key):
    parts = []
    for facet, val in zip(group_facets, key):
        if facet == "SSD_cm":
            txt = "?" if val is None or not np.isfinite(val) else f"{val:g} cm"
            parts.append(f"SSD={txt}")
        elif facet == "FieldSize_cm":
            txt = "?" if val is None or not np.isfinite(val) else f"{val:g} cm"
            parts.append(f"FS={txt}")
        elif facet == "NominalDoseRate_MUmin":
            txt = "?" if val is None else f"{int(val)} MU/min"
            parts.append(f"Nominal={txt}")
        elif facet == "Detector_SN":
            parts.append(f"SN={val if val not in (None,'') else '?'}")
        elif facet == "Linac":
            parts.append(f"LINAC={val if val not in (None,'') else '?'}")
        elif facet == "ProfileDepth_cm":
            txt = "?" if val is None or not np.isfinite(val) else f"{val:g} cm"
            parts.append(f"Depth={txt}")

        else:
            parts.append(str(val) if val not in (None, "") else "?")
    return " — ".join(parts)

def _group_scans(scans, facets):
    groups = {}
    for item in scans:
        v, arr, det, ssd, energy_mv, filt = item[:6]
        fs_cm           = item[6] if len(item) >= 7 else None
        nominal_mu_min  = item[7] if len(item) >= 8 else None
        det_sn          = item[8] if len(item) >= 9 else None
        linac           = item[9] if len(item) >= 10 else None
        scan_type       = item[10] if len(item) >= 11 else "PDD"
        prof_depth_cm   = item[11] if len(item) >= 12 else None     # ← ADD THIS

        meta = {
            "Energy": _energy_label(energy_mv, filt),
            "Detector": det or "",
            "SSD_cm": ssd,
            "FieldSize_cm": fs_cm,
            "NominalDoseRate_MUmin": nominal_mu_min,
            "Detector_SN": det_sn,
            "Linac": linac,
            "ScanType": scan_type,            # ← ADD THIS
            "ProfileDepth_cm": prof_depth_cm, # ← AND THIS
        }


        # stabilize float-ish keys
        if isinstance(meta["SSD_cm"], float) and np.isfinite(meta["SSD_cm"]):
            meta["SSD_cm"] = round(meta["SSD_cm"], 3)
        if isinstance(meta["FieldSize_cm"], float) and np.isfinite(meta["FieldSize_cm"]):
            meta["FieldSize_cm"] = round(meta["FieldSize_cm"], 3)
        if isinstance(meta["NominalDoseRate_MUmin"], (int, float)) and np.isfinite(float(meta["NominalDoseRate_MUmin"])):
            meta["NominalDoseRate_MUmin"] = int(round(meta["NominalDoseRate_MUmin"]))

        key = tuple(meta[f] for f in facets)
        key = key + (meta["ScanType"],)  # always include ScanType at end
        groups.setdefault(key, {}).setdefault(float(v) if v is not None else np.nan, []).append(arr)
    return groups


# ------------------------------
# Pion math + alignment
# ------------------------------
def calc_pion(depths_cm, m_low, m_high, v_low, v_high):
    mratio = m_high / m_low
    vratio = v_high / v_low
    return (1 - vratio) / (mratio - vratio)

def _align_pair(arr_low, arr_high, tol=1e-4):
    dL, yL = arr_low[:,0],  arr_low[:,1]
    dH, yH = arr_high[:,0], arr_high[:,1]
    iL = np.argsort(dL); dL, yL = dL[iL], yL[iL]
    iH = np.argsort(dH); dH, yH = dH[iH], yH[iH]

    if len(dL) == len(dH) and np.allclose(dL, dH, atol=tol):
        return dL, yL, yH

    stepL = np.median(np.diff(dL)) if len(dL) > 1 else np.inf
    stepH = np.median(np.diff(dH)) if len(dH) > 1 else np.inf
    base_d = dL if stepL <= stepH else dH

    lo = max(dL.min(), dH.min())
    hi = min(dL.max(), dH.max())
    mask = (base_d >= lo) & (base_d <= hi)
    base_d = base_d[mask]

    yL2 = np.interp(base_d, dL, yL)
    yH2 = np.interp(base_d, dH, yH)
    return base_d, yL2, yH2
# --- HV uncertainty (conservative) ---
# Spec: ±0.1% per voltage. Treat spec as a rectangular limit → u_rel = 0.1% / sqrt(3).
# For the ratio β = V_high / V_low with independent errors → u(β)/β = sqrt(2) * u_rel.
HV_REL_SPEC = 0.001                         # 0.1%
HV_REL_STD  = HV_REL_SPEC / np.sqrt(3.0)    # one-sigma (rectangular → standard)
# One-sigma relative uncertainty of a single charge reading (editable)
EM_REL_STD = 0.0005   # 0.05% → change to your instrument/repeatability
def _u_beta(beta):
    """One-sigma uncertainty of β due to HV accuracy."""
    return beta * np.sqrt(2.0) * HV_REL_STD

def _u_pion_from_beta(R, beta):
    """Propagate HV uncertainty into pion at the evaluation point (R, beta)."""
    # dπ/dβ = (1 - R) / (R - β)^2
    return np.abs((1.0 - R) / ((R - beta)**2)) * _u_beta(beta)

# ------------------------------
# GUI actions
# ------------------------------
def open_grouping_dialog():
    win = tk.Toplevel(root)
    win.title("Group by…")
    win.transient(root); win.grab_set()

    ttk.Label(win, text="Select facets (top-to-bottom order is fixed):").grid(row=0, column=0, sticky="w", padx=10, pady=(10,6))

    vars_map = {}
    for i, facet in enumerate(available_facets, start=1):
        v = tk.BooleanVar(value=(facet in group_facets))
        ttk.Checkbutton(win, text=facet, variable=v).grid(row=i, column=0, sticky="w", padx=14)
        vars_map[facet] = v

    def apply_and_close():
        chosen = [f for f in available_facets if vars_map[f].get()]
        if not chosen:
            chosen = ["Energy", "Detector", "SSD_cm", "FieldSize_cm"]
        set_group_facets(chosen)
        win.destroy()

    btns = ttk.Frame(win); btns.grid(row=len(available_facets)+1, column=0, pady=10, padx=10, sticky="e")
    ttk.Button(btns, text="OK", command=apply_and_close).pack(side="left", padx=6)
    ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="left")

def set_group_facets(chosen):
    """Update grouping facets. 
    For profile scans, automatically include ProfileDepth_cm as mandatory.
    """
    global group_facets, last_groups, last_group_keys, last_label_map

    path = file_entry.get().strip()
    if not path:
        return

    # Parse scans from the selected file or folder
    scans = parse_mcc_path(path)

    # Check if any scan is a profile (e.g., INPLANE_PROFILE, CROSSPLANE_PROFILE, etc.)
    has_profile = any(
        len(s) >= 11 and isinstance(s[10], str) and "PROFILE" in s[10].upper()
        for s in scans
    )

    # If it's a profile dataset, make ProfileDepth_cm mandatory
    if has_profile and "ProfileDepth_cm" not in chosen:
        chosen = chosen + ["ProfileDepth_cm"]

    group_facets = chosen

    # Regroup data based on chosen facets
    last_groups = _group_scans(scans, group_facets)
    last_group_keys = list(last_groups.keys())
    last_label_map = {_group_label(k): k for k in last_group_keys}

    # Update listbox
    group_listbox.delete(0, tk.END)
    for k in last_group_keys:
        group_listbox.insert(tk.END, _group_label(k))
    group_listbox.select_set(0, tk.END)


def browse_file():
    global last_groups, last_group_keys, last_label_map
    file_path = filedialog.askopenfilename(
        filetypes=[("PTW MCC files", "*.mcc"), ("All files", "*.*")]
    )
    if not file_path:
        return
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)
    set_group_facets(group_facets)  # ensures ProfileDepth_cm is auto-added for profiles

def browse_folder():
    global last_groups, last_group_keys, last_label_map
    folder_path = filedialog.askdirectory()
    if not folder_path:
        return
    file_entry.delete(0, tk.END)
    file_entry.insert(0, folder_path)
    set_group_facets(group_facets)  # ensures ProfileDepth_cm is auto-added for profiles


def _energy_from_key(key):
    try:
        idx = group_facets.index("Energy")
        return key[idx]
    except ValueError:
        return None

# ------------------------------
# Plotting + export
# ------------------------------
def plot_pion():
    global last_df, last_meta, last_groups, last_group_keys, last_label_map

    filepath = file_entry.get()
    if not filepath:
        print("No file selected")
        return

    # Ensure groups are parsed (handles manual path entry)
    if last_groups is None:
        set_group_facets(group_facets)  # builds groups + listbox and auto-adds ProfileDepth_cm for profiles


    groups = last_groups

    # Read selection from Listbox
    sel_idxs = group_listbox.curselection()
    if not sel_idxs:
        chosen_keys = set(last_group_keys)  # default to ALL
    else:
        labels = [group_listbox.get(i) for i in sel_idxs]
        chosen_keys = { last_label_map[lbl] for lbl in labels if lbl in last_label_map }

    if not chosen_keys:
        print("No groups selected.")
        return

    print(f"Plotting {len(chosen_keys)} selected groups out of {len(groups)}.")

    rows = []
    curves_summary = []

    # Figure 1: Pion vs Depth
    fig1 = plt.figure(figsize=(7.6, 5.2))
    ax1 = fig1.add_subplot(111)
    main_lines_fig1, labels_fig1 = [], []

    for key in last_group_keys:
        if key not in chosen_keys:
            continue
        energy_label = _get_facet(key, "Energy")
        det_name     = _get_facet(key, "Detector")
        ssd_cm       = _get_facet(key, "SSD_cm")
        fs_cm        = _get_facet(key, "FieldSize_cm")
        by_v = groups[key]

        uniq_vs = sorted(by_v.keys())
        hv_pairs = [(v, u) for v in uniq_vs for u in uniq_vs if abs(u - 2*v) <= 1e-6]
        if not hv_pairs:
            print(f"Skipping {key} — no HV pairs found: {list(by_v.keys())}")

            continue

        for (v_low, v_high) in hv_pairs:
            lows, highs = by_v.get(v_low, []), by_v.get(v_high, [])
            n_pairs = min(len(lows), len(highs))
            if n_pairs == 0:
                continue

            # Align each pair to a common grid
            pair_ds, pair_mL, pair_mH = [], [], []
            for i in range(n_pairs):
                d_common, mL, mH = _align_pair(lows[i][:,[0,1]], highs[i][:,[0,1]])
                idx = np.argsort(d_common)
                d_common, mL, mH = d_common[idx], mL[idx], mH[idx]
                if d_common.size < 2:
                    continue
                pair_ds.append(d_common); pair_mL.append(mL); pair_mH.append(mH)
            if not pair_ds:
                continue

            lo = max(d[0] for d in pair_ds)
            hi = min(d[-1] for d in pair_ds)
            # choose the DENSER grid as base
            base_idx = int(np.argmin([np.median(np.diff(d)) if d.size>1 else np.inf for d in pair_ds]))
            base_grid = pair_ds[base_idx]
            mask = (base_grid >= lo) & (base_grid <= hi)
            base_grid = base_grid[mask]
            if base_grid.size < 2:
                continue
            
            # --- NEW: clip to >= dmax if requested ---
            if only_beyond_dmax_var.get():
                en_label = energy_label  # already pulled from facets
                dmax_cm = _dmax_for_energy(en_label)
                if dmax_cm is not None and np.isfinite(dmax_cm):
                    dm_mask = base_grid >= float(dmax_cm)
                    base_grid = base_grid[dm_mask]
                    if base_grid.size < 2:
                        # nothing to plot for this pair beyond dmax
                        continue
                    # We'll apply dm_mask to yL/yH after interpolation below
                else:
                    print(f"[dmax] No dmax mapping for energy '{en_label}'. Using full depth range.")
            
            aligned_lows, aligned_highs, pion_list = [], [], []
            for i in range(len(pair_ds)):
                yL_full = np.interp(base_grid, pair_ds[i], pair_mL[i])
                yH_full = np.interp(base_grid, pair_ds[i], pair_mH[i])
            
                # If we created dm_mask above, it's already reflected in base_grid;
                # yL_full/yH_full are on the same (possibly clipped) grid.
                aligned_lows.append(np.column_stack([base_grid, yL_full]))
                aligned_highs.append(np.column_stack([base_grid, yH_full]))
                pion_list.append(calc_pion(base_grid, yL_full, yH_full, v_low, v_high))


            depths   = base_grid
            pion_arr = np.vstack(pion_list)

            # Label from facets + HV tag
            base_lbl = _group_label(key)
            hv_tag   = f"{int(v_high)}/{int(v_low)} V"
            label    = f"{base_lbl}, {hv_tag}"

            if pion_arr.shape[0] == 1:
                y = pion_arr[0]
                y = np.maximum(y, 1.0)  # clamp
                # HV-only uncertainty (one-sigma), then plot 95% (~±2σ)
                beta = v_high / v_low
                mL = aligned_lows[0][:,1]
                mH = aligned_highs[0][:,1]
                with np.errstate(divide='ignore', invalid='ignore'):
                    R0 = np.divide(mH, mL, out=np.full_like(mH, np.nan), where=(mL!=0))
                u_beta_pi = _u_pion_from_beta(R0, beta)  # one-sigma from HV spec

                # ADD: charge-ratio term from Mh,Ml (one-sigma)
                u_R   = R0 * np.sqrt(EM_REL_STD**2 + EM_REL_STD**2)
                den  = np.where(np.abs(beta - R0) < 1e-12, np.nan, np.abs(beta - R0))  # guard R≈s
                dP_dR = (beta - 1.0) / (den**2)  # = (s-1)/(s-R)^2
                u_pion_from_R = np.abs(dP_dR) * u_R
                
                # Total (one-sigma) → plot ~95% as ±2σ
                u_total = np.sqrt(u_beta_pi**2 + u_pion_from_R**2)
                yerr = 2.0 * u_total
           
                cont = ax1.errorbar(depths, y, yerr=yerr, fmt='-o', capsize=3, label=label)
                mainline = cont.lines[0] if hasattr(cont, "lines") else cont[0]
                main_lines_fig1.append(mainline); labels_fig1.append(label)
            
                scan_type = str(key[-1]).upper()

                curves_summary.append({
                    "x": depths, "mean": y, "yerr": yerr,
                    "label": label, "energy": energy_label, "scan_type": scan_type
                })

            else:
                mean = pion_arr.mean(axis=0)
                mean = np.maximum(mean, 1.0)  # clamp
                std  = pion_arr.std(axis=0, ddof=1)  # replicate scatter (one-sigma)
            
                # HV contribution, evaluated on mean charges at each depth
                mL_stack = np.vstack([a[:,1] for a in aligned_lows])   # (n_pairs, n_depths)
                mH_stack = np.vstack([a[:,1] for a in aligned_highs])  # (n_pairs, n_depths)
                mL_mean  = np.nanmean(mL_stack, axis=0)
                mH_mean  = np.nanmean(mH_stack, axis=0)
            
                beta = v_high / v_low
                with np.errstate(divide='ignore', invalid='ignore'):
                    R0 = np.divide(mH_mean, mL_mean, out=np.full_like(mH_mean, np.nan), where=(mL_mean!=0))
                u_beta_pi = _u_pion_from_beta(R0, beta)  # one-sigma from HV spec

                # ADD: charge-ratio term from Mh,Ml (one-sigma)
                u_R   = R0 * np.sqrt(EM_REL_STD**2 + EM_REL_STD**2)
                den  = np.where(np.abs(beta - R0) < 1e-12, np.nan, np.abs(beta - R0))
                dP_dR = (beta - 1.0) / (den**2)  # = (s-1)/(s-R)^2
                u_pion_from_R = np.abs(dP_dR) * u_R
                
                # Combine uncertainties in quadrature; plot 95% (~±2σ)
                u_total = np.sqrt(std**2 + u_beta_pi**2 + u_pion_from_R**2)
                yerr = 2.0 * u_total

            
                cont = ax1.errorbar(depths, mean, yerr=yerr, fmt='-o', capsize=3, label=label)
                mainline = cont.lines[0] if hasattr(cont, "lines") else cont[0]
                main_lines_fig1.append(mainline); labels_fig1.append(label)
            
                scan_type = str(key[-1]).upper() if len(key) > len(group_facets) else "PDD"
                curves_summary.append({
                    "x": depths, "mean": mean, "yerr": yerr,
                    "label": label, "energy": energy_label, "scan_type": scan_type
                })




            # Excel rows
            for pi in range(pion_arr.shape[0]):
                d_aln  = aligned_lows[pi][:,0]
                yL_aln = aligned_lows[pi][:,1]
                yH_aln = aligned_highs[pi][:,1]
                pion_k = np.maximum(pion_arr[pi], 1.0)  # clamp
                for j, depth in enumerate(d_aln):
                    rows.append({"Energy": energy_label,"Detector": det_name,"SSD_cm": ssd_cm,
                                 "HV_low_V": float(v_low),"HV_high_V": float(v_high),"Pair": pi+1,
                                 "Depth_cm": float(depth),"Voltage_V": float(v_low),"Charge_C": float(yL_aln[j]),"Pion": float(pion_k[j])})
                    rows.append({"Energy": energy_label,"Detector": det_name,"SSD_cm": ssd_cm,
                                 "HV_low_V": float(v_low),"HV_high_V": float(v_high),"Pair": pi+1,
                                 "Depth_cm": float(depth),"Voltage_V": float(v_high),"Charge_C": float(yH_aln[j]),"Pion": float(pion_k[j])})

    if not rows:
        print("No plottable groups found.")
        return

    ax1.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    # ---- Axis labels and title based on scan type ----
    # Use the scan type of the first plotted group as representative
    # ---- Axis labels and title based on scan type ----
    first_key = next(iter(chosen_keys))
    scan_type = ""
    if len(first_key) > len(group_facets):
        scan_type = str(first_key[-1]).upper()
    
    if "PROFILE" in scan_type:
        ax1.set_xlabel("Off-axis Position [cm]")
        ax1.set_title("Pion vs Off-axis Position")
    elif scan_type == "PDD":
        ax1.set_xlabel("Depth [cm]")
        ax1.set_title("Pion vs Depth")
    else:
        ax1.set_xlabel("Depth / Off-axis Position [cm]")
        ax1.set_title("Pion vs Depth / Off-axis Position")

    
    ax1.set_ylabel("Pion")
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0,30)
    ax1.set_ylim(.98, 1.06)
    _attach_click_legend(fig1, main_lines_fig1, labels_fig1)
    fig1.tight_layout()



    # Figure 2: CF (normalized differently for PDD vs Profile)
    fig2 = plt.figure(figsize=(7.6, 5.2))
    ax2 = fig2.add_subplot(111)
    main_lines_fig2, labels_fig2 = [], []

    for item in curves_summary:
        #print("Unique scan types in curves_summary:", {i.get("scan_type") for i in curves_summary})

        x   = item["x"]
        y   = item["mean"]
        err = item["yerr"]
        lbl = item["label"]
        en  = item.get("energy")
        scan_type = ""
        if "scan_type" in item:
            scan_type = str(item["scan_type"]).upper()
        elif "PROFILE" in str(item.get("label", "")).upper():
            scan_type = "INPLANE_PROFILE"
        else:
            scan_type = "PDD"


        if y is None or y.size < 2:
            continue

        if "PROFILE" in scan_type:
            # Normalize at CAX (x=0)
            x0 = float(min(max(0.0, np.nanmin(x)), np.nanmax(x)))
            denom = np.interp(x0, x, y)
        else:
            # PDD: normalize at fixed depth
            d = NORM_DEPTH_CM.get(en)
            if d is None:
                print(f"[CF norm] No depth found for energy '{en}'. Skipping.")
                continue
            d = float(min(max(d, float(np.nanmin(x))), float(np.nanmax(x))))
            denom = np.interp(d, x, y)

        if not np.isfinite(denom) or denom == 0:
            print(f"[CF norm] Bad normalization point for '{lbl}'. Skipping.")
            continue

        yn = y / denom
        if err is not None:
            err_n = err / denom
            cont = ax2.errorbar(x, yn, yerr=err_n, fmt='-o', capsize=3, label=lbl)
            mainline = cont.lines[0] if hasattr(cont, "lines") else cont[0]
        else:
            mainline, = ax2.plot(x, yn, marker='o', label=lbl)
        main_lines_fig2.append(mainline)
        labels_fig2.append(lbl)

    # Axis labeling based on scan type (first curve)
    if curves_summary:
        st0 = str(curves_summary[0].get("scan_type", "PDD")).upper()
    else:
        st0 = "PDD"

    if "PROFILE" in st0:
        ax2.set_xlabel("Position [cm]")
        ax2.set_title("Profile CF (normalized at CAX)")
    else:
        ax2.set_xlabel("Depth [cm]")
        ax2.set_title("PDD CF (normalized at fixed depth per energy)")

    ax2.set_ylabel("CF (normalized)")
    ax2.set_xlim(0, 30)
    ax2.set_ylim(.95, 1.02)
    ax2.grid(True, alpha=0.3)
    _attach_click_legend(fig2, main_lines_fig2, labels_fig2)
    fig2.tight_layout()


    plt.show()

    # Save for Excel
    last_df = pd.DataFrame(rows, columns=[
        "Energy","Detector","SSD_cm","HV_low_V","HV_high_V",
        "Pair","Depth_cm","Voltage_V","Charge_C","Pion"
    ])

    # Meta that matches current facets
    last_meta = {
        "GroupFacets": " / ".join(group_facets),
        "Num_groups": len(groups),
        "HV_pairs": ", ".join(sorted({f"{int(r['HV_high_V'])}/{int(r['HV_low_V'])}" for r in rows})),
    }
    facet_names = {"Energy":"Energies","Detector":"Detectors","SSD_cm":"SSDs_cm","FieldSize_cm":"FieldSizes_cm"}
    for facet in group_facets:
        idx = group_facets.index(facet)
        vals = { k[idx] for k in groups.keys() }
        if facet in ("SSD_cm","FieldSize_cm"):
            vals_fmt = sorted({ "?" if (v is None or not np.isfinite(v)) else f"{v:g} cm" for v in vals })
        else:
            vals_fmt = sorted({ str(v) if v not in (None,"") else "?" for v in vals })
        last_meta[facet_names.get(facet, facet+"s")] = ", ".join(vals_fmt)

def export_excel():
    if last_df is None:
        print("Run Plot Pion first to generate results.")
        return
    save_path = filedialog.asksaveasfilename(
        defaultextension=".xlsx",
        initialfile="pion_results.xlsx",
        filetypes=[("Excel", "*.xlsx")]
    )
    if not save_path:
        return
    with pd.ExcelWriter(save_path, engine="xlsxwriter") as xw:
        last_df.to_excel(xw, index=False, sheet_name="Pion")
        pd.DataFrame([last_meta]).to_excel(xw, index=False, sheet_name="Meta")
    print(f"Saved: {save_path}")
    try:
        if sys.platform.startswith("win"):
            os.startfile(save_path)
        elif sys.platform == "darwin":
            subprocess.run(["open", save_path], check=False)
        else:
            subprocess.run(["xdg-open", save_path], check=False)
    except Exception as e:
        print(f"Saved: {save_path} (auto-open failed: {e})")
    else:
        print(f"Saved and opened: {save_path}")

# ------------------------------
# GUI setup
# ------------------------------
root = tk.Tk()
root.title("Pion vs Depth")

def on_close():
    root.destroy()
    root.quit()

root.protocol("WM_DELETE_WINDOW", on_close)

top = ttk.Frame(root, padding=10)
top.grid(row=0, column=0, sticky="ew")
root.columnconfigure(0, weight=1)
top.columnconfigure(1, weight=1)

ttk.Label(top, text="MCC File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(top, width=50)
file_entry.grid(row=0, column=1, padx=5, sticky="ew")
ttk.Button(top, text="Browse", command=browse_file).grid(row=0, column=2, padx=5)
ttk.Button(top, text="Browse Folder", command=browse_folder).grid(row=0, column=3, padx=5)

sel_frame = ttk.Frame(root, padding=(10, 0, 10, 10))
sel_frame.grid(row=1, column=0, sticky="nsew")
root.rowconfigure(1, weight=1)
sel_frame.rowconfigure(1, weight=1)
sel_frame.columnconfigure(0, weight=1)

ttk.Label(sel_frame, text="Select Groups to Plot:").grid(row=0, column=0, sticky="w")

lb_frame = ttk.Frame(sel_frame)
lb_frame.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
lb_frame.rowconfigure(0, weight=1)
lb_frame.columnconfigure(0, weight=1)

group_listbox = tk.Listbox(lb_frame, selectmode="extended", height=8)
lb_scroll = ttk.Scrollbar(lb_frame, orient="vertical", command=group_listbox.yview)
group_listbox.configure(yscrollcommand=lb_scroll.set)

group_listbox.grid(row=0, column=0, sticky="nsew")
lb_scroll.grid(row=0, column=1, sticky="ns")

btns = ttk.Frame(sel_frame)
btns.grid(row=2, column=0, sticky="w", pady=6)

def select_all_groups():
    group_listbox.select_set(0, tk.END)

def clear_groups():
    group_listbox.selection_clear(0, tk.END)

ttk.Button(btns, text="Select All", command=select_all_groups).pack(side="left", padx=(0,6))
ttk.Button(btns, text="Clear",      command=clear_groups).pack(side="left")

action_frame = ttk.Frame(root, padding=10)
action_frame.grid(row=2, column=0, sticky="ew")
ttk.Button(action_frame, text="Plot Pion", command=plot_pion).pack(side="left", padx=5)
ttk.Button(action_frame, text="Export Excel", command=export_excel).pack(side="left", padx=5)
ttk.Button(action_frame, text="Grouping…", command=open_grouping_dialog).pack(side="left", padx=5)
only_beyond_dmax_var = tk.BooleanVar(value=False)
ttk.Checkbutton(action_frame, text="Beyond dmax only", variable=only_beyond_dmax_var).pack(side="left", padx=12)
root.mainloop()
