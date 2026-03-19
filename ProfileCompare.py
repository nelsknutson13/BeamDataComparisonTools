#Created by nels knutson. Current code does not report total pass rate for COMP analysis correctly

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import interpolate as interp
from scipy import signal as sig
from gamma import gamma as g
from center import center as center_function
from comp import dosedif as dosedif
from comp import dta as dtafunc
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
import tkinter as tk
from tkinter import filedialog, ttk
from time import perf_counter
import re
from collections import defaultdict  # if not already imported
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import os, platform, subprocess
import io, contextlib  # for capturing text into report_text

xf = None;sheet_cache={}
fig = None
report_text = ""
DEPTH_ROUND_CM = 0.1   # 1 mm


# Function to handle file selection and populate sheet names
class PassAgg:
    def __init__(self, name=""):
        self.name = name
        self.total_pass = 0
        self.total      = 0
        self.by_fs      = defaultdict(lambda: [0, 0])      # FS -> [pass, total]
        self.by_depth   = defaultdict(lambda: [0, 0])      # depth -> [pass, total]
        self.matrix     = defaultdict(lambda: [0, 0])      # (FS, depth) -> [pass, total]
        # Overall region totals: [pass, total]
        self.region_overall = {
            "in":  [0, 0],
            "pen": [0, 0],
            "tail":[0, 0],
            "ovr_hi":  [0, 0],   
            "ovr_lo":  [0, 0],
            }
    @staticmethod
    def _rate(p, t): 
        return (100.0 * p / t) if t else float('nan')

    def add(self, fs, depth, passed, total):
        # fs, depth should be numbers (float ok)
        p, t = int(passed), int(total)
        self.total_pass += p; self.total += t
        self.by_fs[float(fs)][0]    += p; self.by_fs[float(fs)][1]    += t
        self.by_depth[float(depth)][0] += p; self.by_depth[float(depth)][1] += t
        self.matrix[(float(fs), float(depth))][0] += p
        self.matrix[(float(fs), float(depth))][1] += t
        

    def add_region_totals(self, in_pass, in_tot, pen_pass, pen_tot, tail_pass, tail_tot,
                      ovr_hi_pass=0, ovr_hi_tot=0, ovr_lo_pass=0, ovr_lo_tot=0):

        # core regions
        self.region_overall["in"][0]   += int(in_pass)
        self.region_overall["in"][1]   += int(in_tot)
        self.region_overall["pen"][0]  += int(pen_pass)
        self.region_overall["pen"][1]  += int(pen_tot)
        self.region_overall["tail"][0] += int(tail_pass)
        self.region_overall["tail"][1] += int(tail_tot)
    
        # overlaps (these lines were missing)
        self.region_overall["ovr_hi"][0] += int(ovr_hi_pass)
        self.region_overall["ovr_hi"][1] += int(ovr_hi_tot)
        self.region_overall["ovr_lo"][0] += int(ovr_lo_pass)
        self.region_overall["ovr_lo"][1] += int(ovr_lo_tot)


    def build_matrix(self):
        fs_sorted    = sorted(self.by_fs.keys())
        depth_sorted = sorted(self.by_depth.keys())
        import numpy as _np
        mat = _np.full((len(depth_sorted), len(fs_sorted)), _np.nan)
        for r, d in enumerate(depth_sorted):
            for c, fs in enumerate(fs_sorted):
                p, t = self.matrix[(fs, d)]
                mat[r, c] = self._rate(p, t)
        # row/col totals
        row_totals = []
        for d in depth_sorted:
            p_row = sum(self.matrix[(fs, d)][0] for fs in fs_sorted)
            t_row = sum(self.matrix[(fs, d)][1] for fs in fs_sorted)
            row_totals.append(self._rate(p_row, t_row))
        col_totals = []
        for fs in fs_sorted:
            p_col = sum(self.matrix[(fs, d)][0] for d in depth_sorted)
            t_col = sum(self.matrix[(fs, d)][1] for d in depth_sorted)
            col_totals.append(self._rate(p_col, t_col))
        overall = self._rate(self.total_pass, self.total)
        return fs_sorted, depth_sorted, mat, row_totals, col_totals, overall

    def print_summary(self):
        fs_sorted, depth_sorted, mat, row_totals, col_totals, overall = self.build_matrix()
        fail = self.total - self.total_pass
        # NEW: energy + SSD first
        if hasattr(self, "energy") or hasattr(self, "ssd_cm"):
            parts = []
            if hasattr(self, "energy"):
                parts.append(f"Energy: {self.energy}")
            if hasattr(self, "ssd_cm"):
                parts.append(f"SSD: {self.ssd_cm:.1f} cm SSD")
            print("\n" + " | ".join(parts))
    
        print(f"=== {self.name.upper()} PASS RATE SUMMARY ===")
        print(f"\n=== {self.name.upper()} PASS RATE SUMMARY ===")
    
        # Optional sheet names (if you set self.sheets = (sn1, sn2))
        if hasattr(self, "sheets") and self.sheets:
            s1, s2 = self.sheets
            print(f"Sheets compared: {s1}  vs  {s2}")
    
        # Optional criteria display
        crit = {}
        if hasattr(self, "criteria") and isinstance(getattr(self, "criteria"), dict):
            crit = dict(self.criteria)
        else:
            if hasattr(self, "crit_gamma"):
                crit["gamma"] = self.crit_gamma
            if hasattr(self, "crit_comp"):
                crit["comp"] = self.crit_comp
            if hasattr(self, "crit_mppg"):
                crit["mppg"] = self.crit_mppg
    
        def _fmt_gamma(t):
            try:
                dd_pct, dta_mm = t
                return f"{dd_pct:.1f}% / {dta_mm:.1f}mm"
            except Exception:
                return str(t)
    
        def _fmt_mppg(t):
            try:
                dd_pct, dta_mm, ddtail_pct = t
                return f"{dd_pct:.1f}% / {dta_mm:.1f}mm, tails {ddtail_pct:.1f}% of Dmax"
            except Exception:
                return str(t)
    
        crit_lines = []
        if "gamma" in crit: crit_lines.append(f"Gamma criteria: {_fmt_gamma(crit['gamma'])}")
        if "comp"  in crit: crit_lines.append(f"Composite criteria: {_fmt_gamma(crit['comp'])}")
        if "mppg"  in crit: crit_lines.append(f"MPPG criteria: {_fmt_mppg(crit['mppg'])}")
        if crit_lines:
            print("Criteria:", " | ".join(crit_lines))
    
        print(f"Overall pass rate: {overall:.2f}%  (n={self.total})")
    
        print("\nPer-Field Size pass rate (%):")
        for fs in fs_sorted:
            p, t = self.by_fs[fs]
            rate = self._rate(p, t)
            print(f"  FS {fs:g} cm : {rate:.2f}%  (n={t})")
    
        print("\nPer-Depth pass rate (%):")
        for d in depth_sorted:
            p, t = self.by_depth[d]
            rate = self._rate(p, t)
            print(f"  Depth {d:g} cm : {rate:.2f}%  (n={t})")
    
        # Matrix
        hdr = "Depth\\FS".ljust(10) + "".join(f"{fs:>8g}" for fs in fs_sorted) + "   |  RowTot"
        print("\nDepth × FS matrix (pass rate %):")
        print(hdr)
        for r, d in enumerate(depth_sorted):
            row_vals = "".join(f"{mat[r,c]:8.2f}" for c in range(len(fs_sorted)))
            print(f"{d:>10g}{row_vals}   |  {row_totals[r]:6.2f}")
        print(" " * 10 + "-" * (8*len(fs_sorted) + 9))
        print(" ColTot   " + "".join(f"{ct:8.2f}" for ct in col_totals) + f"   |  {overall:6.2f}")
        # --- Region breakdown (MPPG only) ---
        if hasattr(self, "region_overall"):
            ro = self.region_overall
        
            def _fmt(p, t):
                if t:
                    rate = self._rate(p, t)
                    return f"{rate:6.2f}%  (pass={p}/{t}, fail={t - p})"
                else:
                    return "n/a"
        
            print("\nRegion breakdown (MPPG rules):")
            print(f"  In-field (DD)              : {_fmt(*ro['in'])}")
            print(f"  Penumbra (DTA)             : {_fmt(*ro['pen'])}")
            print(f"  Tails (DD as % of Dmax)    : {_fmt(*ro['tail'])}")
            print(f"  Overlap (in-field side)    : {_fmt(*ro['ovr_hi'])}")
            print(f"  Overlap (tail side)        : {_fmt(*ro['ovr_lo'])}")
        
            # Optional: consistency check (regions should sum to overall counted points)
            sum_pass = sum(ro[k][0] for k in ro)
            sum_tot  = sum(ro[k][1] for k in ro)
            if sum_tot != self.total or sum_pass > self.total_pass:
                print(f"  [note] region totals (pass={sum_pass}, n={sum_tot}) "
                      f"!= overall counted (pass={self.total_pass}, n={self.total}).")

        # Footer (now uses defined vars)
        print(f"\n=== {self.name.upper()} PASS RATE SUMMARY ===")
        if hasattr(self, "sheets") and self.sheets:
            s1, s2 = self.sheets
            print(f"Sheets compared: {s1}  vs  {s2}")
        if crit_lines:
            print("Criteria:", " | ".join(crit_lines))
        print(f"Overall pass rate: {overall:.2f}%  (n={self.total}), fails={fail}")
        print(f"=== END {self.name.upper()} SUMMARY ===\n")




def downsample_to_native(x_interp, vals_interp, x_native):
    """
    Downsample interpolated data (x_interp, vals_interp)
    to roughly match the spacing of x_native.

    Enhancement:
      - If x_native is a list/tuple of arrays (e.g., [y1, y2]),
        it matches the COARSER (larger) median step size.
    """
    x_interp = np.asarray(x_interp, dtype=float)
    vals_interp = np.asarray(vals_interp)

    # accept x_native as single array OR list/tuple of arrays
    if isinstance(x_native, (list, tuple)):
        natives = [np.asarray(x, dtype=float) for x in x_native if x is not None]
    else:
        natives = [np.asarray(x_native, dtype=float)]

    if x_interp.size < 2 or len(natives) == 0:
        return x_interp, vals_interp

    def _median_step(x):
        x = x[np.isfinite(x)]
        if x.size < 2:
            return 0.0
        xu = np.unique(x); xu.sort()
        diffs = np.diff(xu)
        diffs = diffs[diffs > 0]
        return float(np.median(diffs)) if diffs.size else 0.0

    steps = [s for s in (_median_step(x) for x in natives) if s > 0]
    if not steps:
        return x_interp, vals_interp

    native_step = max(steps)  # coarser = larger step

    xi = np.unique(x_interp[np.isfinite(x_interp)])
    xi.sort()
    if xi.size < 2:
        return x_interp, vals_interp
    interp_step = float(np.median(np.diff(xi)))
    if interp_step <= 0:
        return x_interp, vals_interp

    factor = max(1, int(round(native_step / interp_step)))
    return x_interp[::factor], vals_interp[::factor]




def select_file():
    global xf, sheet_cache, SSD_CM, ENERGY
    file_path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
    if not file_path:
        print("No file selected.")
        return
        
    SSD_CM, ENERGY = parse_ssd_energy(file_path, default_ssd=100.0, default_energy="6X")
    print(f"Context parsed from path → SSD: {SSD_CM} cm | Energy: {ENERGY}")

    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)

    try:
        xf = pd.ExcelFile(file_path, engine="openpyxl")  # or "calamine" when you switch
        sheets = xf.sheet_names
        print(f"Selected File: {file_path}")
        print(f"Available Sheets: {sheets}")

        sheet1_combo['values'] = sheets
        sheet2_combo['values'] = sheets
        sheet1_combo.set(sheets[0])
        sheet2_combo.set(sheets[1] if len(sheets) > 1 else sheets[0])
        print(f"Sheet 1 selected: {sheet1_combo.get()}")
        print(f"Sheet 2 selected: {sheet2_combo.get()}")

        # (Re)build cache for the two shown sheets, then populate lists
        sheet_cache.clear()
        preload_two_sheets()

        def on_combobox_change(event):
            print(f"Sheet 1 selected: {sheet1_combo.get()}")
            print(f"Sheet 2 selected: {sheet2_combo.get()}")
            preload_two_sheets()
            populate_fsl()
            populate_scl()
            populate_depths()

        sheet1_combo.bind("<<ComboboxSelected>>", on_combobox_change)
        sheet2_combo.bind("<<ComboboxSelected>>", on_combobox_change)

        populate_fsl()
        populate_scl()
        populate_depths()
        print("Ready: field sizes, scan types, and depths populated.")

    except Exception as e:
        print(f"Error reading sheets: {e}")
def parse_ssd_energy(path: str, default_ssd=100.0, default_energy="6X"):
    """
    Extract SSD (cm) and Energy token from the file path, tolerating
    underscores, dashes, and spaces (e.g., '6FFF_90SSD_*.xlsx', '...-10X-...').
    """
    s = str(path)

    # --- SSD: handle 'SSD100', 'SSD_100', '100SSD', 'SSD100cm'
    m_ssd = (re.search(r'ssd[_\-\s]?(\d+(?:\.\d+)?)', s, re.I)
             or re.search(r'(\d+(?:\.\d+)?)\s*ssd', s, re.I))
    ssd_val = float(m_ssd.group(1)) if m_ssd else float(default_ssd)

    # --- ENERGY: normalize separators so word boundaries work
    s_norm = re.sub(r'[^A-Za-z0-9.]+', ' ', s)  # turn '_', '-', etc. into spaces
    m_energy = re.search(r'(?:^|\s)(6X|6FFF|8FFF|10X|10FFF|15X)(?=\s|$)', s_norm, re.I)
    energy_val = m_energy.group(1).upper() if m_energy else default_energy

    return ssd_val, energy_val
def preload_two_sheets():
    """Read only FS, Axis, Depth from each selected sheet ONCE and cache them."""
    global sheet_cache, xf
    usecols = ['FS', 'Axis', 'Depth']
    for sn in (sheet1_combo.get(), sheet2_combo.get()):
        if not sn:
            continue
        if sn in sheet_cache:
            continue
        t0 = perf_counter()
        # Try faster engine if available
        try:
            df = pd.read_excel(xf, sheet_name=sn, usecols=usecols, engine="pyarrow")
        except Exception:
            df = pd.read_excel(xf, sheet_name=sn, usecols=usecols)
        df["Depth"] = ((df["Depth"] / DEPTH_ROUND_CM).round() * DEPTH_ROUND_CM).round(3)



        
        sheet_cache[sn] = df
        #print(f"[T] Cache '{sn}' {usecols}: {perf_counter()-t0:.2f}s (rows={len(df)})")

# --- PDD lookup by depth (depth-only; no energy yet) ---
# --- PDD lookup table (SSD = 100 cm, normalized at dmax, 10.5cm field) ---
PDD_TABLE = {
    "6X":   {1.42:1.000,  5.0:0.865, 10.0:0.668, 20.0:0.382, 30.0:0.218},
    "6FFF": {1.32:1.000,  5.0:0.848, 10.0:0.636, 20.0:0.346, 30.0:0.190},
    "8FFF": {1.88:1.000,  5.0:0.891, 10.0:0.693, 20.0:0.406, 30.0:0.240},
    "10X":  {2.22:1.000,  5.0:0.925, 10.0:0.744, 20.0:0.469, 30.0:0.296},
    "10FFF":{2.16:1.000,  5.0:0.904, 10.0:0.708, 20.0:0.423, 30.0:0.256},
    "15X":  {2.72:1.000,  5.0:0.949, 10.0:0.774, 20.0:0.502, 30.0:0.325},
}

def pdd_lookup_nearest(energy, depth_cm):
    table = PDD_TABLE.get(energy, {})
    if not table:
        return 1.0, None  # no table for that energy

    keys = np.array(list(table.keys()), dtype=float)
    k = float(keys[np.argmin(np.abs(keys - float(depth_cm)))])
    return float(table[k]), k


def populate_fsl():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        df2 = sheet_cache[sheet2_combo.get()]
        common_fsl = sorted(set(df1['FS'].dropna().unique()).intersection(df2['FS'].dropna().unique()))

        fsl_listbox.delete(0, tk.END)
        for v in common_fsl:
            fsl_listbox.insert(tk.END, v)
        fsl_listbox.selection_clear(0, tk.END)   # default: no selection

        
        print("Common Field Sizes:", [float(v) for v in common_fsl])

    except Exception as e:
        print(f"Error populating field sizes: {e}")


def populate_scl():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        df2 = sheet_cache[sheet2_combo.get()]
        common_scl = sorted(set(df1['Axis'].dropna().unique()).intersection(df2['Axis'].dropna().unique()))

        scl_listbox.delete(0, tk.END)
        for s in common_scl:
            scl_listbox.insert(tk.END, s)
        scl_listbox.selection_clear(0, tk.END)   # default: no selection

        print("Common Scan Types:", common_scl)
    except Exception as e:
        print(f"Error loading scan types: {e}")


def populate_depths():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        df2 = sheet_cache[sheet2_combo.get()]
        # >>> ADD THESE TWO PRINTS <<<
        #print("Depths sheet 1:", sorted(df1['Depth'].dropna().unique()))
        #print("Depths sheet 2:", sorted(df2['Depth'].dropna().unique()))
      # >>> END ADD <<<
        common_depths = sorted(set(df1['Depth'].dropna().unique()).intersection(df2['Depth'].dropna().unique()))

        
        depth_listbox.delete(0, tk.END)
        for v in common_depths:
            depth_listbox.insert(tk.END, v)
        depth_listbox.selection_clear(0, tk.END) # default: no selection

        print("Common Depths:", common_depths)
    except Exception as e:
        print(f"Error populating depths: {e}")

fsl=[];scl=[];dl=[];
 


#Plotting
# Function to run the comparison based on the selected sheets
def run_comparison():
    global fig, report_text

    # Retrieve dd and dta values from the input fields
    try:
        msize = float(marker_size_entry.get())
    except ValueError:
        msize = 2
    dd = float(dd_entry.get())/100
    ddtail=float(ddtail_entry.get())/100
    dta = float(dta_entry.get())/10
    print(f"Selected Dose Difference: {dd} ")
    print(f"Selected DTA Criteria: {dta} cm")
    crit_gamma = f"{dd*100:.1f}%/{dta*10:.1f}mm"
    crit_comp  = f"{dd*100:.1f}%/{dta*10:.1f}mm"
    crit_mppg  = f"{dd*100:.1f}% in-field / {dta*10:.1f}mm penumbra / {ddtail*100:.1f}% of dmax tails"
    agg= None
    global dl    
    global fsl
    global scl
    # Retrieve selected depths
    selected_depth_indices = depth_listbox.curselection()
    dl.clear()  # Clear the current depth list to avoid appending old values
    dl.extend([depth_listbox.get(i) for i in selected_depth_indices])
    dl = [float(x) for x in dl]  # Convert to float if necessary
    # Build PDD list aligned to the order of dl
    # NEW (nearest)
    PDD = []
    for d in dl:
        pdd, d_used = pdd_lookup_nearest(ENERGY, d)
        PDD.append(pdd)
        if d_used is not None and abs(d_used - d) > 1e-9:
            print(f"PDD lookup: depth {d:.2f} cm → using nearest {d_used:.2f} cm for {ENERGY}")


    selected_indices = fsl_listbox.curselection()
    fsl.clear()  # Clear the current fsl list to avoid appending old values
    fsl.extend([fsl_listbox.get(i) for i in selected_indices])
    fsl = [float(x) for x in fsl]
    sclindexs=scl_listbox.curselection()
    scl.clear()
    scl.extend([scl_listbox.get(i) for i in sclindexs])
    # Reorder scl to ensure XY is processed last
    scl = sorted(scl, key=lambda x: (x == 'XY', x))  # This moves 'XY' to the end
    
    mppg=0;
    
    #all of these set in gui intializing as 0
    scale=0;#Scales field size for offset scans
    dif=0;dist=0;comp=0;gam=0;mppg=0;plot=0#iniates as unselected selected by gui
    smooth=0;cent=0;norm =0;#Assumes norm to CAX set to 2 scale by PDD#Center for all measured profiles
       
    #still hard coded
    save=0;
    
    #gs= gridspec.GridSpec(1,1)   
    
    fn = file_entry.get()
    sn1 = sheet1_combo.get()
    sn2 = sheet2_combo.get()
    selected_analysis = analysis_var.get()
    # Debug: Print the selected file and sheets
    print(f"Selected File: {fn}")
    print(f"Sheet 1: {sn1}")
    print(f"Sheet 2: {sn2}")
    # Set the selected analysis variable to 1
    if selected_analysis == "gam":
        gam = 1
    elif selected_analysis == "comp":
        comp = 1
    elif selected_analysis == "dif":
        dif = 1
    elif selected_analysis == "dist":
        dist = 1
    elif selected_analysis == "mppg":
        mppg = 1
    elif selected_analysis == "plot":
        plot = 1
    # Set data processing variables based on user selection
    cent = int(cent_combo.get())
    smooth = int(smooth_combo.get())
    conv= int(conv_combo.get())
    scale = scale_var.get()
    norm = int(norm_combo.get())
    
    
   
    if mppg ==1:
        plt.rcParams.update({'font.size': 22})
        gs = gridspec.GridSpec(4,1,height_ratios=[1.5,1,1,1])
        fig, (ax0, ax1, ax2, ax3) = plt.subplots(
            4, 1, figsize=(15, 12),
            gridspec_kw={'height_ratios': [1.5, 1, 1, 1]}
        )
    elif comp ==1 or gam==1:
        plt.rcParams.update({'font.size': 20})
        fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 11),
            gridspec_kw={'height_ratios': [1.5, 1, 1]})
    elif dif==1 or dist==1:
        plt.rcParams.update({'font.size': 18})
        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15, 9),
            gridspec_kw={'height_ratios': [1.5, 1]})
    elif plot == 1:
        plt.rcParams.update({'font.size': 16})
        fig, ax0 = plt.subplots(1, 1, figsize=(15, 8))
    else:
        plt.rcParams.update({'font.size': 16})
        gs = gridspec.GridSpec(1, 1)
        fig, ax0 = plt.subplots(1, 1, figsize=(10, 5))
       
    #ax0=plt.subplot(gs[0])
    ax0.plot([],'+r',ms=10,label=str(sn1))
    ax0.plot([],'.k',ms=10,label=sn2)
    ax0.legend()
    ax0.set_ylabel('Off Axis Ratio');
    #ax0.set_title('1DS Scanning Simulation')
    ax0.set_xlabel(' Off Axis Position [cm]')
    #cmap=[];
    if gam == 1:
        agg = PassAgg(f"Gamma ({crit_gamma})")
        agg.sheets = (sn1, sn2)
        agg.energy = ENERGY
        agg.ssd_cm = SSD_CM
    
    if comp == 1:
        agg = PassAgg(f"Composite ({crit_comp})")
        agg.sheets = (sn1, sn2)
        agg.energy = ENERGY
        agg.ssd_cm = SSD_CM
    
    if mppg == 1:
        agg_mppg = PassAgg(f"MPPG ({crit_mppg})")
        agg_mppg.sheets = (sn1, sn2)
        agg_mppg.energy = ENERGY
        agg_mppg.ssd_cm = SSD_CM

        # give print_summary a machine-readable tuple, it will format as: "dd% / dta mm, tails ddtail% of Dmax"
        agg_mppg.crit_mppg = (dd*100, dta*10, ddtail*100)    
    prtot='n/a';gvmax='n/a';gvmean='n/a'
    centl=[.35,.5,.5,.5,.5,.4,.3]#Where to center around
    prl=[];gvl=[];gvtot=[];
    sf=1;s1=0;s2=0 
    end=2;
    blist=.0000924# 30cm depth coeficcent
    tf=0;tot=0;totpf=0;failpf=0;tf_dif=0;tot_dif=0;tf_dist=0;tot_dist=0
    dif_results=[];dist_results=[]
    print("Comparison started...")
    df1 = pd.read_excel(fn, sheet_name=sn1, header=0)
    df2 = pd.read_excel(fn, sheet_name=sn2, header=0)
    df1["Depth"] = ((df1["Depth"] / DEPTH_ROUND_CM).round() * DEPTH_ROUND_CM).round(3)
    df2["Depth"] = ((df2["Depth"] / DEPTH_ROUND_CM).round() * DEPTH_ROUND_CM).round(3)
    df1=df1.sort_values(by=['FS','Depth','Axis','Pos'])
    df2=df2.sort_values(by=['FS','Depth','Axis','Pos'])
    # Before starting the loop, check the type and contents of fsl
    print("Type of fsl:", type(fsl))
    print("Contents of fsl:", fsl)
    # Example of using these in your loop
    for j in range(len(fsl)):  # Your existing loop over field sizes
        prfs=[];gvfs=[];totpf=0;failpf=0;
        print(f"Processing Field Size: {fsl[j]}")

        for k in range(len(scl)):  # Your existing loop over scan types
            # Handle special case for 'XY' only on the largest field size
            if (scl[k] == 'XY' or scl[k] == 'YX') and fsl[j] != max(fsl):
                continue  # Skip processing if 'XY' or 'YX' is not on the largest field size

            
            print(f"Processing Scan Type: {scl[k]}")
                       
            for i in range(len(dl)):#depth list 5 max for 30cm# added a sort by Depth
                print(f"Processing Depth: {dl[i]}")

                y1=df1.loc[(df1['Depth']==dl[i]) & (df1['FS']==fsl[j])&(df1['Axis']==scl[k]),'Pos'];
                d1=df1.loc[(df1['Depth']==dl[i])& (df1['FS']==fsl[j])&(df1['Axis']==scl[k]),'Dose']
                y2=df2.loc[(df2['Depth']==dl[i])&(df2['FS']==fsl[j])&(df2['Axis']==scl[k]),'Pos'];
                d2=df2.loc[(df2['Depth']==dl[i])&(df2['FS']==fsl[j])&(df2['Axis']==scl[k]),'Dose']
                
                if cent==1:
                    s1=center_function(y1,d1,0.3) 
                    s2=0
                    #print(s1,s2)
                if cent==2:
                    s1=0
                    s2=center_function(y2,d2,0.3)
                    #print(s1,s2)
                if cent ==3:
                    s1=center_function(y1,d1,.3);s2=center_function(y2,d2,.3)
                    #s2=s2-.05
                    #print(s1,s2)
                    #print(f"Analyzing {scl[k]} Scan {fsl[j]} Field Size Depth {dl[i]} Shift 1 = {s1:.2f} Shift 2 = {s2:.2f}")
                if cent in (1, 3):
                    y1 = y1 + s1
                if cent in (2, 3):
                    y2 = y2 + s2

                
                
                detector_diameter = float(det_diam_var.get())
                y2_spacing =  np.mean(np.diff(y2))
                sigma_physical = detector_diameter / (2 * np.sqrt(2 * np.log(2)))
                sigma_data_points = sigma_physical / y2_spacing
                
                if conv ==1:
                    
                    # Convolve dose data directly with a Gaussian filter
                    d1 = gaussian_filter1d(d1, sigma=sigma_data_points)
                    
                if conv == 2:
                                        
                    # Apply Gaussian convolution to simulate the detector response
                    d2 = gaussian_filter1d(d2, sigma=sigma_data_points)
                    

                if conv ==3:
                                        
                    d1 = gaussian_filter1d(d1, sigma=sigma_data_points)
                    d2 = gaussian_filter1d(d2, sigma=sigma_data_points)
                if smooth == 1:
                    d1=sig.savgol_filter(d1,3,1)
                if smooth == 2:
                    d2=sig.savgol_filter(d2,3,1)
                if smooth == 3:
                    d1=sig.savgol_filter(d1,3,1)
                    d2=sig.savgol_filter(d2,3,1)
             
                
             
                
                if norm ==1:# Normalize to +/- 3mm of CAX
                    d1=d1/df1.loc[(df1['Depth']==dl[i])&(np.abs(df1['Pos'])<=0.3) & (df1['FS']==fsl[j])&(df1['Axis']==scl[k]),'Dose'].mean()
                    d2=d2/df2.loc[(df2['Depth']==dl[i])&(np.abs(df2['Pos'])<=0.3)&(df2['FS']==fsl[j])&(df2['Axis']==scl[k]),'Dose'].mean()
                if norm ==2:# Normalize to +/- 3mm of CAX and scale by PDD as function of depth.
                    d1=d1*PDD[i]/df1.loc[(df1['Depth']==dl[i])&(np.abs(df1['Pos'])<=0.3) & (df1['FS']==fsl[j])&(df1['Axis']==scl[k]),'Dose'].mean()
                    d2=d2*PDD[i]/df2.loc[(df2['Depth']==dl[i])&(np.abs(df2['Pos'])<=0.3)&(df2['FS']==fsl[j])&(df2['Axis']==scl[k]),'Dose'].mean()
                
                
                ax0.plot(y1,d1,'+r',ms=msize)
                ax0.plot(y2,d2,'.k',ms=msize)
    
                if gam ==1:
                    
                    gx , gv = g(y1,d1,y2,d2,dd,dta,1,0,.01);# norm 1 is to central 5 pixel,norm 2 is to dmax, 0 threshold InterpThreshold=0.01
                    # --- NEW: per-profile pass/total and add to aggregator ---
                    #down sample to native input
                    gx, gv = downsample_to_native(gx, gv, y1)
                    gv_a  = np.asarray(gv)
                    valid = np.isfinite(gv_a)
                    tot_p = int(valid.sum())
                    pas_p = int(np.count_nonzero(gv_a[valid] <= 1.0))
                    agg.add(fsl[j], dl[i], pas_p, tot_p)   # NEW

                    gx1=np.extract(gv>1,gx);gv1=np.extract(gv>1,gv);
                    gx2=np.extract(gv<1,gx);gv2=np.extract(gv<1,gv);
                    pr=(np.size((np.where(np.asarray(gv)<=1)))/np.size(np.where(np.asarray(gv)>=0))*100)
                    pr=np.round(pr,1)
                    #print(pr,max(np.round(gv,2)),fsl[j],dl[i],scl[k])
                    gvtot.extend(gv)
                    gvfs.extend(gv)
                    prfs=(np.size((np.where(np.asarray(gvfs)<=1)))/np.size(np.where(np.asarray(gvfs)>=0))*100)
                    prtot=(np.size((np.where(np.asarray(gvtot)<=1)))/np.size(np.where(np.asarray(gvtot)>=0))*100)
                    prtot=np.round(prtot,2)
                    prl.append(pr)
                    #ax0=plt.subplot(gs[0])
                    ax0.set_ylim(0,1.1)
                    ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    #ax1=fig.add_subplot(gs[1],sharex=ax0)
                    ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    #ax2=fig.add_subplot(gs[2])
                    ax1.plot(gx1[1::1],gv1[1::1],'.r',ms=msize)
                    ax1.plot(gx2[1::1],gv2[1::1],'.g',ms=msize)
                    #ax1.set_ylabel(r'$\Gamma$'+'['+ str(dd*100)+'%/'+ str(dta*10) + 'mm]')
                    #ax1.set_xlabel( 'Points with ' + r'$\Gamma$' +' > 1 : '+str(np.size((np.where(np.asarray(gvtot)>1))))+'/'+str(len(gvtot)) + '  Pass Rate : '+ str(prtot)+'%' )
                    gvmean=np.round(np.mean(gvtot),3);gvmax=np.round(np.max(gvtot),3); 
                    bins=0.1;
                    weights = np.ones_like(np.asarray(gvtot))/float(len(np.asarray(gvtot)))
                    
                    #plt.title('Gamma Histogram Data:  '  + str(np.round(100-sum(histd[0][0:(1/bins)]),1)) + '% of Data < Gamma = 1')
                    ax2.set_xlim(0,np.max(gvtot))
                    ax2.set_xticks(np.arange(0,np.max(gvtot),0.5))
                    
                    #print(pr,prtot,gvmean,gvmax,fsl[j],dl[i],scl[k])
    
                if dif ==1:
                    difx, ddif= dosedif(y1,d1,y2,d2,1)
                    difx, ddif = downsample_to_native(difx, ddif, y1)
                    ddx1=np.extract(np.abs(ddif)>dd,difx);ddv1=np.extract(np.abs(ddif)>dd,ddif);
                    ddx2=np.extract(np.abs(ddif)<dd,difx);ddv2=np.extract(np.abs(ddif)<dd,ddif);
                    tf_dif += len(ddx1)
                    tot_dif += len(ddif)
                    fs_dif_pr = 100 * (len(ddif) - len(ddx1)) / len(ddif) if len(ddif) > 0 else 0.0
                    dif_results.append((fsl[j], len(ddx1), len(ddif), fs_dif_pr))
                    ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    ax1.plot(ddx1,ddv1*100,'.r',ms=msize)
                    ax1.plot(ddx2,ddv2*100,'.g',ms=msize)
                    ax1.set_ylabel('Dose Difference [%]')
                if dist ==1:
                    dtax, dtav= dtafunc(y1,d1,y2,d2,dta)
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    dtax1=np.extract(dtav>dta,dtax);dtav1=np.extract(dtav>dta,dtav);
                    dtax2=np.extract(dtav<dta,dtax);dtav2=np.extract(dtav<dta,dtav);
                    tf_dist += len(dtax1)
                    tot_dist += len(dtav)
                    fs_dist_pr = 100 * (len(dtav) - len(dtax1)) / len(dtav) if len(dtav) > 0 else 0.0
                    dist_results.append((fsl[j], len(dtax1), len(dtav), fs_dist_pr))
                    #ax0=plt.subplot(gs[0])
                    ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    #ax1=fig.add_subplot(gs[1],sharex=ax0)
                    ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
                    ax1.plot(dtax1,dtav1,'.r')
                    ax1.plot(dtax2,dtav2,'.g')
                    ax1.set_ylabel('DTA [cm]: ')
                if comp ==1:
                    dtax, dtav= dtafunc(y1,d1,y2,d2,dta)
                    difx, difv= dosedif(y1,d1,y2,d2,1)
                    #resample to input data resolution
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    difx, difv = downsample_to_native(difx, difv, y1)
                    dtainterp=interp.pchip(dtax,dtav)
                    dtaxthresh=np.extract((d1>.01*max(d1)),y1);#find thrshold based on y1
                    dtaxnew=np.extract((dtax>min(dtaxthresh))& (dtax<max(dtaxthresh)),dtax)
                    dtavthresh=dtainterp(dtaxnew);
                    dtafail=np.where(np.abs(dtav)>dta);
                    ddiffail=np.where(np.abs(difv)>dd);
                    totalfailloc=set(dtafail[0])&set(ddiffail[0])
                    totalfail=len(set(dtafail[0])&set(ddiffail[0]))
                    tf=tf+totalfail
                    total=len(dtav)
                    totpf=totpf+total
                    failpf= failpf+totalfail
                    tot=tot+total
                    # add this line to feed the summary:
                    if 'agg' in locals() and isinstance(agg, PassAgg):
                        agg.add(fsl[j], dl[i], total - totalfail, total)

                    #print(dd,dta,fsl[j],totalfail,total,totpf, tot)
                    prtot=(tot-tf)/tot*100
                    prfs=(1-failpf/totpf)*100
                    #prtot=(total-totalfail)/total*100
                    #ax0=plt.subplot(gs[0])
                    ax0.set_xlim(max(min(y1),min(y2))-1,min(max(y1),max(y2))+1)
                    #ax1=fig.add_subplot(gs[1],sharex=ax0)
                    ax1.set_xlim(max(min(y1),min(y2))-1,min(max(y1),max(y2))+1)
                    ax1.plot(dtax,dtav*10,'.g',ms=msize)
                    #ax1.set_ylabel('DTA [mm]')
                    #ax2=fig.add_subplot(gs[2],sharex=ax0)
                    ax2.set_xlim(max(min(y1),min(y2))-1,min(max(y1),max(y2))+1)
                    #ax2.set_xlabel(( 'Points outside of ' + str(dd*100)+'% & '+str(dta*10) +'mm ' +str(tf)+'/'+str(tot) + '  Pass Rate : '+'%.2f' % prtot+'%' ))
                    ax2.plot(difx,difv*100,'.g',ms=msize)#convert to %
                    #ax2.set_ylabel('Dose Difference [%] ')
                    for n in totalfailloc:#plots the failing points
                        ax1.plot(dtax[n],dtav[n]*10,'.r',ms=msize)
                        ax2.plot(difx[n],difv[n]*100,'.r',ms=msize)
                if mppg == 1:
                    # --- metrics ---
                    dtax, dtav = dtafunc(y1, d1, y2, d2, dta)          # DTA (cm)
                    difx, difv = dosedif(y1, d1, y2, d2, 1)            # ΔDose (fraction)
                    #downsample to native resolution
                    dtax, dtav = downsample_to_native(dtax, dtav, y1)
                    difx, difv = downsample_to_native(difx, difv, y1)
                    difv_dmax  = np.asarray(difv) * PDD[i]             # ΔDose × PDD (fraction of dmax)
                
                    # --- common x-limits ---
                    xlo = max(min(y1), min(y2))
                    xhi = min(max(y1), max(y2))
                    ax0.set_xlim(xlo, xhi); ax1.set_xlim(xlo, xhi); ax2.set_xlim(xlo, xhi); ax3.set_xlim(xlo, xhi)
                
                    # --- geometric region masks (SSD-aware) ---
                    # effective field edge at the measurement plane (cm)
                    # get Diagonal Field Size clipping factor from GUI, fall back to global DIAG_Factor if bad input
                    try:
                        diag_factor = float(diag_factor_var.get())
                        
                    except Exception:
                        diag_factor = 1
                    diag     = np.sqrt(2.0* diag_factor) if scl[k] in ('XY', 'YX') else 1.0
                    plane_cm = SSD_CM + dl[i]          # source-to-plane distance at depth
                    mag      = plane_cm / 100      # magnification from SSD plane
                    edge     = (fsl[j] / 2.0) * diag * mag
                    #print(f"diag_factor={diag_factor:.3f}, diag={diag:.3f}, edges=±{edge:.2f} cm")
                    # Get user-entered penumbra region sizes
                    try:
                        PEN_CM = float(pupper_entry.get())/2  # cm, half-width of penumbra
                    except ValueError:
                        PEN_CM = 1.0  # fallback default if user enters bad value
                    
                    try:
                        OVR_CM = float(pulower_entry.get())  # cm, buffer overlap zone
                    except ValueError:
                        OVR_CM = 0.5
        
                
                    # ---- side-aware edges (cm) ----
                    left_edge  = -edge
                    right_edge =  edge
                    
                    # x-coords on each grid
                    x_dif = np.asarray(difx)
                    x_dta = np.asarray(dtax)
                    
                    # ---- core regions (symmetric, true ±PEN_CM around each edge) ----
                    # ΔDose grid
                    pen_left_x   = (x_dif >= left_edge  - PEN_CM) & (x_dif <= left_edge  + PEN_CM)
                    pen_right_x  = (x_dif >= right_edge - PEN_CM) & (x_dif <= right_edge + PEN_CM)
                    pen_core_x   = pen_left_x | pen_right_x
                    
                    in_core_x    = (x_dif >  left_edge  + (PEN_CM + OVR_CM)) & (x_dif <  right_edge - (PEN_CM + OVR_CM))
                    tail_core_x  = (x_dif <= left_edge  - (PEN_CM + OVR_CM)) | (x_dif >= right_edge + (PEN_CM + OVR_CM))
                    
                    # overlap bands (either–or zones)
                    #  - in-field side of each edge
                    near_hi_left_x  = (x_dif >  left_edge + PEN_CM)  & (x_dif <= left_edge + PEN_CM + OVR_CM)
                    near_hi_right_x = (x_dif >= right_edge - (PEN_CM + OVR_CM)) & (x_dif < right_edge - PEN_CM)
                    near_hi_x       = near_hi_left_x | near_hi_right_x
                    
                    #  - tail side of each edge
                    near_lo_left_x  = (x_dif >= left_edge - (PEN_CM + OVR_CM)) & (x_dif < left_edge - PEN_CM)
                    near_lo_right_x = (x_dif >  right_edge + PEN_CM) & (x_dif <= right_edge + PEN_CM + OVR_CM)
                    near_lo_x       = near_lo_left_x | near_lo_right_x
                                                            
                    # ---- DTA grid (penumbra evaluated on dtax) ----
                    pen_left_dta    = (x_dta >= left_edge  - PEN_CM) & (x_dta <= left_edge  + PEN_CM)
                    pen_right_dta   = (x_dta >= right_edge - PEN_CM) & (x_dta <= right_edge + PEN_CM)
                    pen_core_dta    = pen_left_dta | pen_right_dta
                    
                    # ---- minimal debug: region x-coordinates ----
                    def _rng(mask, x):
                        m = np.asarray(mask)
                        if m.size and m.any():
                            xs = x[m]; return f"[{xs.min():.3f}, {xs.max():.3f}] (n={m.sum()})"
                        return "∅ (n=0)"
                    
                    #print(f"Edges: L={left_edge:.3f} cm, R={right_edge:.3f} cm | PEN={PEN_CM:.3f} cm | OVR={OVR_CM:.3f} cm")
                    #print(f"ΔDose pen windows: [{left_edge-PEN_CM:.3f},{left_edge+PEN_CM:.3f}]  "
                    #      f"[{right_edge-PEN_CM:.3f},{right_edge+PEN_CM:.3f}]")
                    #print("ΔDose pen ranges :", "L", _rng(pen_left_x,  x_dif), " | R", _rng(pen_right_x, x_dif))
                    #print("DTA   pen ranges :", "L", _rng(pen_left_dta, x_dta), " | R", _rng(pen_right_dta, x_dta))
                    
                    # keep compat names used later
                    in_field_x   = in_core_x
                    tails_x      = tail_core_x
                    penumbra_dta = pen_core_dta
                
                    # --- pass/fail masks (tolerances unchanged) ---
                    dtav_a       = np.asarray(dtav)            # cm
                    difv_a       = np.asarray(difv)            # fraction
                    pass_dta     = np.abs(dtav_a)    <= dta
                    pass_dd      = np.abs(difv_a)    <= dd
                    pass_ddtail  = np.abs(difv_dmax) <= ddtail
                
                    # DTA evaluated on ΔDose grid for overlap OR checks
                    if len(dtax) > 3:
                        dtainterp_on_difx = interp.pchip(dtax, dtav)(difx)      # cm
                        pass_dta_on_difx  = np.abs(dtainterp_on_difx) <= dta
                    else:
                        pass_dta_on_difx  = np.zeros_like(difx, dtype=bool)
                    # overlap pass/fail using either–or rule
                    overlap_hi_pass = near_hi_x & (pass_dd | pass_dta_on_difx)
                    overlap_hi_fail = near_hi_x & ~(pass_dd | pass_dta_on_difx)
                    overlap_lo_pass = near_lo_x & (pass_ddtail | pass_dta_on_difx)
                    overlap_lo_fail = near_lo_x & ~(pass_ddtail | pass_dta_on_difx)
                    # ===== NEW: MPPG total pass/total on ΔDose grid =====
                    # Evaluate only where MPPG defines criteria (union of all regions on ΔDose grid)
                    domain_mask = in_core_x | tail_core_x | pen_core_x | near_hi_x | near_lo_x
                    
                    # Per-point pass rule per region
                    pass_any = np.zeros_like(domain_mask, dtype=bool)
                    pass_any |= in_core_x   & pass_dd          # in-field → DD
                    pass_any |= tail_core_x & pass_ddtail      # tails    → DD_tails (% of dmax)
                    pass_any |= pen_core_x  & pass_dta_on_difx # penumbra → DTA (on ΔDose grid)
                    pass_any |= overlap_hi_pass                # overlap (in-field side): OR(DD, DTA)
                    pass_any |= overlap_lo_pass                # overlap (tail side):     OR(DD_tails, DTA)
                    
                    t_mp = int(domain_mask.sum())
                    p_mp = int(pass_any.sum())
                    
                    # Record counts for this FS×Depth cell
                    if mppg == 1:
                        agg_mppg.add(fsl[j], dl[i], p_mp, t_mp)
                        # --- Region-wise counts on the ΔDose grid (consistent with overall) ---
                        in_pass   = int((in_core_x   & pass_dd).sum())
                        in_total  = int(in_core_x.sum())
                        
                        # NEW (on ΔDose grid) ✅
                        pen_pass  = int((pen_core_x & pass_dta_on_difx).sum())
                        pen_total = int(pen_core_x.sum())

                        
                        tail_pass  = int((tail_core_x & pass_ddtail).sum())
                        tail_total = int(tail_core_x.sum())
                        
                        ovr_hi_pass = int(overlap_hi_pass.sum())
                        ovr_hi_tot  = int(near_hi_x.sum())
                        
                        ovr_lo_pass = int(overlap_lo_pass.sum())
                        ovr_lo_tot  = int(near_lo_x.sum())
                        
                        agg_mppg.add_region_totals(in_pass, in_total,
                                                   pen_pass, pen_total,
                                                   tail_pass, tail_total,
                                                   ovr_hi_pass, ovr_hi_tot,
                                                   ovr_lo_pass, ovr_lo_tot)



                    # ===== END NEW =====


                    # ========================= PLOTTING =========================
                    # ax1 (DTA)
                    ax1.plot(dtax, dtav_a * 10, '.k', ms=msize)  # base: all values (mm)
                    ax1.plot(dtax[pen_core_dta &  pass_dta], (dtav_a[pen_core_dta &  pass_dta] * 10), '.g', ms=msize)
                    ax1.plot(dtax[pen_core_dta & ~pass_dta], (dtav_a[pen_core_dta & ~pass_dta] * 10), '.r', ms=msize)
                    ax1.set_ylabel('DTA [mm]')
                    
                    # ax2 (ΔDose)
                    ax2.plot(difx, difv_a * 100, '.k', ms=msize)  # base: all values (%)
                    # in-field core
                    ax2.plot(difx[in_core_x &  pass_dd], (difv_a[in_core_x &  pass_dd] * 100), '.', ms=msize, color='g')
                    ax2.plot(difx[in_core_x & ~pass_dd], (difv_a[in_core_x & ~pass_dd] * 100), '.', ms=msize, color='r')
                    # upper overlap (either–or)
                    ax2.plot(difx[overlap_hi_pass], (difv_a[overlap_hi_pass] * 100), '.', ms=msize, color='g')
                    ax2.plot(difx[overlap_hi_fail], (difv_a[overlap_hi_fail] * 100), '.', ms=msize, color='r')
                    ax2.set_ylabel('$\Delta$Dose [%]')
                    
                    
                    # ax3 (ΔDose × PDD)
                    ax3.plot(difx, difv_dmax * 100, '.k', ms=msize)  # base: all values (%)
                    # tail core
                    ax3.plot(difx[tail_core_x &  pass_ddtail], (difv_dmax[tail_core_x &  pass_ddtail] * 100), '.', ms=msize, color='g')
                    ax3.plot(difx[tail_core_x & ~pass_ddtail], (difv_dmax[tail_core_x & ~pass_ddtail] * 100), '.', ms=msize, color='r')
                    # lower overlap (either–or)
                    ax3.plot(difx[overlap_lo_pass], (difv_dmax[overlap_lo_pass] * 100), '.', ms=msize, color='g')
                    ax3.plot(difx[overlap_lo_fail], (difv_dmax[overlap_lo_fail] * 100), '.', ms=msize, color='r')
                    ax3.set_ylabel('$\Delta$Dose [%Dmax]')
                    
                    
                    # keep ΔDose scales matched
                    ax3.set_ylim(ax2.get_ylim())



                 
               


        if comp ==1:
            print('dd: ',dd, 'dta: ',dta, fsl[j],'%.2f' % prfs, failpf, totpf)
        if gam ==1:
            print(fsl[j],'%.2f' % prfs,'%.3f' % np.mean(gvfs),'%.3f' % max(gvfs))    
        
        if comp==1:
            print('dd: ',dd,'dta: ', dta,'All FS: ','%.2f' % prtot)    
            ax1.set_ylabel('DTA [mm]')
            ax2.set_xlabel(
                'Points outside of ' + str(dd*100)+'% & '+str(dta*10) +'mm '
                + str(tf)+'/'+str(tot) + '  Pass Rate : ' + '%.2f' % prtot + '%'
            )
            ax2.set_ylabel('Dose Difference [%] ')
    
            # --- capture composite summary into report_text, still print to terminal ---
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                agg.print_summary()
            report_text = buf.getvalue()
            print(report_text)

    if dif == 1:
        pr_dif = (tot_dif - tf_dif) / tot_dif * 100 if tot_dif > 0 else 0.0
        header = "Field Size [cm] | Fails / Total | Pass Rate [%]"
        sep = "-" * len(header)
        print("\nDose Difference Analysis Results")
        print(f"Criteria: {dd*100:.1f}%")
        print(sep); print(header); print(sep)
        for fs_val, fails, total, pr in dif_results:
            print(f"{fs_val:<15.1f} | {fails:>5}/{total:<7} | {pr:>12.2f}")
        print(sep)
        print(f"{'Overall':<15} | {tf_dif:>5}/{tot_dif:<7} | {pr_dif:>12.2f}")
        ax1.set_xlabel(
            f'Points outside {dd*100:.1f}%  {tf_dif}/{tot_dif}  Pass Rate: {pr_dif:.2f}%'
        )

    if dist == 1:
        dist_prtot = 100 * (tot_dist - tf_dist) / tot_dist if tot_dist > 0 else 0.0
        header = "Field Size [cm] | Fails / Total | Pass Rate [%]"
        sep = "-" * len(header)
        print("\nDTA Analysis Results")
        print(f"Criteria: {dta*10:.1f} mm")
        print(sep); print(header); print(sep)
        for fs_val, fails, total, pr in dist_results:
            print(f"{fs_val:<15.1f} | {fails:>5}/{total:<7} | {pr:>12.2f}")
        print(sep)
        print(f"{'Overall':<15} | {tf_dist:>5}/{tot_dist:<7} | {dist_prtot:>12.2f}")
        ax1.set_xlabel(
            f'Points outside {dta*10:.1f} mm : {tf_dist}/{tot_dist}  Pass Rate: {dist_prtot:.2f}%'
        )

    if mppg == 1:
        _, _, _, _, _, overall = agg_mppg.build_matrix()
        tf_mppg  = agg_mppg.total - agg_mppg.total_pass   # fails
        tot_mppg = agg_mppg.total                          # total
        prtot    = overall                                 # match COMP's naming
    
        print(f'dd: {dd:.4f} ddtail: {ddtail:.4f} dta: {dta:.4f}  All FS overall: {prtot:.2f}')
        ax1.set_ylabel('DTA [mm]')
        ax3.set_xlabel(
            f'Points outside of {dd*100:.1f}% (in-field) | {dta*10:.1f} mm (penumbra) | {ddtail*100:.1f}% of Dmax (out-of-field)  '
            f'{tf_mppg}/{tot_mppg}  Pass Rate : {prtot:.2f}%'
        )

        # --- capture MPPG summary into report_text, still print to terminal ---
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            agg_mppg.print_summary()
        report_text = buf.getvalue()
        print(report_text)



    if gam ==1:
        print('All FS: ','%.2f' % prtot,'%.3f' % np.mean(gvtot),'%.3f' % max(gvtot))

        # --- capture Gamma summary into report_text, still print to terminal ---
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            agg.print_summary()
        report_text = buf.getvalue()
        print(report_text)

        histd=np.histogram(gvtot,bins=np.arange(0,np.max(gvtot),bins),weights=weights);
        ax2.hist(gvtot,bins=np.arange(0,np.max(gvtot),bins),weights=weights)
        ax1.set_ylabel(r'$\Gamma$'+'['+ str(dd*100)+'%/'+ str(dta*10) + 'mm]')
        ax1.set_xlabel( 'Points with ' + r'$\Gamma$' +' > 1 : '+str(np.size((np.where(np.asarray(gvtot)>1))))+'/'+str(len(gvtot)) + '  Pass Rate : '+ str(prtot)+'%' )
        ax2.set_xlabel(r'$\Gamma$'+ '['+ str(dd*100)+'%/'+ str(dta*10) + 'mm]')
        ax2.set_ylabel('Normalized Incidence')

    figManager = plt.get_current_fig_manager()
    try:
        figManager.window.showMaximized()
    except AttributeError:
        figManager.window.state('zoomed')
    fig.tight_layout()
    #fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.16, hspace=0.34)

    if save ==1:
        plt.savefig(str(sn1)+'_correctedvs'+str(sn2)+'_'+str(fsl[j])+'_'+str(int(dd*100))+'_'+str(int(dta*10))+'.png',dpi=400)
    if save ==2:
        plt.savefig('cfOutoffield_'+str(fsl[j])+'_diag',dpi=400)

    plt.pause(0.001)


def save_report():
    global report_text, fig  # Need access to fig for saving
    print("[save_report] called")

    if fig is None:
        print("[save_report] No figure to save. Run an analysis first.")
        return

    if not report_text:
        print("[save_report] No report to save. Run an analysis first.")
        return

    sn1 = sheet1_combo.get()
    sn2 = sheet2_combo.get()
    analysis = analysis_var.get()
    dd = dd_entry.get()
    dta = dta_entry.get()

    print(f"[save_report] sheets: {sn1} vs {sn2}, analysis={analysis}, DD={dd}, DTA={dta}")

    analysis_map = {
        "gam": "Gamma",
        "comp": "Composite",
        "dif": "DoseDiff",
        "dist": "DTA",
        "mppg": "MPPG",
        "none": "NoAnalysis"
    }
    analysis_label = analysis_map.get(analysis, "Analysis")

    default_name = f"ProfileComparison_{sn1}_vs_{sn2}_{analysis_label}_DD{dd}pct_DTA{dta}mm.pdf"

    file_path = filedialog.asksaveasfilename(
        initialfile=default_name,
        defaultextension=".pdf",
        filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
    )

    print(f"[save_report] file_path selected: {file_path!r}")

    if not file_path:
        print("[save_report] User cancelled save dialog.")
        return  # User cancelled

    try:
        if file_path.endswith('.pdf'):
            from reportlab.lib.pagesizes import letter
            from reportlab.pdfgen import canvas
            from reportlab.lib.utils import ImageReader

            c = canvas.Canvas(file_path, pagesize=letter)
            width, height = letter

            # Monospace font for nice alignment
            c.setFont("Courier", 10)

            # Write report_text lines
            lines = report_text.split('\n')
            y = height - 40
            for line in lines:
                c.drawString(40, y, line)
                y -= 14
                if y < 40:
                    c.showPage()
                    c.setFont("Courier", 10)
                    y = height - 40

            # Save figure as temporary PNG
            fig_path = os.path.splitext(file_path)[0] + "_figure.png"
            print(f"[save_report] saving temp figure: {fig_path}")
            fig.savefig(fig_path, dpi=300)

            # Add figure to same PDF
            img = ImageReader(fig_path)
            iw, ih = img.getSize()
            aspect = ih / float(iw)
            desired_width = width - 80
            desired_height = desired_width * aspect

            image_y = y - desired_height - 10
            if image_y < 40:
                c.showPage()
                c.setFont("Courier", 10)
                image_y = height - desired_height - 40

            c.drawImage(img, 40, image_y, width=desired_width, height=desired_height)
            c.save()

            # cleanup temp image
            try:
                os.remove(fig_path)
                print(f"[save_report] removed temp figure: {fig_path}")
            except Exception as e:
                print(f"[save_report] warning: could not remove temp file: {e}")

            print(f"[save_report] Report saved as PDF with figure: {file_path}")

            # Open the PDF automatically
            def open_file(filepath):
                print(f"[save_report] opening file: {filepath}")
                if platform.system() == 'Windows':
                    os.startfile(filepath)
                elif platform.system() == 'Darwin':
                    subprocess.call(['open', filepath])
                else:
                    subprocess.call(['xdg-open', filepath])

            open_file(file_path)
        else:
            print("[save_report] file_path does not end with .pdf, nothing done.")
    except Exception as e:
        print(f"[save_report] ERROR: {e}")


# Initialize main window
root = tk.Tk()
root.title("Profile Compare Tool")

# --- Scrollable shell (Canvas + vertical scrollbar) ---
container = ttk.Frame(root)
container.grid(row=0, column=0, sticky="nsew")
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

canvas = tk.Canvas(container, highlightthickness=0)
vsb = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=vsb.set)

vsb.grid(row=0, column=1, sticky="ns")
canvas.grid(row=0, column=0, sticky="nsew")
container.grid_rowconfigure(0, weight=1)
container.grid_columnconfigure(0, weight=1)

# all your content goes in here
content = ttk.Frame(canvas)
win_id = canvas.create_window((0, 0), window=content, anchor="nw")

def _on_frame_configure(_):
    canvas.configure(scrollregion=canvas.bbox("all"))
content.bind("<Configure>", _on_frame_configure)

def _on_canvas_configure(event):
    # make inner frame match the canvas width
    canvas.itemconfigure(win_id, width=event.width)
canvas.bind("<Configure>", _on_canvas_configure)

# mouse wheel (Windows)
def _on_mousewheel(event):
    canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
canvas.bind_all("<MouseWheel>", _on_mousewheel)

# optional: sensible starting size that respects taskbar
root.geometry("1100x750+80+60")
root.resizable(True, True)

# ------------- BUILD YOUR UI UNDER `content` (only parent changes) -------------

# File selection section
file_frame = ttk.Frame(content, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

# Sheet selection section
sheet_frame = ttk.Frame(content, padding="10")
sheet_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))
ttk.Label(sheet_frame, text="Select Sheet 1:").grid(row=0, column=0, sticky="w")
sheet1_combo = ttk.Combobox(sheet_frame, width=30)
sheet1_combo.grid(row=0, column=1, padx=5)
ttk.Label(sheet_frame, text="Select Sheet 2:").grid(row=1, column=0, sticky="w")
sheet2_combo = ttk.Combobox(sheet_frame, width=30)
sheet2_combo.grid(row=1, column=1, padx=5)

def swap_sheets():
    s1 = sheet1_combo.get()
    s2 = sheet2_combo.get()
    sheet1_combo.set(s2)
    sheet2_combo.set(s1)
    preload_two_sheets()

ttk.Button(sheet_frame, text="Swap \u21c5", command=swap_sheets).grid(row=0, column=2, rowspan=2, padx=5)

# Data Selection Frame
data_selection_frame = ttk.LabelFrame(content, text="Data Selection", padding="10")
data_selection_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))

# Field Size Selection
fsl_label = ttk.Label(data_selection_frame, text="Select Field Sizes")
fsl_label.grid(row=0, column=0, sticky="w")
fsl_listbox = tk.Listbox(data_selection_frame, selectmode="extended", width=20, height=10, exportselection=0)
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
scl_label = ttk.Label(data_selection_frame, text="Select Scan Types")
scl_label.grid(row=0, column=1, sticky="w", padx=20)
scl_listbox = tk.Listbox(data_selection_frame, selectmode="extended", width=20, height=10, exportselection=0)
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

# Depth Selection Listbox
depth_label = ttk.Label(data_selection_frame, text="Select Depths")
depth_label.grid(row=0, column=2, sticky="w", padx=20)
depth_listbox = tk.Listbox(data_selection_frame, selectmode="extended", width=20, height=10, exportselection=0)
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

# Analysis type section with radio buttons
analysis_frame = ttk.LabelFrame(content, text="Analysis Type", padding="10")
analysis_frame.grid(row=3, column=0, sticky="ew")
analysis_var = tk.StringVar(master=root,value="none")
ttk.Radiobutton(analysis_frame, text="Raw Data Only", variable=analysis_var, value="none").grid(row=0, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Plots Only", variable=analysis_var, value="plot").grid(row=1, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Gamma Analysis", variable=analysis_var, value="gam").grid(row=2, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Composite Analysis", variable=analysis_var, value="comp").grid(row=3, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Dose Difference", variable=analysis_var, value="dif").grid(row=4, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Distance to Agreement", variable=analysis_var, value="dist").grid(row=5, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="MPPG5.A Analysis", variable=analysis_var, value="mppg").grid(row=6, column=0, sticky="w")

# Data processing section with dropdowns and checkboxes
processing_frame = ttk.LabelFrame(content, text="Data Processing", padding="10")
processing_frame.grid(row=4, column=0, sticky="ew")
ttk.Label(processing_frame, text="Centering:").grid(row=0, column=0, sticky="w")
cent_combo = ttk.Combobox(processing_frame, values=[0, 1, 2, 3], width=5)
cent_combo.grid(row=0, column=1, padx=5)
cent_combo.set(3)
ttk.Label(processing_frame, text="0: No centering, 1: Center profile 1, 2: Center Profile 2, 3: Center both").grid(row=0, column=2, sticky="w")

ttk.Label(processing_frame, text="Normalization:").grid(row=1, column=0, sticky="w")
norm_combo = ttk.Combobox(processing_frame, values=[0, 1, 2], width=5)
norm_combo.grid(row=1, column=1, padx=5)
norm_combo.set(2)
ttk.Label(processing_frame, text="0: No norm, 1: Normalize to AVG CAX +/- 3mm, 2: Normalize to PDD").grid(row=1, column=2, sticky="w")

ttk.Label(processing_frame, text="Smooth:").grid(row=2, column=0, sticky="w")
smooth_combo = ttk.Combobox(processing_frame, values=[0, 1, 2, 3], width=5)
smooth_combo.grid(row=2, column=1, padx=5)
smooth_combo.set(0)
ttk.Label(processing_frame, text="0: No smoothing, 1: Smooth profile 1, 2: Smooth profile 2, 3: Smooth both").grid(row=2, column=2, sticky="w")

ttk.Label(processing_frame, text=" Dectector Convolution:").grid(row=3, column=0, sticky="w")
conv_combo = ttk.Combobox(processing_frame, values=[0, 1, 2, 3], width=5)
conv_combo.grid(row=3, column=1, padx=5)
conv_combo.set(0)
ttk.Label(processing_frame, text="0: No convolution, 1: Convolve profile 1, 2: Convolve profile 2 3: Convolve both").grid(row=3, column=2, sticky="w")

scale_var = tk.IntVar(value=0)
# Detector diameter input (cm)
det_diam_var = tk.StringVar(master=root)   # bind to current root explicitly
ttk.Label(processing_frame, text="Detector diameter (cm):").grid(row=4, column=0, sticky="w", padx=0)

det_diam_entry = ttk.Entry(processing_frame, textvariable=det_diam_var, width=6)
det_diam_entry.grid(row=4, column=1, sticky="w", padx=5)

# ensure default text appears even in quirky environments
root.after(0, lambda: det_diam_var.set("0.6"))
# Diagonal Field Size clipping factor (user input for DIAG_Factor)
diag_factor_var = tk.StringVar(master=root)
ttk.Label(processing_frame, text="Diagonal Field Size clipping factor:").grid(row=5, column=0, sticky="w", padx=0)
diag_factor_entry = ttk.Entry(processing_frame, textvariable=diag_factor_var, width=6)
diag_factor_entry.grid(row=5, column=1, sticky="w", padx=5)

# default to 0.90
root.after(0, lambda: diag_factor_var.set("0.80"))

ttk.Checkbutton(processing_frame, text="Scale", variable=scale_var).grid(row=6, column=0, sticky="w")
ttk.Label(processing_frame, text="Marker size:").grid(row=7, column=0, sticky="w")
marker_size_entry = ttk.Entry(processing_frame, width=6)
marker_size_entry.grid(row=7, column=1, sticky="w", padx=5)
marker_size_entry.insert(0, "6")

# Criteria frame
criteria_frame = ttk.LabelFrame(content, text="Dose Difference and DTA Criteria Settings: Default:1% and 1mm", padding="10")
criteria_frame.grid(row=5, column=0, sticky="ew")
ttk.Label(criteria_frame, text="Dose Difference (%):").grid(row=0, column=0, sticky="w")
dd_entry = ttk.Entry(criteria_frame, width=10); dd_entry.grid(row=0, column=1, padx=5); dd_entry.insert(0, "2")
ttk.Label(criteria_frame, text="DTA Criteria [mm]:").grid(row=1, column=0, sticky="w")
dta_entry = ttk.Entry(criteria_frame, width=10); dta_entry.grid(row=1, column=1, padx=5); dta_entry.insert(0, "2")
ttk.Label(criteria_frame, text="Dose Difference in tails (%):").grid(row=2, column=0, sticky="w")
ddtail_entry = ttk.Entry(criteria_frame, width=10); ddtail_entry.grid(row=2, column=1, padx=5); ddtail_entry.insert(0, "3")
ttk.Label(criteria_frame, text="Penumbra Width[cm]").grid(row=3, column=0, sticky="w")
pupper_entry = ttk.Entry(criteria_frame, width=10); pupper_entry.grid(row=3, column=1, padx=5); pupper_entry.insert(0, "0.5")
ttk.Label(criteria_frame, text="Penumbra buffer [cm]:").grid(row=4, column=0, sticky="w")
pulower_entry = ttk.Entry(criteria_frame, width=10); pulower_entry.grid(row=4, column=1, padx=5); pulower_entry.insert(0, "1")

# Run comparison button
run_button = ttk.Button(content, text="Run Comparison", command=run_comparison)
run_button.grid(row=6, column=0, pady=10)
#save Repoort Button
save_report_button = ttk.Button(content, text="Save Report (PDF)", command=save_report)
save_report_button.grid(row=7, column=0, pady=5)


def _on_close():
    # unbind global mousewheel so the app doesn't hold input focus
    try:
        canvas.unbind_all("<MouseWheel>")
    except Exception:
        pass
    try:
        root.quit()      # stop Tk's loop (frees Spyder/IPython console)
    finally:
        if root.winfo_exists():
            root.destroy()   # close the window

root.protocol("WM_DELETE_WINDOW", _on_close)

# Start the GUI loop
root.mainloop()
