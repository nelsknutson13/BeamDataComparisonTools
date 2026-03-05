# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 08:37:33 2025

@author: nknutson
"""

#Created by nels knutson. Current code does not report total pass rate for COMP analysis correctly

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import interpolate as interp
from scipy import signal as sig
from matplotlib.patches import Rectangle  # put this with your imports
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from matplotlib.offsetbox import AnchoredText
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
import tkinter as tk
from tkinter import filedialog, ttk
import re
xf = None;sheet_cache={}
# Function to handle file selection and populate sheet names
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
        sheet1_combo.set(sheets[0])
        
        print(f"Sheet 1 selected: {sheet1_combo.get()}")
       

        # (Re)build cache for the two shown sheets, then populate lists
        sheet_cache.clear()
        preload_two_sheets()

        def on_combobox_change(event):
            print(f"Sheet 1 selected: {sheet1_combo.get()}")
            
            preload_two_sheets()
            populate_fsl()
            populate_scl()
            populate_depths()

        sheet1_combo.bind("<<ComboboxSelected>>", on_combobox_change)
        
        populate_fsl()
        populate_scl()
        populate_depths()
        print("Ready: field sizes, scan types, and depths populated.")

    except Exception as e:
        print(f"Error reading sheets: {e}")

def preload_two_sheets():
    """Read only FS, Axis, Depth from the selected sheet ONCE and cache it."""
    global sheet_cache, xf
    usecols = ['FS', 'Axis', 'Depth', 'Pos', 'Dose']
    sn = sheet1_combo.get()
    if not sn or xf is None:
        return
    if sn in sheet_cache:
        return
    try:
        df = pd.read_excel(xf, sheet_name=sn, usecols=usecols, engine="pyarrow")
    except Exception:
        df = pd.read_excel(xf, sheet_name=sn)  # fallback
        # subset if the columns exist
        missing = [c for c in usecols if c not in df.columns]
        if not missing:
            df = df[usecols]
        else:
            print(f"Warning: missing columns in '{sn}': {missing}")
        sheet_cache[sn] = df       
def parse_ssd_energy(path: str, default_ssd=100.0, default_energy="6X"):
    """
    Extract SSD (cm) and Energy (6X, 6FFF, 8FFF, 10X, 10FFF, 15X) from the full path.
    """
    s = str(path)

    # SSD patterns: 'SSD100', 'SSD_100', '100SSD', 'SSD100cm'
    m_ssd = (re.search(r'ssd[_\-\s]?(\d+(?:\.\d+)?)', s, re.IGNORECASE)
             or re.search(r'(\d+(?:\.\d+)?)\s*ssd', s, re.IGNORECASE))
    if m_ssd:
        try:
            ssd_val = float(m_ssd.group(1))
        except:
            ssd_val = float(default_ssd)
    else:
        ssd_val = float(default_ssd)

    # Energy: exactly these labels
    m_energy = re.search(r'\b(6X|6FFF|8FFF|10X|10FFF|15X)\b', s, re.IGNORECASE)
    energy_val = m_energy.group(1).upper() if m_energy else default_energy

    return ssd_val, energy_val
def populate_fsl():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        common_fsl = sorted(set(df1['FS'].dropna().unique()))

        fsl_listbox.delete(0, tk.END)
        for v in common_fsl:
            fsl_listbox.insert(tk.END, v)

        # Auto-select 10.5 if available
        if 10.5 in common_fsl:
            idx = common_fsl.index(10.5)
            fsl_listbox.selection_set(idx)

        print("Field Sizes:", common_fsl)
    except Exception as e:
        print(f"Error populating field sizes: {e}")

# --- TRS-483 Eq.(58): SC_field from the *full* selected profile (assume X≡Y) --

def _Fr_from_full_profile(x_cm, y_rel, S_cm, dr=None):
    """
    Build radial F(r) from a full ±x profile (actual user-selected profile).
    - Normalizes at CAX so F(0)=1.
    - Ensures r-grid extends to the square half-diagonal r_need = (S/2)*sqrt(2).
    - Symmetric average: F(r) = 0.5*(y(+r) + y(-r)).
    """
    x = np.asarray(x_cm, float)
    y = np.asarray(y_rel, float)

    # Normalize to CAX = 1
    i0 = int(np.argmin(np.abs(x)))
    if y[i0] != 0:
        y = y / y[i0]

    # Sort by position
    o = np.argsort(x)
    x = x[o]; y = y[o]

    # Target radial extent (cover the corners of the square)
    r_data_max = float(np.max(np.abs(x)))
    r_need = 0.5 * float(S_cm) * np.sqrt(2.0)
    r_max = max(r_data_max, r_need)

    # Choose dr ~ median dx / 2
    ux = np.unique(x)
    dx = np.diff(ux)
    dx = dx[dx > 0]
    if dr is None:
        dr = float(np.median(dx)/2) if dx.size else 0.01

    r = np.arange(0.0, r_max + 1e-9, dr)

    # Interpolate on both sides and average (assume axial symmetry)
    y_plus  = np.interp(r,  x, y, left=y[i0], right=y[-1])
    y_minus = np.interp(-r, x, y, left=y[0],  right=y[i0])
    F = 0.5 * (y_plus + y_minus)

    # Enforce F(0)=1 and no negatives
    if F[0] != 0:
        F = F / F[0]
    F[F < 0] = 0.0
    return r, F


def _sc_square_eq58_from_Fr(S_cm, r, F, lam=0.18, mu=0.5, n_theta=1440):
    """
    TRS-483 Eq.(58) for a square field of side S_cm using radial F(r):
    SC = (1/2π) ∫_0^{2π} ∫_0^{rmax(θ)} [ λ(1−μ)e^{−λr} + μλ^2 e^{−λr} r ] F(r) dr dθ
    (after the r from the area element cancels the 1/r in the first kernel term).
    """
    r = np.asarray(r, float)
    F = np.asarray(F, float)

    exp = np.exp(-lam * r)
    integrand = (lam*(1.0 - mu) * exp + mu * (lam**2) * exp * r) * F

    # Cumulative ∫_0^R integrand(r) dr via trapezoid on the r grid
    cum = np.concatenate([[0.0], np.cumsum(0.5*(integrand[1:] + integrand[:-1]) * np.diff(r))])

    def I_of_R(R):
        R = np.asarray(R, float)
        Rc = np.clip(R, 0.0, r[-1])
        return np.interp(Rc, r, cum)

    # rmax(θ) for a square with half-side h = S/2
    h = 0.5 * float(S_cm)
    theta = np.linspace(0.0, 2.0*np.pi, n_theta, endpoint=False)
    eps = 1e-12
    rmax = np.minimum(h/np.maximum(np.abs(np.cos(theta)), eps),
                      h/np.maximum(np.abs(np.sin(theta)), eps))

    # (1/2π) ∫ I(rmax(θ)) dθ  ==  average over θ
    return I_of_R(rmax).mean()

def _sc_square_uniform_eq58(S_cm, lam=0.18, mu=0.5, n_theta=2048):
    """
    TRS-483 Eq.(58) for a UNIFORM square (F(r) = 1).
    Uses the closed-form radial primitive:
      I(R) = (1 - e^{-λR}) - μ λ R e^{-λR}
    and averages it over θ via R(θ) for a square.
    """
    S = float(S_cm)
    h = 0.5 * S
    theta = np.linspace(0.0, 2.0*np.pi, n_theta, endpoint=False)
    eps = 1e-12
    R = np.minimum(h/np.maximum(np.abs(np.cos(theta)), eps),
                   h/np.maximum(np.abs(np.sin(theta)), eps))
    # Closed form for ∫ kernel dr with F=1
    exp = np.exp(-lam * R)
    I = (1.0 - exp) - mu * lam * R * exp
    return I.mean()  # = (1/2π)∫ I(R(θ)) dθ


def equivalent_uniform_square_from_profile(S0_cm, r, F, lam=0.18, mu=0.5):
    """
    Solve for S_eq such that SC_uniform(S_eq) = SC_profile(S0).
    Monotone in S, so bisection is robust.
    """
    target = _sc_square_eq58_from_Fr(S0_cm, r, F, lam=lam, mu=mu, n_theta=1440)

    def f(S):
        return _sc_square_uniform_eq58(S, lam=lam, mu=mu) - target

    lo, hi = 1e-3, max(2.0*float(S0_cm), 40.0)  # broad bracket
    flo, fhi = f(lo), f(hi)
    # Expand hi if needed (rare)
    k = 0
    while flo * fhi > 0 and k < 20:
        hi *= 2.0
        fhi = f(hi)
        k += 1

    # Bisection
    for _ in range(64):
        mid = 0.5*(lo + hi)
        fm = f(mid)
        if fm == 0 or (hi - lo) < 1e-6:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo + hi)

def populate_scl():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        raw = df1['Axis'].dropna().unique()
        common_scl = sorted({str(s).strip().upper() for s in raw if str(s).strip().upper() in ('X','Y')})

        scl_listbox.delete(0, tk.END)
        for s in common_scl:
            scl_listbox.insert(tk.END, s)

        # Auto-select X if available
        if "X" in common_scl:
            idx = common_scl.index("X")
            scl_listbox.selection_set(idx)

        print("Scan Types (X/Y only):", common_scl)
    except Exception as e:
        print(f"Error loading scan types: {e}")


def populate_depths():
    try:
        df1 = sheet_cache[sheet1_combo.get()]
        common_depths = sorted(set(df1['Depth'].dropna().unique()))

        depth_listbox.delete(0, tk.END)
        for v in common_depths:
            depth_listbox.insert(tk.END, v)

        # Auto-select 10.0 if available
        if 10.0 in common_depths:
            idx = common_depths.index(10.0)
            depth_listbox.selection_set(idx)

        print("Depths:", common_depths)
    except Exception as e:
        print(f"Error populating depths: {e}")

def _second_deriv_at_zero(x, y, win=0.5):
    """Return f''(0) via local quadratic fit within ±win cm; fallback to 3 nearest pts."""
    x = np.asarray(x); y = np.asarray(y)
    m = (x >= -win) & (x <= win)
    if m.sum() >= 3:
        xf, yf = x[m], y[m]
    else:
        idx = np.argsort(np.abs(x))[:3]
        xf, yf = x[idx], y[idx]
    a, b, c = np.polyfit(xf, yf, 2)  # y ≈ a x^2 + b x + c
    return 2.0 * a


def run_comparison():
    msize = 5  # marker size

    # --- read single selections ---
    sel = fsl_listbox.curselection()
    if not sel:
        print("Select a Field Size.")
        return
    fs = float(fsl_listbox.get(sel[0]))

    sel = scl_listbox.curselection()
    if not sel:
        print("Select a Scan Type.")
        return
    axis_sel = scl_listbox.get(sel[0])

    sel = depth_listbox.curselection()
    if not sel:
        print("Select a Depth.")
        return
    depth_sel = float(depth_listbox.get(sel[0]))

    # --- file/sheet ---
    fn  = file_entry.get()
    sn1 = sheet1_combo.get()
    print(f"Selected File: {fn}")
    print(f"Sheet 1: {sn1}")
    print(f"Plotting FS={fs}, Axis={axis_sel}, Depth={depth_sel}")

    df1 = sheet_cache.get(sn1)
    if df1 is None:
        preload_two_sheets()
        df1 = sheet_cache.get(sn1)
        if df1 is None:
            print(f"No cached data for '{sn1}'.")
            return

    # --- slice selected profile ---
    mask = (df1['FS'] == fs) & (df1['Axis'] == axis_sel) & (df1['Depth'] == depth_sel)
    if mask.sum() == 0:
        print(f"No rows for FS={fs}, Axis={axis_sel}, Depth={depth_sel}")
        return
    x = df1.loc[mask, 'Pos'].to_numpy()
    y = df1.loc[mask, 'Dose'].to_numpy()
    
    # --- normalize to center (CAX) ---
    CAX_WIN_CM = 0.1  # ±1 mm
    cax_mask = np.abs(x) <= CAX_WIN_CM
    if np.any(cax_mask):
        cax_val = float(y[cax_mask].mean())
    else:
        cax_val = float(y[np.argmin(np.abs(x))])  # fallback: nearest to x=0
    if cax_val != 0:
        y = y / cax_val
    else:
        print("Warning: central value is zero; skipping normalization.")
    # --- detector inputs ---
    det_diam_cm = float(det_diam_var.get())
    det_len_cm  = float(det_len_var.get())
    det_orient  = det_orient_var.get() if 'det_orient_var' in globals() else "Assume Symmetry"
    
    
    
    # --- TRS-483 lateral magnification → map to 100 cm plane ---
    # M = (SSD + depth) / 100
    M = (float(SSD_CM) + float(depth_sel)) / 100.0
    print(M)
    # Use positions at 100 cm plane for Eq.(58)
    #x_100 = x / M
    S_cm = float(fs)  # square side = selected field size
    r_rad, F_rad = _Fr_from_full_profile(x, y, S_cm=S_cm)  # uses full ±x profile
    SC_field = _sc_square_eq58_from_Fr(S_cm, r_rad, F_rad, lam=0.18, mu=0.5, n_theta=1440)
    S_eq_uniform = (1/M)*equivalent_uniform_square_from_profile(S_cm, r_rad, F_rad, lam=0.18, mu=0.5)
    print(f"Equivalent UNIFORM square (F=1): S_eq = {S_eq_uniform:.3f} cm")

    print(f"TRS-483 SC_field (Eq.58, profile-weighted, S={S_cm:.2f} cm): {SC_field:.6f}")


    

    # --- curvature along selected axis (needed for VAE) ---
    d2_sel = _second_deriv_at_zero(x, y, win=0.5)  # f''(0) along the plotted axis
    
    # helper to get curvature from the *other* axis if available
    def _curv_for_axis(axis_letter):
        if axis_letter == axis_sel:
            return d2_sel
        mask_other = (df1['FS'] == fs) & (df1['Axis'] == axis_letter) & (df1['Depth'] == depth_sel)
        if mask_other.sum() >= 3:
            xo = df1.loc[mask_other, 'Pos'].to_numpy()
            yo = df1.loc[mask_other, 'Dose'].to_numpy()
            # normalize same way (±1 mm window)
            CAX_WIN_CM = 0.1
            cmask = np.abs(xo) <= CAX_WIN_CM
            cax = float(yo[cmask].mean()) if np.any(cmask) else float(np.interp(0.0, xo, yo))
            if cax != 0:
                yo = yo / cax
            return _second_deriv_at_zero(xo, yo, win=0.5)
        return np.nan

    # map orientation to axes
    if det_orient == "Short Axis X, Long Axis Y":
        short_dir, long_dir = 'X', 'Y'
    elif det_orient == "Short Axis Y, Long Axis X":
        short_dir, long_dir = 'Y', 'X'
    else:  # "Assume Symmetry"
        short_dir = long_dir = None
    
    # curvatures for VAE (use symmetry fallback if needed)
    if short_dir is None:   # symmetry assumption
        d2_short = d2_sel
        d2_long  = d2_sel
    else:
        d2_short = _curv_for_axis(short_dir)
        d2_long  = _curv_for_axis(long_dir)
        if np.isnan(d2_short): d2_short = d2_sel
        if np.isnan(d2_long):  d2_long  = d2_sel
    
    # --- VAE & Prp for a cylindrical chamber about CAX ---
    # Short axis (circular cross-section, radius r=b/2):  ⟨x²⟩ = r²/4 → coeff = b²/32
    # Long axis (average along a line of length a):       ⟨x²⟩ = a²/12 → coeff = a²/24
    eps_short = d2_short * (det_diam_cm**2) / 32.0   # fraction
    eps_long  = d2_long  * (det_len_cm**2)  / 24.0   # fraction
    eps_total = eps_short + eps_long                  # fraction
    
    # optional: report in %
    vae_short_pct = eps_short * 100.0
    vae_long_pct  = eps_long  * 100.0
    vae_total_pct = eps_total * 100.0
    
    # correction factor (true/measured at CAX)
    Prp = 1.0 / (1.0 + eps_total)
    
    print(f"VAE short-axis (cyl)  ±{det_diam_cm/2:.2f} cm: {vae_short_pct:.3f}%")
    print(f"VAE long-axis  (line) ±{det_len_cm/2:.2f} cm:  {vae_long_pct:.3f}%")
    print(f"VAE combined: {vae_total_pct:.3f}%")
    print(f"Prp (point correction factor): {Prp:.5f}")


    # --- plot ---
    #gs = gridspec.GridSpec(1, 1)
    fig, ax0 = plt.subplots(1, 1, figsize=(10, 5))
    plt.rcParams.update({'font.size': 16})

   #ax0 = plt.subplot(gs[0])
    ax0.plot(x, y, '.k', ms=msize)                 # selected profile
    ax0.legend()
    ax0.set_ylabel('Off Axis Ratio')
    ax0.set_xlabel('Off Axis Position [cm]')
    ax0.set_title(f"FS {fs} cm | Axis {axis_sel} | Depth {depth_sel} cm")
    ax0.set_xlim(x.min(), x.max())
    # --- visualize detector averaging window(s) about CAX as rectangles ---
    ymin, ymax = ax0.get_ylim()
    # center value after your CAX normalization (≈1.0); fallback to nearest-to-zero point
    y0 = 1.0 if np.isfinite(y).all() else float(y[np.argmin(np.abs(x))])
    
    # rectangle height: small band around CAX (e.g., 6% of y-range)
    h = 0.06 * (ymax - ymin)
    
    def draw_rect(half_width, label_text, face_alpha=0.18, z=0):
        rect = Rectangle(
            (-half_width, y0 - h/2),   # (x0, y0) lower-left
            2 * half_width,            # width
            h,                         # height
            linewidth=1,
            fill=True,
            alpha=face_alpha,
            zorder=z
        )
        ax0.add_patch(rect)
        # label near the rectangle (bottom-center in axes coords)
        ax0.text(0.5, 0.02, label_text, transform=ax0.transAxes,
                 ha='center', va='bottom')
    


    

    # --- annotation block ---
    det_orient_plot = det_orient_var.get()
    info_lines = []
    
    if det_orient_plot == "Short Axis X, Long Axis Y":
        half_short = det_diam_cm / 2.0
        draw_rect(half_short, f"Short-axis (X) avg window ±{half_short:.2f} cm")
        info_lines.append(f"Long-axis (Y): {det_len_cm:.2f} cm ⟂")
        info_lines.append(f"Short-axis (X): ±{det_diam_cm/2:.2f} cm")
    
    elif det_orient_plot == "Short Axis Y, Long Axis X":
        half_long = det_len_cm / 2.0
        draw_rect(half_long, f"Long-axis (X) avg window ±{half_long:.2f} cm")
        info_lines.append(f"Long-axis (X): ±{det_len_cm/2:.2f} cm")
        info_lines.append(f"Short-axis (Y): {det_diam_cm:.2f} cm ⟂")
    
    else:  # "Assume ..."
        h_short = det_diam_cm / 2.0
        h_long  = det_len_cm  / 2.0
        big, small = (h_long, h_short) if h_long >= h_short else (h_short, h_long)
        draw_rect(big,   "", face_alpha=0.10, z=0)
        draw_rect(small, "", face_alpha=0.20, z=1)
        info_lines.append("Assuming X and Y profiles are same")
        info_lines.append(f"Long Axis:  ±{big:.2f} cm")
        info_lines.append(f"Short Axis: ±{small:.2f} cm")
        info_lines.append(f"TRS-483 SC_field: {SC_field:.6f}")
        info_lines.append(f"S_eq (uniform F=1): {S_eq_uniform:.2f} cm")


    # --- concise contribution line ---
    if det_orient_plot == "Short Axis X, Long Axis Y":
        short_lab, long_lab = "Short (X)", "Long (Y)"
    elif det_orient_plot == "Short Axis Y, Long Axis X":
        short_lab, long_lab = "Short (Y)", "Long (X)"
    else:
        short_lab, long_lab = "Short", "Long"
    
      # --- concise contribution lines ---
    info_lines.append(
        f"VAE Short: {vae_short_pct:+.2f}%, Long: {vae_long_pct:+.2f}%"
    )
    info_lines.append(
        f"VAE Total: {vae_total_pct:+.2f}%"
    )
    info_lines.append(
        f"Prp = {Prp:.5f}"
)
    
   
    
    # AnchoredText box in top right
    info_text = "\n".join(info_lines)
    anchored = AnchoredText(info_text, loc="upper right",
                            prop={'size': 12}, frameon=True, borderpad=0.5)
    anchored.patch.set_alpha(0.7)
    ax0.add_artist(anchored)


    # window niceties
    try:
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
    except Exception:
        pass
    fig.tight_layout()
    fig.subplots_adjust(hspace=.3)


            
# Initialize main window
root = tk.Tk()
root.title("Profile Compare Tool")

# File selection section
file_frame = ttk.Frame(root, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

# Sheet selection section
sheet_frame = ttk.Frame(root, padding="10")
sheet_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))  # Add some padding to separate from the file frame
ttk.Label(sheet_frame, text="Select Sheet 1:").grid(row=0, column=0, sticky="w")
sheet1_combo = ttk.Combobox(sheet_frame, width=30)
sheet1_combo.grid(row=0, column=1, padx=5)

# Data Selection Frame
data_selection_frame = ttk.LabelFrame(root, text="Data Selection", padding="10")
data_selection_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))  # Adjust row number to fit below file and sheet sections

# Field Size Selection
fsl_label = ttk.Label(data_selection_frame, text="Select Field Sizes")
fsl_label.grid(row=0, column=0, sticky="w")
fsl_listbox = tk.Listbox(data_selection_frame, selectmode="browse", width=20, height=10, exportselection=0)
fsl_listbox.grid(row=1, column=0, padx=5, pady=5)

# Scan Types Selection
scl_label = ttk.Label(data_selection_frame, text="Select Scan Types")
scl_label.grid(row=0, column=1, sticky="w", padx=20)
scl_listbox = tk.Listbox(data_selection_frame, selectmode="browse", width=20, height=10, exportselection=0)
scl_listbox.grid(row=1, column=1, padx=5, pady=5)

# Depth Selection Listbox
depth_label = ttk.Label(data_selection_frame, text="Select Depths")
depth_label.grid(row=0, column=2, sticky="w", padx=20)  # Place this next to the other selection lists
depth_listbox = tk.Listbox(data_selection_frame, selectmode="browse", width=20, height=10, exportselection=0)
depth_listbox.grid(row=1, column=2, padx=5, pady=5)

# Data processing section with two entries
processing_frame = ttk.LabelFrame(root, text="Data Processing", padding="10")
processing_frame.grid(row=4, column=0, sticky="ew")

# Detector length input (cm)
det_len_var = tk.StringVar(master=root,value="2.3")
tk.Label(processing_frame, text="Detector length (cm):").grid(row=0, column=0, sticky="w")
tk.Entry(processing_frame, textvariable=det_len_var, width=6).grid(row=0, column=1, sticky="w", padx=5)

# Detector diameter input (cm)
det_diam_var = tk.StringVar(master=root,value="0.61")
tk.Label(processing_frame, text="Detector diameter (cm):").grid(row=1, column=0, sticky="w")
tk.Entry(processing_frame, textvariable=det_diam_var, width=6).grid(row=1, column=1, sticky="w", padx=5)
# Detector orientation
det_orient_var = tk.StringVar(value="Assume X and Y profile is same")
ttk.Label(processing_frame, text="Detector orientation:").grid(row=2, column=0, sticky="w")
det_orient_combo = ttk.Combobox(
    processing_frame,
    textvariable=det_orient_var,
    state="readonly",
    width=28,
    values=[
        "Short Axis X, Long Axis Y",
        "Short Axis Y, Long Axis X",
        "Assume X and Y profile is same",
    ],
)
det_orient_combo.grid(row=2, column=1, sticky="w", padx=5)

# --- NEW: Explicitly set values after widget creation ---
det_len_var.set("2.3")
det_diam_var.set("0.61")
det_orient_combo.set("Assume X and Y profile is same")
root.update_idletasks()
# Run comparison button
run_button = ttk.Button(root, text="Run Comparison", command=run_comparison)
run_button.grid(row=6, column=0, pady=10)


def _on_close():
    root.quit()     # stop Tk’s loop so console is released
    root.destroy()  # close the GUI window

root.protocol("WM_DELETE_WINDOW", _on_close)




# Start the GUI loop
root.mainloop()