import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import interpolate as interp
from scipy import signal as sig
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import BSpline, PPoly
from gamma import gamma as g
from center import center as c
from comp import dosedif as dosedif
from comp import dta as dtafunc
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import filedialog, ttk
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import sys, os, signal

import subprocess
import platform
import enum
from typing import NamedTuple, Literal, Union
from scipy.ndimage import gaussian_filter1d
# Function to handle file selection and populate sheet names




# ---------------- Ion chamber models ----------------

class IonChamberData(NamedTuple):
    radius_cavity_mm: float  # mm
    wall_thickness_mm: float # mm

class IonChambers(str, enum.Enum):
    PTW_31021 = "PTW31021"
    IBA_CC13  = "IBA CC13"
    PTW_31010 ='PTW31010'
    def _parameters(self) -> IonChamberData:
        if self is IonChambers.PTW_31021:
            # PTW 31021: rcav=2.4 mm, wall=0.57+0.09 mm
            return IonChamberData(2.4, 0.57 + 0.09)
        if self is IonChambers.PTW_31010:
            # PTW 31010: rcav=2.75 mm, wall=0.55+0.15 mm
            return IonChamberData(2.75, 0.55 + 0.15)
        if self is IonChambers.IBA_CC13:
            # IBA CC13: rcav=3.0 mm, wall=0.4 mm
            return IonChamberData(3.0, 0.4)
        # fallback (should not hit)
        return IonChamberData(2.4, 0.66)

    @property
    def outer_radius(self) -> float:
        p = self._parameters()
        return p.radius_cavity_mm + p.wall_thickness_mm  # mm

    @property
    def simple_photon_effective_point_of_measurement(self) -> float:
        return 0.6 * self._parameters().radius_cavity_mm  # mm

    @property
    def simple_electron_effective_point_of_measurement(self) -> float:
        return 0.5 * self._parameters().radius_cavity_mm  # mm


class DistanceUnit(str, enum.Enum):
    cm = "cm"
    mm = "mm"

class Modality(str, enum.Enum):
    PHOTON  = "PHOTON"
    ELECTRON = "ELECTRON"

_Unit = Union[Literal["cm","mm"], DistanceUnit]

# ---------------- PDD curve & shifts ----------------

class PDDCurve:
    """
    Cubic smoothing-spline style PDD with EPOM/center shifts.
    """
    def __init__(self, z, dose, z_distance_unit: _Unit = "cm", lam: float = 0.0):
        z = np.asarray(z, dtype=float)
        dose = np.asarray(dose, dtype=float)
        if z.shape != dose.shape:
            raise ValueError("z and dose must have the same shape")

        # sort by depth
        idx = np.argsort(z)
        self.z = z[idx]
        self.dose = dose[idx]

        # store smoothness & map to UnivariateSpline's 's'
        self.lam = float(lam)
        N = max(1, self.z.size)
        s = max(0.0, self.lam * N * float(np.var(self.dose)))  # scale-invariant heuristic

        # choose degree: cubic if we have enough points, else degrade gracefully
        k_used = min(3, max(1, N - 1))  # need ≥4 points for cubic
        self._k_used = k_used
        self.curve = UnivariateSpline(self.z, self.dose, s=s, k=k_used)

        # units
        if isinstance(z_distance_unit, DistanceUnit):
            z_distance_unit = z_distance_unit.value
        if z_distance_unit not in ("cm", "mm"):
            raise ValueError("z_distance_unit must be 'cm' or 'mm'")
        self.z_distance_unit = z_distance_unit

    # ---- internal: convert FITPACK spline -> PPoly (works for roots of any order) ----
    def _ppoly(self) -> PPoly:
        t, c, k = self.curve._eval_args  # FITPACK representation
        return PPoly.from_spline(BSpline(t, c, k))

    # ---- derivatives/roots using PPoly ----
    def derivative(self, order: int = 1):
        # Return a callable derivative like before, backed by PPoly
        return self._ppoly().derivative(order)

    def derivative_roots(self, order: int = 1) -> np.ndarray:
        return self._ppoly().derivative(order).roots()

    # ---- unit conversion ----
    def _mm_to_z_units(self, mm: float) -> float:
        return mm / (10.0 if self.z_distance_unit == "cm" else 1.0)

    # ---- surface locator (original method with safe fallback) ----
    def getIonChamberAtWaterSurfaceLocation(self) -> float:
        """
        Original method:
          - Find roots of the 2nd derivative; choose the root where d¹ is maximal.
          - If roots missing, fall back to dense-grid argmax near the entrance.
        Returns depth in the same units as self.z.
        """
        try:
            roots2 = self.derivative_roots(2)
            if roots2.size:
                d1_vals = self.derivative(1)(roots2)
                print("[Auto-shift] using cubic roots method (via PPoly)")
                return float(roots2[int(np.argmax(d1_vals))])
        except Exception:
            pass  # fall through to grid fallback
        print("[Auto-shift] using dense-grid fallback")

        # Fallback: dense-grid argmax over first ~2 cm (or 20 mm)
        zmin, zmax = float(self.z.min()), float(self.z.max())
        win = 2.0 if self.z_distance_unit == "cm" else 20.0
        zhi = min(zmax, zmin + win)
        n = max(1000, 10 * len(self.z))
        zfine = np.linspace(zmin, zhi, n)
        d1 = self.derivative(1)(zfine)
        d1 = np.nan_to_num(d1, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
        return float(zfine[int(np.argmax(d1))])

    # ---- shifts ----
    def getShiftToIonChamberCenter(self, ion_chamber: IonChambers) -> float:
        """Shift (in z units) so the chamber center is at the water surface (z=0)."""
        outer = self._mm_to_z_units(ion_chamber.outer_radius)
        return outer - self.getIonChamberAtWaterSurfaceLocation()

    def getShiftToIonChamberEPOM(self, ion_chamber: IonChambers, modality: Modality) -> float:
        """Shift (in z units) so the EPOM is at the water surface (z=0)."""
        if modality == Modality.PHOTON:
            epom_mm = ion_chamber.simple_photon_effective_point_of_measurement
        else:
            epom_mm = ion_chamber.simple_electron_effective_point_of_measurement
        offset = self._mm_to_z_units(ion_chamber.outer_radius - epom_mm)  # mm -> z units
        return offset - self.getIonChamberAtWaterSurfaceLocation()

    # ---- shifted copies ----
    def _refit(self):
        """Refit spline after in-place z shift, keeping same smoothness intent."""
        N = max(1, self.z.size)
        s = max(0.0, self.lam * N * float(np.var(self.dose)))
        self.curve = UnivariateSpline(self.z, self.dose, s=s, k=self._k_used)

    def shiftPDDToIonChamberCenter(self, ion_chamber: IonChambers, in_place: bool = False):
        sft = self.getShiftToIonChamberCenter(ion_chamber)
        if in_place:
            self.z = self.z + sft
            self._refit()
            return self
        return PDDCurve(self.z + sft, self.dose, self.z_distance_unit, self.lam)

    def shiftPDDToIonChamberEPOM(self, ion_chamber: IonChambers, modality: Modality, in_place: bool = False):
        sft = self.getShiftToIonChamberEPOM(ion_chamber, modality)
        if in_place:
            self.z = self.z + sft
            self._refit()
            return self
        return PDDCurve(self.z + sft, self.dose, self.z_distance_unit, self.lam)


# ---------------- Optional: tiny helpers for your GUI ----------------

def map_detector_name(name: str) -> IonChambers:
    key = (name or "").strip().lower().replace("-", "").replace("_", "").replace(" ", "")
    if key in ("ptw31021", "ptw021"): return IonChambers.PTW_31021
    if key in ("ibacc13", "cc13", "ibacc-13", "ibacc 13", "ibacc", "ibacc-13"): return IonChambers.IBA_CC13
    if key in ("ptw31010", "ptw010"): return IonChambers.PTW_31010
    return IonChambers.PTW_31021

def compute_epom_shifts(z1, d1, det1, z2, d2, det2, modality="PHOTON", unit="cm", lam=0.0):
    """
    Convenience function: returns (s1_cm, s2_cm) using EPOM for each curve.
    """
    ic1 = map_detector_name(det1); ic2 = map_detector_name(det2)
    mod = Modality(modality) if isinstance(modality, str) else modality
    p1 = PDDCurve(z1, d1, z_distance_unit=unit, lam=lam)
    p2 = PDDCurve(z2, d2, z_distance_unit=unit, lam=lam)
    s1 = p1.getShiftToIonChamberEPOM(ic1, mod)
    s2 = p2.getShiftToIonChamberEPOM(ic2, mod)
    return float(s1), float(s2)
def downsample_to_native(x_interp, vals_interp, x_native):
    """
    Downsample interpolated data (x_interp, vals_interp)
    to roughly match the spacing of x_native.

    Enhancement:
      - If x_native is a list/tuple of arrays (e.g., [y1, y2]),
        it will match the COARSER (larger) of their median step sizes.
    """
    x_interp = np.asarray(x_interp, dtype=float)
    vals_interp = np.asarray(vals_interp)
    
    # Allow x_native to be a single array OR a list/tuple of arrays
    if isinstance(x_native, (list, tuple)):
        natives = [np.asarray(x, dtype=float) for x in x_native if x is not None]
    else:
        natives = [np.asarray(x_native, dtype=float)]

    if x_interp.size < 2 or len(natives) == 0:
        return x_interp, vals_interp  # nothing to do

    def _median_step(x):
        x = x[np.isfinite(x)]
        if x.size < 2:
            return 0.0
        xu = np.unique(x)
        xu.sort()
        diffs = np.diff(xu)
        diffs = diffs[diffs > 0]
        return float(np.median(diffs)) if diffs.size else 0.0

    # COARSER native step = max(step sizes), ignoring zeros
    steps = [s for s in (_median_step(x) for x in natives) if s > 0]
    if not steps:
        return x_interp, vals_interp

    native_step = max(steps)

    # Interp step from x_interp
    xi = np.unique(x_interp[np.isfinite(x_interp)])
    xi.sort()
    if xi.size < 2:
        return x_interp, vals_interp
    interp_step = float(np.median(np.diff(xi)))
    if interp_step <= 0:
        return x_interp, vals_interp

    factor = max(1, int(round(native_step / interp_step)))
    return x_interp[::factor], vals_interp[::factor]

def _surface_resolution(y, unit="cm", window=None, max_pairs=5) -> float:
    """
    Estimate scan resolution near the water surface using the first few
    spacings close to the minimum depth in y.
      - window: span (in same units as y) over which to sample steps.
                Defaults to 0.5 cm or 5 mm depending on unit.
      - Uses median of the first `max_pairs` positive diffs within window.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if y.size < 2:
        return 0.0

    # default window by unit
    if window is None:
        window = 0.5 if (unit or "").lower() == "cm" else 5.0  # cm or mm

    yu = np.unique(y)
    yu.sort()
    diffs = np.diff(yu)

    # consider only diffs starting near the smallest depth
    start = yu[0]
    within = (yu[:-1] - start) <= float(window)
    cand = diffs[within]
    if cand.size == 0:
        cand = diffs[:max_pairs]
    cand = cand[cand > 0]
    return float(np.median(cand)) if cand.size else 0.0


def select_file():
    # Open file dialog to select an Excel file
    file_path = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
    file_entry.delete(0, tk.END)
    file_entry.insert(0, file_path)
    
    # Read the sheet names from the selected Excel file
    if file_path:
        try:
            sheets = pd.ExcelFile(file_path).sheet_names
            # Populate the dropdown menus with the sheet names
            sheet1_combo['values'] = sheets
            sheet2_combo['values'] = sheets
            sheet1_combo.set(sheets[0])  # Defaults to first sheet
            sheet2_combo.set(sheets[1])  # Defaults to second sheet
            
            # Function to handle combobox selection changes
            def on_combobox_change(event):
                # Actions to take when a combobox selection changes
                print(f"Sheet 1 selected: {sheet1_combo.get()}")
                print(f"Sheet 2 selected: {sheet2_combo.get()}")
                # Repopulate field sizes after a new sheet is selected
                populate_fsl()
                populate_scl()

            # Bind the change event to the function
            sheet1_combo.bind("<<ComboboxSelected>>", on_combobox_change)
            sheet2_combo.bind("<<ComboboxSelected>>", on_combobox_change)
            
            # Populate field sizes after file and sheet selections
            populate_fsl()
            populate_scl()
        
        except Exception as e:
            print(f"Error reading Excel file: {e}")



def populate_scl():
    try:
        # Load data from the selected sheets
        df1 = pd.read_excel(file_entry.get(), sheet_name=sheet1_combo.get())
        df2 = pd.read_excel(file_entry.get(), sheet_name=sheet2_combo.get())

        # Extract unique scan types from the "Axis" column of both dataframes
        scl1 = set(df1['Axis'].unique())
        scl2 = set(df2['Axis'].unique())

        # Find common scan types available in both datasets
        common_scl = sorted(list(scl1.intersection(scl2)))

        # Clear the listbox before inserting new values
        scl_listbox.delete(0, tk.END)

        # Populate the listbox with the common scan types
        for scan in common_scl:
            scl_listbox.insert(tk.END, scan)
            scl_listbox.selection_set(0, tk.END)  # Select all by default

        # Debugging: Print to ensure values are correct
        print("Common Scan Types:", common_scl)

    except Exception as e:
        print(f"Error loading scan types: {e}")

# Add this call to populate_scl() at the end of your file and sheet selection logic

# Function to populate the field sizes (fsl) based on selected sheets
def populate_fsl():
   try:
    
   # Example extraction from the sheets (replace with actual data extraction logic)
    # Assuming 'df1' and 'df2' are dataframes loaded from sn1 and sn2 respectively
        df1 = pd.read_excel(file_entry.get(), sheet_name=sheet1_combo.get())
        df2 = pd.read_excel(file_entry.get(), sheet_name=sheet2_combo.get())
    
    # Extract unique field sizes from both dataframes
        fsl1 = set(df1['FS'].unique())
        fsl2 = set(df2['FS'].unique())
    
    # Find common field sizes available in both datasets
        common_fsl = sorted(list(fsl1.intersection(fsl2)))
    
        # Print the common field sizes to verify they are correct
        print("Common Field Sizes:", common_fsl)    
    
    # Clear the listbox before inserting new values
        fsl_listbox.delete(0, tk.END)
    
        # Populate the listbox with the common field sizes
        for value in common_fsl:
            fsl_listbox.insert(tk.END, value)
            #fsl_listbox.selection_set(0, tk.END)  # Select all by default
   except Exception as e:
        print(f"Error populating field sizes: {e}")
        
def save_current_figure():
    try:
        from tkinter import filedialog
        # Build default filename based on selected sheets, analysis type, and criteria
        sn1 = sheet1_combo.get()
        sn2 = sheet2_combo.get()
        analysis = analysis_var.get()
        
        # Dose diff and DTA values from entries
        dd = dd_entry.get()
        dta = dta_entry.get()

        # Friendly label for analysis type
        analysis_map = {
            "gam": "Gamma",
            "comp": "Composite",
            "dif": "DoseDiff",
            "dist": "DTA",
            "none": "NoAnalysis"
        }
        analysis_label = analysis_map.get(analysis, "Analysis")

        default_name = f"PDDComparison_{sn1}_vs_{sn2}_{analysis_label}_DD{dd}pct_DTA{dta}mm.png"

        file_path = filedialog.asksaveasfilename(
            initialfile=default_name,
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        if file_path:
            fig = plt.gcf()   # get current active figure
            fig.savefig(file_path, dpi=400)
            print(f"Figure saved to {file_path}")
    except Exception as e:
        print(f"Error saving figure: {e}")

def save_report():
    global report_text, fig  # Need access to fig for saving
    if fig is None:
        print("No figure to save. Run an analysis first.")
        return

    if not report_text:
        print("No report to save. Run an analysis first.")
        return

    sn1 = sheet1_combo.get()
    sn2 = sheet2_combo.get()
    analysis = analysis_var.get()
    dd = dd_entry.get()
    dta = dta_entry.get()

    analysis_map = {
        "gam": "Gamma",
        "comp": "Composite",
        "dif": "DoseDiff",
        "dist": "DTA",
        "none": "NoAnalysis"
    }
    analysis_label = analysis_map.get(analysis, "Analysis")

    default_name = f"PDDComparison_{sn1}_vs_{sn2}_{analysis_label}_DD{dd}pct_DTA{dta}mm.pdf"

    file_path = filedialog.asksaveasfilename(
        initialfile=default_name,
        defaultextension=".pdf",
        filetypes=[("PDF files", "*.pdf"), ("All files", "*.txt"), ("All files", "*.*")]
    )

    if not file_path:
        return  # User cancelled

    if file_path.endswith('.pdf'):
        c = canvas.Canvas(file_path, pagesize=letter)
        width, height = letter

        # Set monospace font for alignment
        c.setFont("Courier", 10)

        # Write text lines
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
        fig_path = os.path.splitext(file_path)[0] + "_figuretest.png"
        fig.savefig(fig_path, dpi=300)
        print(f"Temp figure saved at: {fig_path}") 
        # Add a new page for the figure
    # Load image and get size to maintain aspect ratio
        img = ImageReader(fig_path)
        iw, ih = img.getSize()
        aspect = ih / float(iw)
        desired_width = width - 80
        desired_height = desired_width * aspect

        # Add a new page and draw the image
        #c.showPage()
        image_y = y - desired_height - 10  # y is last text Y position
        if image_y < 40:
            c.showPage()
            c.setFont("Courier", 10)
            image_y = height - desired_height - 40
        c.drawImage(img, 40, image_y, width=desired_width, height=desired_height)


        c.save()

        #Optionally delete the temporary image file if you want:
        os.remove(fig_path)

        print(f"Report saved as PDF with figure: {file_path}")
        # Open the PDF automatically
        def open_file(filepath):
            if platform.system() == 'Windows':
                os.startfile(filepath)
            elif platform.system() == 'Darwin':  # macOS
                subprocess.call(['open', filepath])
            else:  # Linux and others
                subprocess.call(['xdg-open', filepath])

        open_file(file_path)


def apply_detector_convolution(profile_y, profile_dose, fwhm_cm):
    if fwhm_cm <= 0:
        return profile_dose  # No convolution needed

    # Convert FWHM (cm) to sigma (cm)
    sigma_cm = fwhm_cm / 2.355

    # Estimate step size in cm (assumes uniform spacing)
    step_cm = np.mean(np.diff(profile_y))

    # Convert sigma from cm to number of data points
    sigma_idx = sigma_cm / step_cm

    # Apply Gaussian convolution/smoothing
    dose_smoothed = gaussian_filter1d(profile_dose.values, sigma=sigma_idx)

    return dose_smoothed

def _to_ionchamber(name: str) -> IonChambers:
    key = (name or "").strip().lower().replace("-", "").replace("_", "").replace(" ", "")
    if key in ("ptw31021",): return IonChambers.PTW_31021
    if key in ("ptw31010",): return IonChambers.PTW_31010
    if key in ("ibacc13","cc13","ibacc-13","iba cc13"): return IonChambers.IBA_CC13
    return IonChambers.PTW_31021

def autofill_depth_shifts_epom_ui(y1, d1, y2, d2, unit="cm"):
    try:
        ic1 = _to_ionchamber(det1_combo.get())
        ic2 = _to_ionchamber(det2_combo.get())
        mod = Modality(modality_combo.get()) if modality_combo.get() in ("PHOTON","ELECTRON") else Modality.PHOTON

        p1 = PDDCurve(np.asarray(y1, float), np.asarray(d1, float), z_distance_unit=unit, lam=0.0)
        p2 = PDDCurve(np.asarray(y2, float), np.asarray(d2, float), z_distance_unit=unit, lam=0.0)

        # --- Raw EPOM shifts ---
        s1_raw = p1.getShiftToIonChamberEPOM(ic1, mod)
        s2_raw = p2.getShiftToIonChamberEPOM(ic2, mod)

        # --- Local resolution near surface ---
        res1 = _surface_resolution(y1, unit=unit)
        res2 = _surface_resolution(y2, unit=unit)

        # --- Apply half-step correction (shift shallower by 0.5 * res) ---
        s1 = s1_raw + 0.5 * res1 if res1 > 0 else s1_raw
        s2 = s2_raw + 0.5 * res2 if res2 > 0 else s2_raw

        # --- Update GUI entries ---
        depth_shift1_entry.delete(0, tk.END)
        depth_shift1_entry.insert(0, f"{s1:.3f}")
        depth_shift2_entry.delete(0, tk.END)
        depth_shift2_entry.insert(0, f"{s2:.3f}")

        # --- Debug printout ---
        print(f"[Auto-shift EPOM] Det1={ic1.value}: raw={s1_raw:.3f} cm, res={res1:.3f} cm, corrected={s1:.3f} cm")
        print(f"[Auto-shift EPOM] Det2={ic2.value}: raw={s2_raw:.3f} cm, res={res2:.3f} cm, corrected={s2:.3f} cm")

    except Exception as e:
        print(f"[Auto-shift] failed: {e}")


def run_comparison():
   
    # Retrieve dd and dta values from the input fields Assumes Profile is in scale 100% and in cm
    dd = float(dd_entry.get())
    dta = float(dta_entry.get())/10
    print(f"Selected Dose Difference: {dd} ")
    print(f"Selected DTA Criteria: {dta} cm")
    global dl    
    global fsl
    global scl
    global report_text
    global fig
    fsl=[];scl=[];dl=[];   
    selected_indices = fsl_listbox.curselection()
    fsl.clear()  # Clear the current fsl list to avoid appending old values
    fsl.extend([float(fsl_listbox.get(i)) for i in selected_indices])
    fsl = [(x) for x in fsl]
    sclindexs=scl_listbox.curselection()
    scl.clear()
    scl.extend([scl_listbox.get(i) for i in sclindexs])
    # Reorder scl to ensure XY is processed last
    scl = sorted(scl, key=lambda x: (x == 'XY', x))  # This moves 'XY' to the end
       
    #fss=0;fsf=2;
    #dd=0.02;dta=.1;#dd should be 0.01 for 1%  Dta in cm so 0.1 = 1 mm
    
  
    #all of these set in gui intializing as 0
    dif=0;dist=0;comp=0;gam=0;plot=0;#iniates as unselected selected by gui
  

    
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
    elif selected_analysis == "plot":
        plot = 1
    
    # Set data processing variables based on user selection
   
    try:
        conv_fwhm = float(conv_fwhm_entry.get())
    except ValueError:
        conv_fwhm = 0.0

    try:
        ms = float(marker_size_entry.get())
    except ValueError:
        ms = 6
    
    conv_apply = conv_apply_var.get()
 
    smooth_option = smooth_combo.get()
    
    norm = int(norm_combo.get())
    
     
    if comp ==1 or gam==1:
        plt.rcParams.update({'font.size': 20})
        gs = gridspec.GridSpec(3,1,height_ratios=[1.5,1,1])
        fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 11), gridspec_kw={'height_ratios': [1.5, 1, 1]})
    elif dif==1 or dist==1:
        plt.rcParams.update({'font.size': 18})
        gs = gridspec.GridSpec(2,1,height_ratios=[1.5,1])
        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15, 9), gridspec_kw={'height_ratios': [1.5, 1]})
    else:
        plt.rcParams.update({'font.size': 16})
        #for no analysis
        gs = gridspec.GridSpec(1, 1)
        fig, ax0 = plt.subplots(1, 1, figsize=(10, 5))
    
    
    prtot=0;gvmax=0;gvmean=0
    total_points_cumulative = 0
    total_fails_cumulative = 0
    prl=[];gvl=[];gvtot=[];  comp_results=[]
    dif_total = 0; dif_fails = 0
    dist_total = 0; dist_fails = 0
    dif_results = []
    dist_results = []
    print("Comparison started...")
    df1=pd.read_excel(fn,sheet_name=sn1,header=0);
    df2=pd.read_excel(fn,sheet_name=sn2,header=0);   
    df1=df1.sort_values(by=['FS','Axis','Pos'])
    df2=df2.sort_values(by=['FS','Axis','Pos'])
    # Before starting the loop, check the type and contents of fsl
    print("Type of fsl:", type(fsl))
    print("Contents of fsl:", fsl)
    # Intialize plots
    ax0.plot([],'+r',ms=10,label=str(sn1))
    ax0.plot([],'.k',ms=10,label=sn2)
    ax0.legend()
    ax0.set_ylabel('Percentage Depth Dose [%]');
    ax0.set_xlabel(' Depth [cm]')
           
    for j in range(len(fsl)):#FS List  
        fsl[j]=float(fsl[j])    
        prfs=[];gvfs=[];
        # Extract dose and position data
        y1 = df1.loc[(df1['FS']==fsl[j]),'Pos'].reset_index(drop=True)
        d1 = df1.loc[(df1['FS']==fsl[j]),'Dose'].reset_index(drop=True)
        y2 = df2.loc[(df2['FS']==fsl[j]),'Pos'].reset_index(drop=True)
        d2 = df2.loc[(df2['FS']==fsl[j]),'Dose'].reset_index(drop=True)

        # Check for duplicate depth positions — skip FS and warn if found
        dup1 = y1[y1.duplicated(keep=False)]
        dup2 = y2[y2.duplicated(keep=False)]
        if not dup1.empty:
            print(f"  FS {fsl[j]}: WARNING — Sheet 1 has duplicate Pos values "
                  f"{sorted(dup1.unique().tolist())} — skipping.")
            continue
        if not dup2.empty:
            print(f"  FS {fsl[j]}: WARNING — Sheet 2 has duplicate Pos values "
                  f"{sorted(dup2.unique().tolist())} — skipping.")
            continue
        # Ensure strictly increasing depth order
        s1 = y1.argsort(); y1, d1 = y1.iloc[s1].reset_index(drop=True), d1.iloc[s1].reset_index(drop=True)
        s2 = y2.argsort(); y2, d2 = y2.iloc[s2].reset_index(drop=True), d2.iloc[s2].reset_index(drop=True)

        if auto_shift_var.get():
            autofill_depth_shifts_epom_ui(y1, d1, y2, d2, unit="cm")

        # Retrieve user-specified shifts
        try:
            depth_shift1 = float(depth_shift1_entry.get())
        except:
            depth_shift1 = 0.0
        
        try:
            depth_shift2 = float(depth_shift2_entry.get())
        except:
            depth_shift2 = 0.0
        
        # Apply shifts: Negative = shallower, Positive = deeper
        y1 = y1 + depth_shift1
        y2 = y2 + depth_shift2
        
        #print(f"Applied depth shifts: Curve 1 = {depth_shift1} cm, Curve 2 = {depth_shift2} cm")
        # ---- Optional cutoff filter: blank = off; numeric (including 0) = apply ----
        cutoff_txt = cutoff_depth_entry.get().strip()
        
        if cutoff_txt != "":
            try:
                cutoff = float(cutoff_txt)
            except:
                cutoff = None
                print(f"Cutoff Depth is not a number: '{cutoff_txt}'. No cutoff applied.")
        
            if cutoff is not None:
                m1 = (y1 >= cutoff)
                m2 = (y2 >= cutoff)
                y1 = y1[m1]
                d1 = d1[m1]
                y2 = y2[m2]
                d2 = d2[m2]
        
                # keep things clean for later interpolation/indexing
                try:
                    y1 = y1.reset_index(drop=True); d1 = d1.reset_index(drop=True)
                    y2 = y2.reset_index(drop=True); d2 = d2.reset_index(drop=True)
                except Exception:
                    pass
        
                if len(y1) < 2 or len(y2) < 2:
                    print(f"FS {fsl[j]:.1f} cm: cutoff {cutoff} cm removed too many points; skipping.")
                    continue


        # Normalization
        norm = int(norm_combo.get())
        
        if norm == 0:
            pass  # No normalization
        
        elif norm == 1:
            d1 = d1 / d1.max() *100
            d2 = d2 / d2.max() *100
        
        elif norm == 2:
            try:
                fixed_depth = float(fixed_depth_entry.get())
            except:
                fixed_depth = 1.5  # fallback
        
            # Dataset 1
            if fixed_depth < y1.min() or fixed_depth > y1.max():
                print(f"Warning: Fixed depth {fixed_depth} cm out of range for dataset 1. Using closest value.")
                closest_idx = (y1 - fixed_depth).abs().idxmin()
                d1_fixed = d1.loc[closest_idx]
            else:
                interp1 = interp.pchip(y1, d1)
                d1_fixed = interp1(fixed_depth)
        
            # Dataset 2
            if fixed_depth < y2.min() or fixed_depth > y2.max():
                print(f"Warning: Fixed depth {fixed_depth} cm out of range for dataset 2. Using closest value.")
                closest_idx = (y2 - fixed_depth).abs().idxmin()
                d2_fixed = d2.loc[closest_idx]
            else:
                interp2 = interp.pchip(y2, d2)
                d2_fixed = interp2(fixed_depth)
        
            d1 = d1 / d1_fixed *100
            d2 = d2 / d2_fixed *100
            
        if conv_fwhm > 0:
            
            if conv_apply == "Curve 1 Only":
                d1 = apply_detector_convolution(y1, d1, conv_fwhm)
            elif conv_apply == "Curve 2 Only":
                d2 = apply_detector_convolution(y2, d2, conv_fwhm)
            elif conv_apply == "Both Curves":
                d1 = apply_detector_convolution(y1, d1, conv_fwhm)
                d2 = apply_detector_convolution(y2, d2, conv_fwhm)

    
        # Apply smoothing based on user selection
        if smooth_option == "Curve 1 Only":
            d1 = sig.savgol_filter(d1, 3, 1)
        
        elif smooth_option == "Curve 2 Only":
            d2 = sig.savgol_filter(d2, 3, 1)
        
        elif smooth_option == "Both":
            d1 = sig.savgol_filter(d1, 3, 1)
            d2 = sig.savgol_filter(d2, 3, 1)
        
        # If "None", do nothing

       
        ax0.plot(y1,d1,'+r',ms=ms)
        ax0.plot(y2,d2,'.k',ms=ms)  
        if gam ==1:
            gx , gv = g(y1,d1,y2,d2,dd,dta,0,0,.01);#norm 0 is no norm in gamma code, 0 threshold InterpThreshold=0.01
            gx, gv = downsample_to_native(gx, gv, y1)
            gx1=np.extract(gv>1,gx);gv1=np.extract(gv>1,gv);
            gx2=np.extract(gv<1,gx);gv2=np.extract(gv<1,gv);
            pr=(np.size((np.where(np.asarray(gv)<=1)))/np.size(np.where(np.asarray(gv)>=0))*100)
            pr=np.round(pr,1)
            #print(pr,max(np.round(gv,2)),fsl[j],dl[i],scl[k])
            gvtot.extend(gv)
            gvfs.extend(gv)
            gvl.append(gvfs.copy())  # Save gamma values for this field size

            prfs=(np.size((np.where(np.asarray(gvfs)<=1)))/np.size(np.where(np.asarray(gvfs)>=0))*100)
            prtot=(np.size((np.where(np.asarray(gvtot)<=1)))/np.size(np.where(np.asarray(gvtot)>=0))*100)
            prtot=np.round(prtot,2)
            prl.append(pr)
            ax0.set_ylim(0, max(1.05 * np.max(d1), 1.05 * np.max(d2)))
            ax0.set_xlim(0,30)
            ax1.set_xlim(0,30)
            ax1.plot(gx1[1::1],gv1[1::1],'.r',ms=ms)
            ax1.plot(gx2[1::1],gv2[1::1],'.g',ms=ms)
            ax1.set_ylim(0, max(1.5, np.max(gv) + 0.1))
            ax1.set_ylabel(r'$\Gamma$'+'['+ str(dd)+'%/'+ str(dta*10) + 'mm]')
            ax1.set_xlabel( 'Points with ' + r'$\Gamma$' +' > 1 : '+str(np.size((np.where(np.asarray(gvtot)>1))))+'/'+str(len(gvtot)) + '  Pass Rate : '+ str(prtot)+'%' )
            gvmean=np.round(np.mean(gvtot),3);gvmax=np.round(np.max(gvtot),3); 
            bins=0.1;
            weights = np.ones_like(np.asarray(gvtot))/float(len(np.asarray(gvtot)))
            
            #plt.title('Gamma Histogram Data:  '  + str(np.round(100-sum(histd[0][0:(1/bins)]),1)) + '% of Data < Gamma = 1')
            ax2.set_xlim(0,1.5)
            ax2.set_xticks(np.arange(0,1.5,bins))
            ax2.set_xlabel(r'$\Gamma$'+ '['+ str(dd)+'%/'+ str(dta*10) + 'mm]')
            ax2.set_ylabel('Normalized Incidence')
            #print(pr,prtot,gvmean,gvmax,fsl[j],dl[i],scl[k])
    
        if dif ==1:
            
            difx, ddif= dosedif(y1,d1,y2,d2,0)#no normalization in dose dif
            difx, ddif = downsample_to_native(difx, ddif, y1)
            ddx1=np.extract(np.abs(ddif)>dd,difx);ddv1=np.extract(np.abs(ddif)>dd,ddif);
            ddx2=np.extract(np.abs(ddif)<dd,difx);ddv2=np.extract(np.abs(ddif)<dd,ddif);
            dif_total += len(ddif)
            dif_fails += len(ddx1)
            fs_pr = 100 * (len(ddif) - len(ddx1)) / len(ddif) if len(ddif) > 0 else 0.0
            dif_results.append((fsl[j], len(ddx1), len(ddif), fs_pr))
            ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.plot(ddx1,ddv1,'.r',ms=ms)
            ax1.plot(ddx2,ddv2,'.g',ms=ms)
            ax1.set_ylabel('Dose Difference [%]')
        if dist ==1:
            dtax, dtav= dtafunc(y1,d1,y2,d2,dta)
            dtax, dtav = downsample_to_native(dtax, dtav, y1)
            dtax1=np.extract(dtav>dta,dtax);dtav1=np.extract(dtav>dta,dtav);
            dtax2=np.extract(dtav<dta,dtax);dtav2=np.extract(dtav<dta,dtav);
            dist_total += len(dtav)
            dist_fails += len(dtax1)
            fs_dist_pr = 100 * (len(dtav) - len(dtax1)) / len(dtav) if len(dtav) > 0 else 0.0
            dist_results.append((fsl[j], len(dtax1), len(dtav), fs_dist_pr))
            ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.plot(dtax1,dtav1*10,'.r',ms=ms)
            ax1.plot(dtax2,dtav2*10,'.g',ms=ms)
            ax1.set_ylabel('DTA [mm]')
        if comp ==1:
          
            dtax, dtav= dtafunc(y1,d1,y2,d2,dta)
            difx, difv= dosedif(y1,d1,y2,d2,0)#No normalization in dose difference
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
            # Calculate mean values for this FS
            mean_dose_diff = np.mean(np.abs(difv)) if len(difv) > 0 else 0.0
            mean_dta_val   = np.mean(np.abs(dtav)) if len(dtav) > 0 else 0.0
          
            total=len(dtav)
            total_points_cumulative += total
            total_fails_cumulative += totalfail

            prtot=(total-totalfail)/total*100
            comp_results.append((fsl[j], prtot, mean_dose_diff, mean_dta_val*10,totalfail,total))

            print(f"FS {fsl[j]:.1f} cm : PassRate={prtot:.2f}%, MeanDoseDiff={mean_dose_diff:.2f}%, MeanDTA={mean_dta_val*10:.2f} mm, Failing points = {totalfail}/{total}")


            
            ax0.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax2.set_xlim(max(min(y1),min(y2)),min(max(y1),max(y2)))
            ax1.plot(dtax,dtav*10,'.g',ms=ms)
            ax1.set_ylabel('DTA [mm]')
            #ax2.set_xlabel(( 'Points outside of ' + str(dd)+'% & '+str(dta*10) +'mm ' +str(totalfail)+'/'+str(total) + '  Pass Rate : '+'%.2f' % prtot+'%' ))
            ax2.plot(difx,difv,'.g',ms=ms)#convert to %
            ax2.set_ylabel('Dose Difference [%] ')
            for n in totalfailloc:#plots the failing points
                ax1.plot(dtax[n],dtav[n]*10,'.r',ms=ms)
                ax2.plot(difx[n],difv[n],'.r',ms=ms)
                
        if gam ==1:
            print(f"FS {fsl[j]:<5}: PassRate={pr:.2f}%, Mean={np.mean(gvfs):.3f}, Max={max(gvfs):.3f}") 
#    if gam==1:
#        print('All FS: ','%.2f' % prtot,'%.3f' % np.mean(gvtot),'%.3f' % max(gvtot))
    
    if comp == 1:    
        report_lines = []
        report_lines.append("PDD Composite Analysis Results")
        report_lines.append(f"Criteria: {dd:.1f}% / {dta*10:.1f} mm")
        report_lines.append("---------------------------------------------------------------------------------")
        report_lines.append("Field Size [cm] | Pass Rate [%] | Mean Dose Diff [%] | Mean DTA [mm] | Failing Points / Total Points")
        report_lines.append("---------------------------------------------------------------------------------")
        
        for fs_val, passrate, mean_diff, mean_dta_val, fails, total in comp_results:
            report_lines.append(f"{fs_val:<15.1f} | {passrate:>12.2f} | {mean_diff:>18.2f} | {mean_dta_val:>12.2f} | {fails}/{total}")
        
        # Calculate overall sums for failing and total points
        total_fails_cumulative = sum(row[4] for row in comp_results)
        total_points_cumulative = sum(row[5] for row in comp_results)
        
        overall_pass = 100 * (1 - total_fails_cumulative / total_points_cumulative) if total_points_cumulative > 0 else 0
        overall_diff = np.mean([row[2] for row in comp_results])
        overall_dta  = np.mean([row[3] for row in comp_results])
        
        report_lines.append("---------------------------------------------------------------------------------")
        report_lines.append(f"{'Overall':<15} | {overall_pass:>12.2f} | {overall_diff:>18.2f} | {overall_dta:>12.2f} | {total_fails_cumulative}/{total_points_cumulative}")
        
        report_text = "\n".join(report_lines)
        print("\n" + report_text)
        
        # Update xlabel with the cumulative summary
        ax2.set_xlabel(
            f'Points outside of {dd:.1f}% & {dta*10:.1f}mm '
            f'{total_fails_cumulative}/{total_points_cumulative}  Pass Rate : {overall_pass:.2f}%'
        )
        figManager = plt.get_current_fig_manager()
        try:
            figManager.window.showMaximized()
        except AttributeError:
            figManager.window.state('zoomed')
        fig.tight_layout()
        fig.subplots_adjust(hspace=.3)

    if gam == 1:
        results = []
        
        # Collect results for each field size
        for j in range(len(fsl)):
            fs_val = float(fsl[j])
            prfs_val = prl[j] if j < len(prl) else 0.0
            current_gvfs = gvl[j] if j < len(gvl) else []
            mean_val = np.mean(current_gvfs) if len(current_gvfs) > 0 else 0.0
            max_val = max(current_gvfs) if len(current_gvfs) > 0 else 0.0
            results.append((fs_val, prfs_val, mean_val, max_val, current_gvfs))

        # Build report lines
        report_lines = []
        report_lines.append("Gamma Analysis Results")
        report_lines.append(f"Criteria: {dd:.1f}% / {dta*10:.1f} mm")
        header = "Field Size [cm] | Pass Rate [%] | Gamma Mean | Gamma Max | Failing Points / Total Points"
        sep = "-" * len(header)
        report_lines.append(sep)
        report_lines.append(header)
        report_lines.append(sep)
        
        # Add each FS result
        for fs_val, pr_val, mean_val, max_val, gvfs_data in results:
            fails = int(np.sum(np.asarray(gvfs_data) > 1))
            total = len(gvfs_data)
            report_lines.append(f"{fs_val:<15.1f} | {pr_val:>12.2f} | {mean_val:>10.3f} | {max_val:>9.3f} | {fails}/{total}")

        # Add overall result
        report_lines.append(sep)
        fails_overall = int(np.sum(np.asarray(gvtot) > 1))
        total_overall = len(gvtot)
        report_lines.append(f"{'Overall':<15} | {prtot:>12.2f} | {np.mean(gvtot):>10.3f} | {max(gvtot):>9.3f} | {fails_overall}/{total_overall}")
        # Convert list to one string (for report)
        report_text = "\n".join(report_lines)

        # Print to console
        print(report_text)


        if gam ==1:
            r=[0,1.5]
            histd=np.histogram(gvtot,bins=np.arange(0,np.max(gvtot),bins),weights=weights,range=r);
            ax2.hist(gvtot,bins=np.arange(0,np.max(gvtot),bins),weights=weights,range=r)
        figManager = plt.get_current_fig_manager()#this makes full screen by default
        try:
            figManager.window.showMaximized()
        except AttributeError:
            figManager.window.state('zoomed')
        fig.tight_layout()
        fig.subplots_adjust(hspace=.3)

    if dif == 1:
        dif_prtot = 100 * (dif_total - dif_fails) / dif_total if dif_total > 0 else 0
        header = "Field Size [cm] | Fails / Total | Pass Rate [%]"
        sep = "-" * len(header)
        print("\nDose Difference Analysis Results")
        print(f"Criteria: {dd:.1f}%")
        print(sep); print(header); print(sep)
        for fs_val, fails, total, pr in dif_results:
            print(f"{fs_val:<15.1f} | {fails:>5}/{total:<7} | {pr:>12.2f}")
        print(sep)
        print(f"{'Overall':<15} | {dif_fails:>5}/{dif_total:<7} | {dif_prtot:>12.2f}")
        ax1.set_xlabel(
            f'Points outside {dd:.1f}% : {dif_fails}/{dif_total}  Pass Rate : {dif_prtot:.2f}%'
        )
        figManager = plt.get_current_fig_manager()
        try:
            figManager.window.showMaximized()
        except AttributeError:
            figManager.window.state('zoomed')
        fig.tight_layout()
        fig.subplots_adjust(hspace=.3)

    if dist == 1:
        dist_prtot = 100 * (dist_total - dist_fails) / dist_total if dist_total > 0 else 0
        header = "Field Size [cm] | Fails / Total | Pass Rate [%]"
        sep = "-" * len(header)
        print("\nDTA Analysis Results")
        print(f"Criteria: {dta*10:.1f} mm")
        print(sep); print(header); print(sep)
        for fs_val, fails, total, pr in dist_results:
            print(f"{fs_val:<15.1f} | {fails:>5}/{total:<7} | {pr:>12.2f}")
        print(sep)
        print(f"{'Overall':<15} | {dist_fails:>5}/{dist_total:<7} | {dist_prtot:>12.2f}")
        ax1.set_xlabel(
            f'Points outside {dta*10:.1f} mm : {dist_fails}/{dist_total}  Pass Rate : {dist_prtot:.2f}%'
        )
        figManager = plt.get_current_fig_manager()
        try:
            figManager.window.showMaximized()
        except AttributeError:
            figManager.window.state('zoomed')
        fig.tight_layout()
        fig.subplots_adjust(hspace=.3)

    if plot == 1:
        ax0.set_xlim(max(min(y1), min(y2)), min(max(y1), max(y2)))
        ax0.set_ylim(0, max(1.05 * np.max(d1), 1.05 * np.max(d2)))
        ax0.set_ylabel('Percentage Depth Dose [%]')
        ax0.set_xlabel('Depth [cm]')
        ax0.legend([sn1, sn2])
        figManager = plt.get_current_fig_manager()
        try:
            figManager.window.showMaximized()
        except AttributeError:
            figManager.window.state('zoomed')
        fig.tight_layout()

    plt.pause(0.001)

    
# Initialize main window
root = tk.Tk()
root.title("PDD Compare ToolV1.0")

# Bind Tk variables to root (portable across Tcl/Tk builds)
analysis_var   = tk.StringVar(master=root, value="none")
auto_shift_var = tk.BooleanVar(master=root, value=False)
conv_apply_var = tk.StringVar(master=root, value="None")

# --- Scrollable wrapper (Canvas + vertical scrollbar) ---
def make_scrollable(master):
    container = ttk.Frame(master)
    container.pack(fill="both", expand=True)

    canvas = tk.Canvas(container, highlightthickness=0)
    vbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vbar.set)

    vbar.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    content = ttk.Frame(canvas)
    win_id = canvas.create_window((0, 0), window=content, anchor="nw")

    def _sync_scrollregion(_event=None):
        canvas.configure(scrollregion=canvas.bbox("all"))

    def _sync_width(_event=None):
        # Keep the embedded frame the same width as the canvas viewport
        canvas.itemconfigure(win_id, width=canvas.winfo_width())
        _sync_scrollregion()

    # When inner content changes size
    content.bind("<Configure>", _sync_scrollregion)
    # When the canvas viewport changes size (this is the key fix)
    canvas.bind("<Configure>", _sync_width)

    # Force one layout pass after everything is created
    canvas.after_idle(_sync_width)

    # Mouse wheel (Win/Mac/Linux) -- keep your existing code here if you want
    def _on_mousewheel(e):
        if e.delta:
            canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
        elif getattr(e, "num", None) == 4:
            canvas.yview_scroll(-1, "units")
        elif getattr(e, "num", None) == 5:
            canvas.yview_scroll(1, "units")
    canvas.bind_all("<MouseWheel>", _on_mousewheel)
    canvas.bind_all("<Button-4>", _on_mousewheel)
    canvas.bind_all("<Button-5>", _on_mousewheel)

    return content

# Use scrollable frame as the parent for all widgets below
main = make_scrollable(root)

# File selection section
file_frame = ttk.Frame(main, padding="10")
file_frame.grid(row=0, column=0, sticky="ew")
ttk.Label(file_frame, text="Select File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(file_frame, width=40)
file_entry.grid(row=0, column=1, padx=5)
ttk.Button(file_frame, text="Browse", command=select_file).grid(row=0, column=2, padx=5)

# Sheet selection section
sheet_frame = ttk.Frame(main, padding="10")
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
    populate_fsl()
    populate_scl()

ttk.Button(sheet_frame, text="Swap \u21c5", command=swap_sheets).grid(row=0, column=2, rowspan=2, padx=5)

# Data Selection Frame
data_selection_frame = ttk.LabelFrame(main, text="Data Selection", padding="10")
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
scl_listbox = tk.Listbox(data_selection_frame, selectmode="multiple", width=20, height=10, exportselection=0)
scl_listbox.grid(row=1, column=1, padx=5, pady=5)

# Analysis type section with radio buttons
analysis_frame = ttk.LabelFrame(main, text="Analysis Type", padding="10")
analysis_frame.grid(row=3, column=0, sticky="ew")
ttk.Radiobutton(analysis_frame, text="Raw Data Only",        variable=analysis_var, value="none").grid(row=0, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Plots Only",           variable=analysis_var, value="plot").grid(row=1, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Gamma",                variable=analysis_var, value="gam").grid(row=2, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Composite Analysis",   variable=analysis_var, value="comp").grid(row=3, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Dose Difference",      variable=analysis_var, value="dif").grid(row=4, column=0, sticky="w")
ttk.Radiobutton(analysis_frame, text="Distance to Agreement",variable=analysis_var, value="dist").grid(row=5, column=0, sticky="w")

# Data processing section with dropdowns and checkboxes
processing_frame = ttk.LabelFrame(main, text="Data Processing", padding="10")
processing_frame.grid(row=4, column=0, sticky="ew")

# Normalization dropdown
ttk.Label(processing_frame, text="Normalization:").grid(row=1, column=0, sticky="w")
norm_combo = ttk.Combobox(processing_frame, values=[0, 1, 2], width=5)
norm_combo.grid(row=1, column=1, padx=5)
norm_combo.set(1)
ttk.Label(processing_frame, text="0: None, 1: Normalize to dmax, 2: Normalize to Fixed Depth").grid(row=1, column=2, sticky="w")

# Fixed Depth Input for normalization option 2
ttk.Label(processing_frame, text="Fixed Depth [cm] for Norm:").grid(row=2, column=0, sticky="w")
fixed_depth_entry = ttk.Entry(processing_frame, width=10)
fixed_depth_entry.grid(row=2, column=1, padx=5)
fixed_depth_entry.insert(0, "1.5")

# Smoothing option dropdown
ttk.Label(processing_frame, text="Smoothing:").grid(row=3, column=0, sticky="w")
smooth_combo = ttk.Combobox(processing_frame, values=["None", "Curve 1 Only", "Curve 2 Only", "Both"], width=15)
smooth_combo.grid(row=3, column=1, padx=5)
smooth_combo.set("None")

# Auto-shift (above manual shift boxes)
ttk.Checkbutton(processing_frame, text="Auto-shift", variable=auto_shift_var).grid(row=4, column=0, columnspan=2, sticky="w")

# Detector / modality (used by Auto-shift)
ttk.Label(processing_frame, text="Detector 1:").grid(row=4, column=2, sticky="w")
det1_combo = ttk.Combobox(processing_frame, width=12, values=("PTW31021", "PTW31010", "IBA CC13"))
det1_combo.grid(row=4, column=3, padx=5)
det1_combo.set("PTW31021")

ttk.Label(processing_frame, text="Detector 2:").grid(row=5, column=2, sticky="w")
det2_combo = ttk.Combobox(processing_frame, width=12, values=("PTW31021", "PTW31010", "IBA CC13"))
det2_combo.grid(row=5, column=3, padx=5)
det2_combo.set("PTW31021")

ttk.Label(processing_frame, text="Modality:").grid(row=6, column=2, sticky="w")
modality_combo = ttk.Combobox(processing_frame, width=12, values=("PHOTON", "ELECTRON"))
modality_combo.grid(row=6, column=3, padx=5)
modality_combo.set("PHOTON")

# Depth Shift entries
ttk.Label(processing_frame, text="Depth Shift Curve 1 [cm]:").grid(row=5, column=0, sticky="w")
depth_shift1_entry = ttk.Entry(processing_frame, width=10)
depth_shift1_entry.grid(row=5, column=1, padx=5)
depth_shift1_entry.insert(0, "0")

ttk.Label(processing_frame, text="Depth Shift Curve 2 [cm]:").grid(row=6, column=0, sticky="w")
depth_shift2_entry = ttk.Entry(processing_frame, width=10)
depth_shift2_entry.grid(row=6, column=1, padx=5)
depth_shift2_entry.insert(0, "0")

# Detector Convolution
ttk.Label(processing_frame, text="Detector Convolution FWHM [cm]:").grid(row=7, column=0, sticky="w")
conv_fwhm_entry = ttk.Entry(processing_frame, width=10)
conv_fwhm_entry.grid(row=7, column=1, padx=5)
conv_fwhm_entry.insert(0, "0")

ttk.Label(processing_frame, text="Apply Convolution To:").grid(row=8, column=0, sticky="w")
conv_apply_combo = ttk.Combobox(processing_frame, width=15, textvariable=conv_apply_var,
                                values=("None", "Curve 1 Only", "Curve 2 Only", "Both Curves"))
conv_apply_combo.grid(row=8, column=1, padx=5)
conv_apply_combo.current(0)

# --- Cutoff depth filter ---
ttk.Label(processing_frame, text="Cutoff Depth [cm] (filter shallower):").grid(row=9, column=0, sticky="w")
cutoff_depth_entry = ttk.Entry(processing_frame, width=10)
cutoff_depth_entry.grid(row=9, column=1, padx=5)
cutoff_depth_entry.insert(0, 0)

# --- Marker size ---
ttk.Label(processing_frame, text="Marker Size:").grid(row=10, column=0, sticky="w")
marker_size_entry = ttk.Entry(processing_frame, width=10)
marker_size_entry.grid(row=10, column=1, padx=5)
marker_size_entry.insert(0, "6")


# Criteria frame
criteria_frame = ttk.LabelFrame(main, text="Dose Difference and DTA Criteria Settings: Default:1% and 1mm", padding="10")
criteria_frame.grid(row=5, column=0, sticky="ew")
ttk.Label(criteria_frame, text="Dose Difference (%):").grid(row=0, column=0, sticky="w")
dd_entry = ttk.Entry(criteria_frame, width=10)
dd_entry.grid(row=0, column=1, padx=5)
dd_entry.insert(0, "1")
ttk.Label(criteria_frame, text="DTA Criteria [mm]:").grid(row=1, column=0, sticky="w")
dta_entry = ttk.Entry(criteria_frame, width=10)
dta_entry.grid(row=1, column=1, padx=5)
dta_entry.insert(0, "1")

# Buttons
run_button = ttk.Button(main, text="Run Comparison", command=run_comparison)
run_button.grid(row=6, column=0, pady=10)
savefig_button = ttk.Button(main, text="Save Figure", command=lambda: save_current_figure())
savefig_button.grid(row=7, column=0, pady=5)
save_report_button = ttk.Button(main, text="Save Report", command=save_report)
save_report_button.grid(row=8, column=0, pady=5)


# Start the GUI loop
root.mainloop()
