# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 11:51:49 2026
Parses Dicom data out of blue ethos downloaded results.

@author: nknutson
"""

import os
import re
from glob import glob
from tkinter import Tk
from tkinter.filedialog import askdirectory

import numpy as np
import pydicom
import matplotlib.pyplot as plt



def ethos_stacked_mlc_metrics_from_rp(rp, beam_index=0, cp_index=0,dy_sample_mm=0.05, gap_min_mm=0.1):
    """
    Ethos/Enhanced RT Plan support:
      - Uses EnhancedRTBeamLimitingDeviceSequence for leaf boundaries (per stack)
      - Uses EnhancedRTBeamLimitingOpeningSequence for leaf positions (per stack), keyed by ReferencedDeviceIndex
    Handles double-stacked MLC with different leaf boundaries by resampling to a common Y grid.

    Returns dict:
      x_width_mm, y_width_mm, area_mm2, n_samples, device_indices_used
    """

    beam = rp.BeamSequence[beam_index]
    cp   = beam.ControlPointSequence[cp_index]

    if "EnhancedRTBeamLimitingDeviceSequence" not in beam:
        raise ValueError("Beam has no EnhancedRTBeamLimitingDeviceSequence (Ethos leaf geometry missing).")
    if "EnhancedRTBeamLimitingOpeningSequence" not in cp:
        raise ValueError("Control point has no EnhancedRTBeamLimitingOpeningSequence (Ethos leaf positions missing).")

    dev_seq = beam.EnhancedRTBeamLimitingDeviceSequence
    op_seq  = cp.EnhancedRTBeamLimitingOpeningSequence

    # ---- Build geometry map: DeviceIndex -> boundaries ----
    geom = {}
    for dev in dev_seq:
        if "DeviceIndex" not in dev:
            continue
        dev_idx = int(dev.DeviceIndex)

        # Ethos stores the leaf boundaries under ParallelRTBeamDelimiterDeviceSequence
        if "ParallelRTBeamDelimiterDeviceSequence" not in dev:
            continue
        p = dev.ParallelRTBeamDelimiterDeviceSequence[0]
        if "ParallelRTBeamDelimiterBoundaries" not in p:
            continue
        bounds = np.array(p.ParallelRTBeamDelimiterBoundaries, dtype=float)  # mm
        geom[dev_idx] = bounds

    if not geom:
        raise ValueError("Could not find any leaf boundaries in EnhancedRTBeamLimitingDeviceSequence.")

    # ---- Build openings map: ReferencedDeviceIndex -> (A,B) per leaf ----
    openings = {}
    for op in op_seq:
        if "ReferencedDeviceIndex" not in op or "ParallelRTBeamDelimiterPositions" not in op:
            continue
        dev_idx = int(op.ReferencedDeviceIndex)
        pos = np.array(op.ParallelRTBeamDelimiterPositions, dtype=float)  # mm
        n = pos.size // 2
        A = pos[:n]      # left edges (mm)
        B = pos[n:]      # right edges (mm)
        openings[dev_idx] = (A, B)

    # Keep only device indices that have BOTH geometry + openings
    dev_ids = [d for d in openings.keys() if d in geom]
    if len(dev_ids) == 0:
        raise ValueError("No matching device indices between geometry and openings.")

    # In your file you have two stacks (dev_ids should be length 2)
    # If there are more, we’ll intersect across all of them.
    # ---- Construct a common Y grid over the overlap of all stacks ----
    ymin = max(np.min(geom[d]) for d in dev_ids)
    ymax = min(np.max(geom[d]) for d in dev_ids)
    if ymax <= ymin:
        raise ValueError("Stacked MLC leaf boundary ranges do not overlap.")

    y = np.arange(ymin, ymax + dy_sample_mm, dy_sample_mm)  # mm

    def piecewise_edges_at_y(bounds, A, B, y_query):
        """
        bounds: (N+1) mm
        A,B: (N) mm leaf edges for each leaf interval
        returns left(y), right(y) arrays
        """
        # leaf interval index for each y (clipped)
        idx = np.searchsorted(bounds, y_query, side="right") - 1
        idx = np.clip(idx, 0, len(A)-1)
        return A[idx], B[idx]

    # ---- Intersect across stacks on the common Y grid ----
    left_eff  = None
    right_eff = None

    for d in dev_ids:
        bounds = geom[d]
        A, B = openings[d]
        # Sanity: A/B length should match leaf intervals count
        n_leaf = len(bounds) - 1
        if len(A) != n_leaf or len(B) != n_leaf:
            # If mismatch, still try a best-effort by clipping to min length
            m = min(len(A), len(B), n_leaf)
            bounds = bounds[:m+1]
            A = A[:m]
            B = B[:m]

        left_d, right_d = piecewise_edges_at_y(bounds, A, B, y)

        if left_eff is None:
            left_eff  = left_d
            right_eff = right_d
        else:
            left_eff  = np.maximum(left_eff, left_d)
            right_eff = np.minimum(right_eff, right_d)

    gap = right_eff - left_eff
    open_mask = gap > gap_min_mm

    if not np.any(open_mask):
        return {
            "device_indices_used": dev_ids,
            "x_width_mm": 0.0,
            "y_width_mm": 0.0,
            "area_mm2": 0.0,
            "n_samples": int(y.size),
        }

    x_min = float(np.min(left_eff[open_mask]))
    x_max = float(np.max(right_eff[open_mask]))
    x_width = x_max - x_min

    y_open = y[open_mask]
    y_width = float(np.max(y_open) - np.min(y_open))

    # Area by integrating gap(y) dy over open region (mm^2)
    area = float(np.sum(gap[open_mask]) * dy_sample_mm)

    return {
        "device_indices_used": dev_ids,
        "x_width_mm": x_width,
        "y_width_mm": y_width,
        "area_mm2": area,
        "n_samples": int(y.size),
    }
    

# -----------------------------
# Select directory containing the 3 files
# -----------------------------
Tk().withdraw()
root = askdirectory(title="Select CalculationResult folder (contains RP*, RD*, and CalculationLog)")
if not root:
    raise SystemExit("No folder selected.")

rp_files = sorted(glob(os.path.join(root, "RP*.dcm")))
rd_files = sorted(glob(os.path.join(root, "RD*.dcm")))
log_files = sorted(glob(os.path.join(root, "*CalculationLog*.txt")))

if not rp_files:
    raise FileNotFoundError("No RP*.dcm found in selected folder.")
if not rd_files:
    raise FileNotFoundError("No RD*.dcm found in selected folder.")

rp_path = rp_files[0]
rd_path = rd_files[0]
log_path = log_files[0] if log_files else None

print("RP:", os.path.basename(rp_path))
print("RD:", os.path.basename(rd_path))
if log_path:
    print("LOG:", os.path.basename(log_path))

# -----------------------------
# Read Plan (RP)
# -----------------------------
rp = pydicom.dcmread(rp_path)
m = ethos_stacked_mlc_metrics_from_rp(rp, beam_index=0, cp_index=0)

print("Device indices used:", m["device_indices_used"])
print(f"Effective MLC field bbox (intersection): "
      f"{m['x_width_mm']/10:.3f} x {m['y_width_mm']/10:.3f} cm")
print(f"Effective area: {m['area_mm2']/100:.3f} cm^2 (samples={m['n_samples']})")
beam = rp.BeamSequence[0]
cp0  = beam.ControlPointSequence[0]

plan_label = getattr(rp, "RTPlanLabel", "")
plan_name  = getattr(rp, "RTPlanName", "")
beam_name  = getattr(beam, "BeamName", "")
nominal_E  = float(getattr(cp0, "NominalBeamEnergy", np.nan))

# MU is in FractionGroupSequence -> ReferencedBeamSequence
mu = float(rp.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset)

# FFF flag lives in Primary Fluence Mode Sequence (Ethos/Varian often uses this)
fluence_mode = None
if "PrimaryFluenceModeSequence" in rp:
    f = rp.PrimaryFluenceModeSequence[0]
    fluence_mode = f"{getattr(f, 'FluenceMode', '')}/{getattr(f, 'FluenceModeID', '')}"

gantry = float(getattr(cp0, "GantryAngle", np.nan))
couch  = float(getattr(cp0, "PatientSupportAngle", np.nan))
coll   = float(getattr(cp0, "BeamLimitingDeviceAngle", np.nan))

print("\n--- PLAN ---")
print("Plan:", plan_label, "|", plan_name)
print("Beam:", beam_name, "| MU:", mu)
print("Nominal energy (MV):", nominal_E, "| FluenceMode:", fluence_mode)
print("Gantry/Couch/Coll (deg):", gantry, couch, coll)



# -----------------------------
# Read Dose (RD)
# -----------------------------
rd = pydicom.dcmread(rd_path)
dose = rd.pixel_array.astype(np.float32) * np.float32(rd.DoseGridScaling)

rows = int(rd.Rows)
cols = int(rd.Columns)
nz   = int(rd.NumberOfFrames)

dose_vol = dose.reshape((nz, rows, cols))  # (z,y,x)

dx_mm = float(rd.PixelSpacing[1])
dy_mm = float(rd.PixelSpacing[0])
gvec  = np.array(rd.GridFrameOffsetVector, dtype=np.float32)
dz_mm = float(np.mean(np.diff(gvec)))

print("\n--- DOSE ---")
print("Grid (mm):", dx_mm, dy_mm, dz_mm)
print("Dims (z,y,x):", dose_vol.shape)

# -----------------------------
# USER SETTINGS
# -----------------------------
phantomdepth = 35              # cm  (set to your actual phantom depth)
origindepth  = phantomdepth / 2.0 # cm
overscan_cm  = 10.0               # cm (± around center)
depths_req_cm = [5, 10, 20, 30]   # cm in phantom

# centers
cx = cols // 2
cy = rows // 2
cz = nz   // 2

# coordinates (cm)
x_cm_full = (np.arange(cols) - cx) * (dx_mm / 10.0)

z_cm_full = (gvec - gvec[cz]) / 10.0  # FORCE GridFrameOffsetVector
x_mask = np.abs(x_cm_full) <= overscan_cm
z_mask = np.abs(z_cm_full) <= overscan_cm

x_idx = np.where(x_mask)[0]
z_idx = np.where(z_mask)[0]
x_cm  = x_cm_full[x_mask]
z_cm  = z_cm_full[z_mask]

def depth_cm_to_y_index(depth_cm):
    y_rel_cm = depth_cm - origindepth
    return int(round(cy + (y_rel_cm * 10.0 / dy_mm)))

depth_pairs = []
for d in depths_req_cm:
    yi = depth_cm_to_y_index(d)
    if 0 <= yi < rows:
        d_actual = (yi - cy) * (dy_mm / 10.0) + origindepth
        depth_pairs.append((float(d), float(d_actual), yi))
    else:
        print(f"Skipping depth {d:.2f} cm (y-index {yi} out of bounds)")

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(20, 5))

plt.subplot(1, 3, 1)
for _, d_act, yi in depth_pairs:
    prof = dose_vol[cz, yi, x_idx]
    plt.plot(x_cm, prof, 'o', ms=3, label=f"{d_act:.2f} cm")
plt.xlabel("X (cm)"); plt.ylabel("Dose (Gy)")
plt.title(f"X profiles (±{overscan_cm:g} cm)")
plt.grid(True); plt.legend(fontsize=8)

plt.subplot(1, 3, 2)
for _, d_act, yi in depth_pairs:
    prof = dose_vol[z_idx, yi, cx]
    plt.plot(z_cm, prof, 'o', ms=3, label=f"{d_act:.2f} cm")
plt.xlabel("Z (cm)"); plt.ylabel("Dose (Gy)")
plt.title(f"Z profiles (±{overscan_cm:g} cm)")
plt.grid(True); plt.legend(fontsize=8)

plt.subplot(1, 3, 3)
y_cm_full = (np.arange(rows) - cy) * (dy_mm / 10.0) + origindepth
y_mask = (y_cm_full >= 0.0) & (y_cm_full <= phantomdepth)
plt.plot(y_cm_full[y_mask], dose_vol[cz, y_mask, cx], 'o', ms=3)
plt.xlabel("Depth in phantom (cm)"); plt.ylabel("Dose (Gy)")
plt.title(f"Central axis depth scan (0–{phantomdepth:g} cm)")
plt.grid(True)

plt.tight_layout()
plt.show()


