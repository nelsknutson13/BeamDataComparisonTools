"""Shared parsing, grouping, and physics helpers for MCC-based dosimetry tools."""
import re
import os
import numpy as np

# ── Depth / dmax tables ────────────────────────────────────────────────────────
NORM_DEPTH_CM = {
    "6MV":    1.5, "6MV FFF":  1.4, "8MV FFF": 1.9,
    "10MV":   2.4, "10MV FFF": 2.2, "15MV":    3.0,
}
DMAX_DEPTH_CM = dict(NORM_DEPTH_CM)

def _dmax_for_energy(energy_label):
    return DMAX_DEPTH_CM.get(energy_label)

# ── Available grouping facets ──────────────────────────────────────────────────
available_facets = [
    "Energy", "Detector", "SSD_cm", "FieldSize_cm",
    "NominalDoseRate_MUmin", "Detector_SN", "Linac",
    "ScanType", "ProfileDepth_cm", "HVPair",
]

# ── Uncertainty constants ──────────────────────────────────────────────────────
HV_REL_SPEC = 0.001
HV_REL_STD  = HV_REL_SPEC / np.sqrt(3.0)
EM_REL_STD  = 0.0005

def _u_beta(beta):
    return beta * np.sqrt(2.0) * HV_REL_STD

def _u_pion_from_beta(R, beta):
    return np.abs((1.0 - R) / ((R - beta)**2)) * _u_beta(beta)

# ── Detector name cleaning ─────────────────────────────────────────────────────
_ROLE_SUFFIX = re.compile(r'\s+(?:Field|Reference|Ref|Monitor|Det)\s*$', re.IGNORECASE)
_DEMO_WORD   = re.compile(r'\bDemo\b\s*', re.IGNORECASE)
_PAREN_NUM   = re.compile(r'\s*\(\d+\)\s*$')
_LETTER_SUF  = re.compile(r'\s*-\s*[A-Z]\s*$', re.IGNORECASE)

def _clean_det_name(name):
    name = (name or '').strip()
    name = _DEMO_WORD.sub('', name)
    name = _PAREN_NUM.sub('', name)
    name = _LETTER_SUF.sub('', name)
    name = _ROLE_SUFFIX.sub('', name)
    return name.strip()

# ── Parsing ───────────────────────────────────────────────────────────────────
def _energy_label(energy_mv, filt):
    if energy_mv is None:
        return None
    return f"{energy_mv:g}MV" + (" FFF" if (filt or "").strip().upper() == "FFF" else "")

def parse_scan(block):
    def _m(pat): return re.search(pat, block)
    m = _m(r'DETECTOR_HV=([0-9.]+)');       voltage       = float(m.group(1)) if m else None
    m = _m(r'DETECTOR_NAME=([^\n\r]+)');     det_name      = m.group(1).strip() if m else ""
    m = _m(r'SSD=([0-9.]+)');               ssd_cm        = float(m.group(1))/10.0 if m else None
    m = _m(r'ENERGY=([0-9.]+)');            energy_mv     = float(m.group(1)) if m else None
    m = _m(r'FILTER=([^\n\r]+)');           filt          = m.group(1).strip().upper() if m else None
    m = _m(r'SCAN_CURVETYPE=([^\n\r]+)');   scan_type     = m.group(1).strip().upper() if m else "PDD"
    m = _m(r'SCAN_DEPTH=([0-9.]+)');        prof_depth_cm = float(m.group(1))/10.0 if m else None
    m = _m(r'DETECTOR_CALIBRATION=([0-9.Ee+-]+)'); det_cal = float(m.group(1)) if m else None
    m = _m(r'MEAS_TIME=([0-9.]+)');         meas_time     = float(m.group(1)) if m else None
    m = _m(r'EXPECTED_MAX_DOSE_RATE=([0-9.]+)'); nominal_mu_min = float(m.group(1))*100 if m else None
    m = _m(r'DETECTOR_SN=([^\n\r]+)');      det_sn        = m.group(1).strip() if m else None
    m = _m(r'LINAC=([^\n\r]+)');            linac         = m.group(1).strip() if m else None

    fs_cm = None
    mfx = re.search(r'FIELD_INPLANE\s*=\s*([0-9.]+)', block)
    mfy = re.search(r'FIELD_CROSSPLANE\s*=\s*([0-9.]+)', block)
    if not (mfx and mfy):
        mfx = mfx or re.search(r'REF_FIELD_INPLANE\s*=\s*([0-9.]+)', block)
        mfy = mfy or re.search(r'REF_FIELD_CROSSPLANE\s*=\s*([0-9.]+)', block)
    if mfx and mfy:
        fx, fy = float(mfx.group(1))/10.0, float(mfy.group(1))/10.0
        fs_cm = (2*fx*fy/(fx+fy)) if (fx+fy) else None

    rows = []
    dm = re.search(r'BEGIN_DATA(.*?)END_DATA', block, flags=re.S)
    if dm:
        for ln in dm.group(1).strip().splitlines():
            parts = [p for p in re.split(r'\s+', ln.strip()) if p]
            if len(parts) >= 2:
                try:
                    d_mm  = float(parts[0])
                    dr_gpm = float(parts[1])
                    charge = (dr_gpm/60.0 * meas_time / det_cal) if (det_cal and meas_time) else np.nan
                    rows.append((d_mm/10.0, charge, dr_gpm))
                except ValueError:
                    pass
    arr = np.array(rows) if rows else np.empty((0, 3), dtype=float)
    return (voltage, arr, det_name, ssd_cm, energy_mv, filt, fs_cm,
            nominal_mu_min, det_sn, linac, scan_type, prof_depth_cm)

def parse_mcc(filepath):
    with open(filepath, "r", errors="ignore") as f:
        text = f.read()
    blocks = [b.split('END_SCAN')[0] for b in re.split(r'BEGIN_SCAN\s+\d+\s*', text)[1:]]
    return [parse_scan(b) for b in blocks if "BEGIN_DATA" in b]

def parse_mcc_path(path, recurse=False):
    if os.path.isdir(path):
        scans, file_count = [], 0
        walker = os.walk(path) if recurse else [(path, [], os.listdir(path))]
        for root, _, files in walker:
            for fn in sorted(files):
                if fn.lower().endswith(".mcc"):
                    file_count += 1
                    try:
                        scans.extend(parse_mcc(os.path.join(root, fn)))
                    except Exception as e:
                        print(f"[skip] {fn}: {e}")
        print(f"No .mcc files found." if file_count == 0
              else f"Loaded {file_count} .mcc file(s); {len(scans)} scan(s).")
        return scans
    scans = parse_mcc(path)
    print(f"Loaded 1 .mcc file; {len(scans)} scan(s).")
    return scans

# ── Grouping helpers ───────────────────────────────────────────────────────────
def _format_facet_val(facet, val):
    if facet == "SSD_cm":
        txt = "?" if val is None or (isinstance(val, float) and not np.isfinite(val)) else f"{val:g} cm"
        return f"SSD={txt}"
    if facet == "FieldSize_cm":
        txt = "?" if val is None or (isinstance(val, float) and not np.isfinite(val)) else f"{val:g} cm"
        return f"FS={txt}"
    if facet == "NominalDoseRate_MUmin":
        return f"Nominal={'?' if val is None else f'{int(val)} MU/min'}"
    if facet == "Detector_SN":
        return f"SN={val if val not in (None, '') else '?'}"
    if facet == "Linac":
        return f"LINAC={val if val not in (None, '') else '?'}"
    if facet == "ProfileDepth_cm":
        txt = "?" if val is None or (isinstance(val, float) and not np.isfinite(val)) else f"{val:g} cm"
        return f"Depth={txt}"
    if facet == "HVPair":
        return val if val not in (None, "") else "?"
    return str(val) if val not in (None, "") else "?"

def _group_label(key, facets):
    return " — ".join(_format_facet_val(f, v) for f, v in zip(facets, key))

def _get_facet(key, facet, facets):
    try:
        return key[list(facets).index(facet)]
    except ValueError:
        return None

def _split_common_varying(keys, facets):
    if not keys:
        return "", {}
    n = len(facets)
    varying = [i for i in range(n) if len({k[i] for k in keys}) > 1]
    ref = keys[0]
    common_str = " — ".join(_format_facet_val(facets[i], ref[i]) for i in range(n) if i not in varying)
    per_key = {
        key: " — ".join(_format_facet_val(facets[i], key[i]) for i in varying) or "—"
        for key in keys
    }
    return common_str, per_key

def _expand_by_hv_pair(groups, facets):
    if "HVPair" not in facets:
        return groups
    hv_pos = list(facets).index("HVPair")
    new_groups = {}
    for key, by_v in groups.items():
        uniq_vs = sorted(by_v.keys())
        hv_pairs = [(v, u) for v in uniq_vs for u in uniq_vs if abs(u - 2*v) <= 1e-6]
        if not hv_pairs:
            new_groups[key[:hv_pos] + ("?",) + key[hv_pos+1:]] = by_v
        else:
            for v_low, v_high in hv_pairs:
                tag = f"{int(v_high)}/{int(v_low)} V"
                new_groups[key[:hv_pos] + (tag,) + key[hv_pos+1:]] = {
                    v_low: by_v.get(v_low, []), v_high: by_v.get(v_high, []),
                }
    return new_groups

def _group_scans(scans, facets):
    groups = {}
    for item in scans:
        v, arr, det, ssd, energy_mv, filt = item[:6]
        fs_cm          = item[6]  if len(item) >= 7  else None
        nominal_mu_min = item[7]  if len(item) >= 8  else None
        det_sn         = item[8]  if len(item) >= 9  else None
        linac          = item[9]  if len(item) >= 10 else None
        scan_type      = item[10] if len(item) >= 11 else "PDD"
        prof_depth_cm  = item[11] if len(item) >= 12 else None
        meta = {
            "Energy": _energy_label(energy_mv, filt), "Detector": _clean_det_name(det),
            "SSD_cm": ssd, "FieldSize_cm": fs_cm, "NominalDoseRate_MUmin": nominal_mu_min,
            "Detector_SN": det_sn, "Linac": linac, "ScanType": scan_type,
            "ProfileDepth_cm": prof_depth_cm, "HVPair": None,
        }
        if isinstance(meta["SSD_cm"], float) and np.isfinite(meta["SSD_cm"]):
            meta["SSD_cm"] = round(meta["SSD_cm"], 3)
        if isinstance(meta["FieldSize_cm"], float) and np.isfinite(meta["FieldSize_cm"]):
            meta["FieldSize_cm"] = round(meta["FieldSize_cm"], 3)
        if isinstance(meta["NominalDoseRate_MUmin"], (int, float)) and np.isfinite(float(meta["NominalDoseRate_MUmin"])):
            meta["NominalDoseRate_MUmin"] = int(round(meta["NominalDoseRate_MUmin"]))
        key = tuple(meta[f] for f in facets) + (meta["ScanType"],)
        groups.setdefault(key, {}).setdefault(float(v) if v is not None else np.nan, []).append(arr)
    return groups

# ── Alignment ─────────────────────────────────────────────────────────────────
def _align_pair(arr_low, arr_high, tol=1e-4):
    dL, yL = arr_low[:,0], arr_low[:,1]
    dH, yH = arr_high[:,0], arr_high[:,1]
    iL = np.argsort(dL); dL, yL = dL[iL], yL[iL]
    iH = np.argsort(dH); dH, yH = dH[iH], yH[iH]
    if len(dL) == len(dH) and np.allclose(dL, dH, atol=tol):
        return dL, yL, yH
    base_d = dL if (np.median(np.diff(dL)) if len(dL)>1 else np.inf) <= (np.median(np.diff(dH)) if len(dH)>1 else np.inf) else dH
    lo, hi = max(dL.min(), dH.min()), min(dL.max(), dH.max())
    base_d = base_d[(base_d >= lo) & (base_d <= hi)]
    return base_d, np.interp(base_d, dL, yL), np.interp(base_d, dH, yH)

def align_all_voltages(by_v):
    """Average repeats per voltage, then align all voltages to a common depth grid.
    Returns (depths, {voltage: charge_array}, {voltage: sem_array})
    or (None, {}, {}) if insufficient data.
    SEM is the standard error of the mean across repeats; zero if only one repeat.
    """
    averaged = {}
    for v, arr_list in by_v.items():
        if not arr_list:
            continue
        ref_d = arr_list[0][:,0]
        iref = np.argsort(ref_d); ref_d = ref_d[iref]
        interped = [np.interp(ref_d, np.sort(a[:,0]), a[np.argsort(a[:,0]),1]) for a in arr_list]
        m = np.nanmean(interped, axis=0)
        n = len(interped)
        sem = (np.nanstd(interped, axis=0, ddof=1) / np.sqrt(n)) if n > 1 else np.zeros_like(m)
        averaged[v] = (ref_d, m, sem)

    if len(averaged) < 2:
        return None, {}, {}

    lo = max(d[0]  for d, _, _ in averaged.values())
    hi = min(d[-1] for d, _, _ in averaged.values())
    step = min((np.median(np.diff(d)) if len(d)>1 else np.inf) for d, _, _ in averaged.values())
    if not np.isfinite(step) or lo >= hi:
        return None, {}, {}

    base_d = sorted(averaged.items(), key=lambda kv: np.median(np.diff(kv[1][0])) if len(kv[1][0])>1 else np.inf)[0][1][0]
    base_d = base_d[(base_d >= lo) & (base_d <= hi)]

    aligned     = {v: np.interp(base_d, d, m)   for v, (d, m, sem) in averaged.items()}
    aligned_sem = {v: np.interp(base_d, d, sem) for v, (d, m, sem) in averaged.items()}
    return base_d, aligned, aligned_sem

# ── Pion math ─────────────────────────────────────────────────────────────────
def calc_pion(depths_cm, m_low, m_high, v_low, v_high):
    mratio = m_high / m_low
    vratio = v_high / v_low
    return (1 - vratio) / (mratio - vratio)

# ── Jaffe analysis ────────────────────────────────────────────────────────────
def jaffe_fit(depths, aligned_charges):
    """
    Fit 1/M = a + b*(1/V) at every depth using all available voltages.

    Parameters
    ----------
    depths : array (n_depths,)
    aligned_charges : dict {voltage_V: charge_array (n_depths,)}

    Returns
    -------
    M_inf : array  — true ionisation extrapolated to infinite voltage
    alpha : array  — slope (recombination coefficient, V·M units)
    R_sq  : array  — goodness-of-fit per depth
    """
    voltages = np.array(sorted(aligned_charges.keys()))
    inv_V = 1.0 / voltages
    n = len(depths)
    M_inf = np.full(n, np.nan)
    alpha = np.full(n, np.nan)
    R_sq  = np.full(n, np.nan)

    for i in range(n):
        inv_M = np.array([1.0 / aligned_charges[v][i] for v in voltages])
        ok = np.isfinite(inv_M) & np.isfinite(inv_V) & (inv_M != 0)
        if ok.sum() < 2:
            continue
        coeffs = np.polyfit(inv_V[ok], inv_M[ok], 1)  # [slope b, intercept a]
        b, a = coeffs
        if a > 0:
            M_inf[i] = 1.0 / a
        alpha[i] = b
        y_pred = np.polyval(coeffs, inv_V[ok])
        ss_res = np.sum((inv_M[ok] - y_pred)**2)
        ss_tot = np.sum((inv_M[ok] - inv_M[ok].mean())**2)
        R_sq[i] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return M_inf, alpha, R_sq

def jaffe_pion(M_inf, alpha, v_op):
    """
    Pion at operating voltage v_op from Jaffe fit parameters.
    pion(V) = M_inf / M(V) = 1 + alpha * M_inf / V
    """
    with np.errstate(invalid='ignore'):
        return 1.0 + alpha * M_inf / v_op
