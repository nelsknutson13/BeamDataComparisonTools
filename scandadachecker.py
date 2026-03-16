import os
import re
import threading
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import sys
import pandas as pd

CACHE_FILE = Path.home() / ".scan_inventory_cache.pkl"

def load_cache():
    try:
        if CACHE_FILE.exists():
            with open(CACHE_FILE, "rb") as f:
                return pickle.load(f)
    except Exception:
        pass
    return {}

def save_cache(cache):
    try:
        with open(CACHE_FILE, "wb") as f:
            pickle.dump(cache, f)
    except Exception:
        pass

# ---------- Core logic ----------

def pick_col(cols, keywords):
    for c in cols:
        cl = str(c).lower()
        if any(kw in cl for kw in keywords):
            return c
    return None

def norm_scan_type(s):
    if pd.isna(s):
        return ""
    t = str(s).strip().upper()
    mapping = {"CROSSLINE":"X", "X-AXIS":"X", "INLINE":"Y", "Y-AXIS":"Y", "DIAGONAL":"XY", "DIAG":"XY"}
    return mapping.get(t, t)
# NEW: robust parse from arbitrary text (e.g., filename)
def parse_energy_ssd_from_string(text):
    s = str(text).upper().replace("\\", "/")
    s_compact = s.replace("_", "").replace(" ", "")

    # Energy: allow a digit after (e.g., 6FFF90SSD...), just forbid a letter (so we don't eat "6FFFX")
    m_e = re.search(r'(?<!\d)(\d{1,2}(?:X|FFF))(?![A-Z])', s_compact)
    energy = m_e.group(1) if m_e else None

    # SSD: support "100SSD", "100 SSD", "SSD100", "SSD 100"
    m_ssd = re.search(r'(?<!\d)(\d{2,3})\s*SSD(?!\d)', s) or re.search(r'SSD\s*(\d{2,3})', s)
    if not m_ssd:
        # fallback to common SSD values if explicit "SSD" pattern isn't present
        m_ssd = re.search(r'(?<!\d)(80|85|88|90|95|96|98|100|105|110|115|120)(?!\d)', s_compact)

    ssd = int(m_ssd.group(1)) if m_ssd else None
    return energy, ssd

def parse_energy_ssd_from_path(pathlike):
    s = str(pathlike).upper().replace("\\", "/")
    m_e = re.search(r'(?<!\d)(\d{1,2}X|\d{1,2}FFF)(?!\d)', s)
    energy = m_e.group(1) if m_e else None
    m_ssd = re.search(r'(?<!\d)(80|85|88|90|95|96|98|100|105|110|115|120)(?!\d)', s)
    ssd = int(m_ssd.group(1)) if m_ssd else None
    return energy, ssd

def parse_energy_ssd_from_sheet(sheet_name):
    s = sheet_name.upper().replace("_", "").replace(" ", "")
    m_e = re.search(r'(\d{1,2}X|\d{1,2}FFF)', s)
    energy = m_e.group(1) if m_e else None
    m_ssd = re.search(r'(80|85|88|90|95|96|98|100|105|110|115|120)', s)
    ssd = int(m_ssd.group(1)) if m_ssd else None
    return energy, ssd

def get_energy_ssd(xlsx_path, sheet_name=None, df=None):
    """Return (energy, ssd) with precedence:
       1. Data columns in sheet
       2. Sheet name
       3. File name
       4. Folder path
    """
    energy_dir, ssd_dir   = parse_energy_ssd_from_path(xlsx_path.parent)
    energy_file, ssd_file = parse_energy_ssd_from_string(xlsx_path.name)

    energy_cell, ssd_cell = None, None
    if df is not None:
        if "Energy" in df.columns:
            v = df["Energy"].dropna()
            if not v.empty:
                energy_cell = str(v.iloc[0]).strip().upper()
        if "SSD" in df.columns:
            v = pd.to_numeric(df["SSD"], errors="coerce").dropna()
            if not v.empty:
                ssd_cell = int(v.iloc[0])

    e_sheet, s_sheet = (None, None)
    if sheet_name:
        e_sheet, s_sheet = parse_energy_ssd_from_sheet(sheet_name)

    energy = energy_cell or e_sheet or energy_file or energy_dir
    ssd    = ssd_cell    or s_sheet  or ssd_file   or ssd_dir
    return energy, ssd

def parse_device(sheet_name):
    s = sheet_name.upper().replace(" ", "")
    m = re.search(r'(TPS-?SN\d+|SN\d+)', s)
    return m.group(1) if m else sheet_name

def load_requirements(want_path):
    xls = pd.ExcelFile(want_path)
    tabs = []
    for sh in xls.sheet_names:
        df = xls.parse(sh)
        df.columns = [str(c).strip() for c in df.columns]
        # If an Energy column is missing or blank, fill with sheet name
        if 'Energy' not in df.columns:
            df['Energy'] = sh
        else:
            df['Energy'] = df['Energy'].where(df['Energy'].notna() & (df['Energy'].astype(str).str.strip()!=''), sh)
        tabs.append(df)
    want = pd.concat(tabs, ignore_index=True)
    colmap = {}
    for c in want.columns:
        cl = str(c).lower()
        if cl in ("energy","beam","modality"):
            colmap[c] = "Energy"
        elif "ssd" in cl:
            colmap[c] = "SSD"
        elif cl in ("scan","scans","axis","profiletype","type"):
            colmap[c] = "ScanType"
    want = want.rename(columns=colmap)
    base_cols = [c for c in ["Energy","SSD","ScanType"] if c in want.columns]
    if base_cols:
        want = want.dropna(subset=base_cols, how="all")
    device_cols = [c for c in want.columns if c not in base_cols]
    return want, device_cols, base_cols

def inventory_for_excel(xlsx_path):
    rows = []
    try:
        xls = pd.ExcelFile(xlsx_path)
    except Exception as e:
        return pd.DataFrame([{"File": str(xlsx_path), "Error": str(e)}])

    for sh in xls.sheet_names:
        try:
            df = xls.parse(sh).dropna(how="all")
        except Exception as e:
            rows.append({"File": str(xlsx_path), "Device": parse_device(sh), "Error": str(e)})
            continue

        cols = list(df.columns)
        energy, ssd = get_energy_ssd(xlsx_path, sheet_name=sh, df=df)
        device = parse_device(sh)

        axis_col  = pick_col(cols, ["axis", "scan", "profiletype", "type"])
        fs_col    = pick_col(cols, ["fs", "field", "field size"])
        depth_col = pick_col(cols, ["depth", "z"])

        if axis_col and axis_col in df:
            scan_types = set(df[axis_col].dropna().astype(str).map(norm_scan_type).unique())
        else:
            scan_types = {"Z"}  # default to depth scan

        for st in scan_types:
            if axis_col and axis_col in df:
                sub = df[df[axis_col].map(norm_scan_type) == st]
            else:
                sub = df.copy()

            # Field sizes & depths
            fs_vals = sorted(set(sub[fs_col].dropna().unique())) if fs_col and fs_col in sub else []
            depth_vals = sorted(set(sub[depth_col].dropna().unique())) if depth_col and depth_col in sub else []

            def _to_scalar_list(vals):
                out = []
                for v in vals:
                    try:
                        f = float(v)
                        out.append(int(f) if f.is_integer() else f)
                    except Exception:
                        pass
                return out

            fs_vals    = _to_scalar_list(fs_vals)
            depth_vals = _to_scalar_list(depth_vals)

            # Count scans
            if fs_col and fs_col in sub:
                if st == "Z":
                    num_scans = sub[fs_col].dropna().nunique()
                else:
                    if depth_col and depth_col in sub:
                        sub2 = sub.dropna(subset=[fs_col, depth_col])
                        unique_pairs = set(zip(
                            pd.to_numeric(sub2[fs_col], errors="coerce"),
                            pd.to_numeric(sub2[depth_col], errors="coerce")
                        ))
                        num_scans = len(unique_pairs)
                    else:
                        num_scans = sub[fs_col].dropna().nunique()
            else:
                num_scans = 0

            rows.append({
                "File": str(xlsx_path),
                "Device": device,
                "Energy": energy,
                "SSD": ssd,
                "ScanType": st,
                "FieldSizes": fs_vals,
                "Depths": [] if st == "Z" else depth_vals,
                "NumScans": num_scans
            })

    inv_df = pd.DataFrame(rows)
    if inv_df.empty:
        return inv_df

    def merge_lists(series):
        vals = set()
        for v in series:
            if isinstance(v, list):
                vals.update(v)
        return sorted(vals)

    
    return (
        inv_df.groupby(["File", "Device", "Energy", "SSD", "ScanType"], dropna=False, as_index=False)
              .agg({"FieldSizes": merge_lists, "Depths": merge_lists, "NumScans": "sum"})
    )




def build_coverage(inventory, want, device_cols, base_cols):
    report_rows = []
    for _, row in want.iterrows():
        energy = str(row.get("Energy")).upper().strip() if "Energy" in row and pd.notna(row.get("Energy")) else None
        ssd_val = row.get("SSD")
        try:
            ssd = int(float(ssd_val)) if pd.notna(ssd_val) else None
        except Exception:
            ssd = None
        scan = norm_scan_type(row.get("ScanType")) if "ScanType" in want.columns else None
        for dev_col in device_cols:
            want_flag = str(row.get(dev_col)).strip().lower() if pd.notna(row.get(dev_col)) else ""
            required = want_flag in ("done","y","yes","x","required","1","true")
            if not required:
                continue
            dev_norm = str(dev_col).upper().replace(" ", "")
            mask = inventory["Device"].str.upper().str.contains(dev_norm, na=False)
            if energy:
                mask &= inventory["Energy"].fillna("").str.upper().str.contains(energy, na=False)
            if ssd:
                mask &= (inventory["SSD"] == ssd)
            candidates = inventory[mask]
            found = False
            where = ""
            if not candidates.empty:
                if scan:
                    for _, c in candidates.iterrows():
                        scans = (c.get("ScanTypes") or "")
                        if scan in scans.split(","):
                            found = True
                            where = f'{Path(c["File"]).name}::{c["Sheet"]}'
                            break
                else:
                    found = True
                    first = candidates.iloc[0]
                    where = f'{Path(first["File"]).name}::{first["Sheet"]}'
            report_rows.append({
                "Energy": energy,
                "SSD": ssd,
                "ScanType": scan,
                "Device": dev_col,
                "Status": "Found" if found else "MISSING",
                "Where": where
            })
    return pd.DataFrame(report_rows)



def run_audit(root, req_path, out_path, log_cb=lambda m: None):
    import time
    _t0 = time.time()
    log_cb("Loading requirements…")

    # Read requirements file (SCANS.xlsx)
    req_xls = pd.ExcelFile(req_path)
    req_tabs = []
    for sh in req_xls.sheet_names:
        # Parse Energy + SSD from sheet name
        energy, ssd = parse_energy_ssd_from_sheet(sh)
        if not energy or not ssd:
            log_cb(f"WARNING: Could not parse Energy/SSD from sheet '{sh}'")
        df = req_xls.parse(sh)
        df.columns = [c.strip() for c in df.columns]

        # Rename columns to standard names
        colmap = {
            "Energy": "Energy",
            "Field Size": "FieldSize",
            "Scan Type": "ScanType",
            "Depth": "Depth",
            "SSD": "SSD",
        }
        for c in list(colmap):
            if c not in df.columns:
                # Create if missing
                df[c] = None
        df = df.rename(columns=colmap)

        # Fill Energy/SSD from sheet name
        df["Energy"] = energy
        df["SSD"] = ssd

        # Normalize values
        df["Energy"]    = df["Energy"].astype(str).str.upper().str.strip()
        df["ScanType"]  = df["ScanType"].map(norm_scan_type)
        df["SSD"]       = pd.to_numeric(df["SSD"], errors="coerce")
        df["FieldSize"] = pd.to_numeric(df["FieldSize"], errors="coerce")
        df["Depth"]     = pd.to_numeric(df["Depth"], errors="coerce")

        req_tabs.append(df)

    req_df = pd.concat(req_tabs, ignore_index=True)

    # Helper: required scan types for an E/SSD
    FALLBACK_SCAN_SET = ["X", "Y", "XY", "YX", "Z"]
    def needed_scans_for(e, ssd):
        sub = req_df[(req_df["Energy"] == e) & (req_df["SSD"] == ssd)]
        scans = sorted(set(v for v in sub["ScanType"] if isinstance(v, str) and v))
        return scans if scans else FALLBACK_SCAN_SET[:]

    # Helper: required FS/Depths for E/SSD/ST
    def required_sets(e, ssd, st):
        sub = req_df[(req_df["Energy"] == e) & (req_df["SSD"] == ssd) & (req_df["ScanType"] == st)]
        rfs = sorted(set(sub["FieldSize"].dropna()))
        rdp = sorted(set(sub["Depth"].dropna())) if st != "Z" else []
        return rfs, rdp

    # Crawl & build inventory
    root = Path(root)
    log_cb(f"Scanning for Excel files under: {root}")
    excel_files = [p for p in root.rglob("*.xlsx")
                   if not any(ex in p.name for ex in ["~$", Path(req_path).name])]
    if not excel_files:
        raise RuntimeError("No .xlsx files found under selected root.")

    total = len(excel_files)
    inv_parts = [None] * total
    cache = load_cache()

    # Separate files that need re-reading from those cached
    to_read = {}
    for i, p in enumerate(excel_files):
        mtime = p.stat().st_mtime
        key = str(p)
        if key in cache and cache[key][0] == mtime:
            inv_parts[i] = cache[key][1]
            log_cb(f"[cached] {p.name}")
        else:
            to_read[i] = p

    # Read only new/modified files in parallel
    if to_read:
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = {executor.submit(inventory_for_excel, p): (i, p) for i, p in to_read.items()}
            for future in as_completed(futures):
                i, p = futures[future]
                result = future.result()
                inv_parts[i] = result
                cache[str(p)] = (p.stat().st_mtime, result)
                log_cb(f"[{sum(x is not None for x in inv_parts)}/{total}] Indexed {p.name} …")
        save_cache(cache)

    inventory = (pd.concat(inv_parts, ignore_index=True)
                 if inv_parts else pd.DataFrame(
                     columns=["File","Device","Energy","SSD","ScanType","FieldSizes","Depths","NumScans"]))

    inventory["Energy"]   = inventory["Energy"].astype(str).str.upper().str.strip()
    inventory["SSD"]      = pd.to_numeric(inventory["SSD"], errors="coerce")
    inventory["ScanType"] = inventory["ScanType"].map(norm_scan_type)

    # Devices present per (E,SSD)
    dev_by_group = (
        inventory.groupby(["Energy", "SSD"])["Device"]
                 .apply(lambda s: sorted(set(v for v in s if isinstance(v, str) and v)))
                 .to_dict()
    )

    # Existing keys (E,SSD,ST,Device)
    inv_keys = set(zip(inventory["Energy"], inventory["SSD"], inventory["ScanType"], inventory["Device"]))

    # Add rows for required but missing scan types
    missing_rows = []
    for (e, ssd), devices in dev_by_group.items():
        for st in needed_scans_for(e, ssd):
            for dev in devices:
                if (e, ssd, st, dev) not in inv_keys:
                    missing_rows.append({
                        "File": "",
                        "Device": dev,
                        "Energy": e,
                        "SSD": ssd,
                        "ScanType": st,
                        "FieldSizes": [],
                        "Depths": [],
                        "NumScans": 0,
                    })
    if missing_rows:
        missing_df = pd.DataFrame(missing_rows)
        # Keep only truly missing combinations
        existing_keys = set(zip(inventory["Energy"], inventory["SSD"],
                                inventory["ScanType"], inventory["Device"]))
        missing_df = missing_df[~missing_df.apply(
            lambda r: (r["Energy"], r["SSD"], r["ScanType"], r["Device"]) in existing_keys,
            axis=1)]
        inventory = pd.concat([inventory, missing_df], ignore_index=True)


    # Compute Missing FS/Depths
    missing_fs, missing_dp, complete = [], [], []
    for _, r in inventory.iterrows():
        energy_val = str(r["Energy"]).upper().strip()
        ssd_val    = pd.to_numeric(r["SSD"], errors="coerce")
        st         = norm_scan_type(r["ScanType"])
    
        req_fs, req_dp = required_sets(energy_val, ssd_val, st)
        found_fs = sorted(set(r["FieldSizes"])) if isinstance(r["FieldSizes"], (list, set, tuple)) else []
        found_dp = sorted(set(r["Depths"]))     if isinstance(r["Depths"],   (list, set, tuple)) else []
    
        mf = [x for x in req_fs if x not in found_fs]
        md = [x for x in req_dp if x not in found_dp]
    
        required_here = st in needed_scans_for(energy_val, ssd_val)
        if required_here and int(r.get("NumScans", 0)) == 0:
            mf = req_fs[:]
            md = req_dp[:]
            comp = "No"
        else:
            comp = "Yes" if (len(mf) == 0 and len(md) == 0) else "No"

        missing_fs.append(mf)
        missing_dp.append(md)
        complete.append(comp)

    inventory["MissingFieldSizes"] = missing_fs
    inventory["MissingDepths"]     = missing_dp
    inventory["Complete"]          = complete

    # Sort & write
    inventory = inventory.sort_values(
        by=["Device", "Energy", "SSD", "ScanType", "File"],
        kind="mergesort"
    ).reset_index(drop=True)

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    out_path = Path(out_path)
    out_path = out_path.parent / f"{out_path.stem}_{timestamp}{out_path.suffix}"
    log_cb(f"Writing: {out_path}")
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as xw:
        inventory.to_excel(xw, index=False, sheet_name="Inventory")

    elapsed = time.time() - _t0
    log_cb(f"Done.  Total time: {elapsed:.1f} s")
    return out_path






# ---------- GUI ----------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Scan Coverage Checker")
        self.geometry("720x460")

        pad = {"padx": 8, "pady": 6}

        ttk.Label(self, text="Root folder to search:").grid(row=0, column=0, sticky="w", **pad)
        self.root_var = tk.StringVar(self)
        ttk.Entry(self, textvariable=self.root_var, width=70).grid(row=0, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_root).grid(row=0, column=2, **pad)

        ttk.Label(self, text="Requirements (SCANS.xlsx):").grid(row=1, column=0, sticky="w", **pad)
        self.req_var = tk.StringVar(self,value=r"C:/Users/nknutson/OneDrive - Washington University in St. Louis/NGDS QA Consortium/Wash U - NGDS QA Folder/SN 21/Measurement/SCANS.xlsx")
        ttk.Entry(self, textvariable=self.req_var, width=70).grid(row=1, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_req).grid(row=1, column=2, **pad)

        ttk.Label(self, text="Output report:").grid(row=2, column=0, sticky="w", **pad)
        self.out_var = tk.StringVar(self,value=str(Path.home() / "ScanBatchCoverage.xlsx"))
        ttk.Entry(self, textvariable=self.out_var, width=70).grid(row=2, column=1, **pad)
        ttk.Button(self, text="Choose…", command=self.pick_out).grid(row=2, column=2, **pad)

        self.run_btn = ttk.Button(self, text="Run Audit", command=self.on_run)
        self.run_btn.grid(row=3, column=1, sticky="e", **pad)

        self.log = tk.Text(self, height=14, width=96, state="disabled")
        self.log.grid(row=4, column=0, columnspan=3, sticky="nsew", **pad)

        self.grid_rowconfigure(4, weight=1)
        self.grid_columnconfigure(1, weight=1)

        self.pb = ttk.Progressbar(self, mode="indeterminate")
        self.pb.grid(row=5, column=0, columnspan=3, sticky="ew", **pad)

    def pick_root(self):
        d = filedialog.askdirectory(title="Select root folder")
        if d:
            self.root_var.set(d)

    def pick_req(self):
        f = filedialog.askopenfilename(title="Select SCANS.xlsx", filetypes=[("Excel files","*.xlsx")])
        if f:
            self.req_var.set(f)

    def pick_out(self):
        f = filedialog.asksaveasfilename(title="Save output Excel as", defaultextension=".xlsx",
                                         filetypes=[("Excel files","*.xlsx")])
        if f:
            self.out_var.set(f)

    def append_log(self, msg):
        self.log.configure(state="normal")
        self.log.insert("end", msg + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def open_file(self, path):
        try:
            if os.name == "nt":
                os.startfile(path)
            elif sys.platform == "darwin":
                os.system(f'open "{path}"')
            else:
                os.system(f'xdg-open "{path}"')
        except Exception:
            pass

    def on_run(self):
        root = self.root_var.get().strip()
        req  = self.req_var.get().strip()
        outp = self.out_var.get().strip()

        if not root or not os.path.isdir(root):
            messagebox.showerror("Missing root folder", "Please choose a valid root folder.")
            return
        if not req or not os.path.isfile(req):
            messagebox.showerror("Missing SCANS.xlsx", "Please choose your requirements workbook (SCANS.xlsx).")
            return
        if not outp:
            messagebox.showerror("Missing output path", "Please choose a path to save the report.")
            return

        # All UI interactions stay on the main thread
        self.run_btn.config(state="disabled")
        self.pb.start(15)
        self.append_log("Starting audit…")

        def worker():
            try:
                # NO direct Tk calls in this thread; use after() for UI updates
                saved = run_audit(root, req, outp,
                          log_cb=lambda m: self.after(0, lambda: self.append_log(m)))
                self.after(0, lambda p=saved: self.append_log(f"Report written to: {p}"))
                self.after(0, lambda p=saved: self.open_file(str(p)))
            except Exception as exc:
                err_msg = str(exc)  # capture now; 'exc' will be cleared after except
                self.after(0, lambda msg=err_msg: self.append_log(f"ERROR: {msg}"))
                self.after(0, lambda msg=err_msg: messagebox.showerror("Error", msg))
            finally:
                self.after(0, self.pb.stop)
                self.after(0, lambda: self.run_btn.config(state="normal"))


        threading.Thread(target=worker, daemon=True).start()


if __name__ == "__main__":
    App().mainloop()

