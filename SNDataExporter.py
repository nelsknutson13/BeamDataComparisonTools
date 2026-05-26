"""
SN Data Exporter
----------------
Scan a combined-data root folder and export every scan belonging to one
machine (SN) into a single workbook.

Source layout (per the NGDS Combined Consortium Data tree):
    <root>/<...>/{Energy}_{SSD}SSD_PDDData.xlsx       (sheet per SN)
    <root>/<...>/{Energy}_{SSD}SSD_ProfileData.xlsx   (sheet per SN)
    <root>/<...>/OF_Data/NGDSOutputFactorData.xlsx    (long format, 'SN' column)

Output: <SN>_AllData.xlsx
    - one tab per (Energy, SSD)  →  combined PDD + Profile rows for that SN
    - one 'OutputFactors' tab    →  all of that SN's OF rows
"""

import os
import re
import glob
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import pandas as pd


DEFAULT_ROOT = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\Combined Consortium Data"

_ENERGY_RE = re.compile(r'(?<![A-Za-z0-9])(6X|6FFF|8FFF|10X|10FFF|15X)(?![A-Za-z0-9])', re.I)
_SSD_RE    = re.compile(r'(?<!\d)(\d{2,3})\s*SSD', re.I)

OF_FILENAME_HINT = 'OutputFactor'   # substring identifying the OF workbook


def _parse_energy_ssd(name):
    e = _ENERGY_RE.search(str(name))
    s = _SSD_RE.search(str(name))
    energy = e.group(1).upper() if e else None
    ssd    = int(s.group(1)) if s else None
    return energy, ssd


def _scan_root(root):
    """Walk the root, returning:
       pdd_profile_files: list of (path, energy, ssd, kind)  kind in {'PDD','Profile'}
       of_files:          list of paths to OF workbooks
    """
    pp_files = []
    of_files = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if not f.lower().endswith('.xlsx') or f.startswith('~$'):
                continue
            full = os.path.join(dirpath, f)
            if OF_FILENAME_HINT.lower() in f.lower():
                of_files.append(full)
                continue
            energy, ssd = _parse_energy_ssd(f)
            if energy is None or ssd is None:
                continue
            if 'pdd' in f.lower():
                pp_files.append((full, energy, ssd, 'PDD'))
            elif 'profile' in f.lower():
                pp_files.append((full, energy, ssd, 'Profile'))
    return pp_files, of_files


def _discover_sns(pp_files, of_files):
    """Return a sorted list of available SN identifiers across all sources."""
    sns = set()
    for path, _e, _s, _k in pp_files:
        try:
            sns.update(str(n).strip() for n in pd.ExcelFile(path).sheet_names)
        except Exception:
            pass
    for path in of_files:
        try:
            df = pd.read_excel(path, sheet_name='OutputFactors')
            if 'SN' in df.columns:
                sns.update(str(n).strip() for n in df['SN'].dropna().unique())
        except Exception:
            pass
    # Keep only SN-looking entries; sort by trailing number when possible.
    sn_list = [s for s in sns if re.search(r'\bSN', s, re.I) or re.fullmatch(r'SN\s*\d+', s, re.I)]
    if not sn_list:
        sn_list = sorted(sns)   # fall back to everything
    def _key(s):
        m = re.search(r'(\d+)', s)
        return (int(m.group(1)) if m else 0, s)
    return sorted(set(sn_list), key=_key)


def export_sn(root, sn, status_cb):
    pp_files, of_files = _scan_root(root)
    if not pp_files and not of_files:
        raise RuntimeError("No recognized data files found under the root folder.")

    # Group PDD + Profile data by (energy, ssd)
    combo_data  = {}   # (energy, ssd) -> list of DataFrames
    combo_kinds = {}   # (energy, ssd) -> {'PDD': rows, 'Profile': rows}
    for path, energy, ssd, kind in pp_files:
        try:
            xl = pd.ExcelFile(path)
        except Exception as e:
            status_cb(f"  skip (open error): {os.path.basename(path)} — {e}")
            continue
        if sn not in xl.sheet_names:
            continue
        try:
            df = pd.read_excel(xl, sheet_name=sn)
        except Exception as e:
            status_cb(f"  skip (read error): {os.path.basename(path)} — {e}")
            continue
        if df.empty:
            continue
        combo_data.setdefault((energy, ssd), []).append(df)
        ck = combo_kinds.setdefault((energy, ssd), {'PDD': 0, 'Profile': 0})
        ck[kind] += len(df)
        status_cb(f"  {energy} {ssd}SSD {kind}: {len(df)} rows")

    # Output factors for this SN
    of_frames = []
    for path in of_files:
        try:
            df = pd.read_excel(path, sheet_name='OutputFactors')
        except Exception:
            continue
        if 'SN' in df.columns:
            sub = df[df['SN'].astype(str).str.strip() == sn]
            if not sub.empty:
                of_frames.append(sub)
                status_cb(f"  OutputFactors: {len(sub)} rows")

    if not combo_data and not of_frames:
        raise RuntimeError(f"No data found for '{sn}'.")

    # Order tabs by preferred energy then SSD
    pref = ['6X', '10X', '15X', '6FFF', '8FFF', '10FFF']
    def _combo_key(k):
        e, s = k
        return (pref.index(e) if e in pref else 99, e, s)
    ordered = sorted(combo_data.keys(), key=_combo_key)

    out_name = f"{sn.replace(' ', '')}_AllData.xlsx"
    out_path = os.path.join(root, out_name)

    # Pre-merge each combo so we can both summarize and write it.
    merged_by_combo = {(e, s): pd.concat(combo_data[(e, s)], ignore_index=True)
                       for (e, s) in ordered}
    of_all = pd.concat(of_frames, ignore_index=True) if of_frames else None

    # ── Build summary rows ──
    def _uniq_str(df, col):
        if col not in df.columns:
            return ''
        vals = sorted(df[col].dropna().unique(), key=lambda v: (str(type(v)), v))
        return ', '.join(str(v) for v in vals)

    summary_rows = []
    for (energy, ssd) in ordered:
        merged = merged_by_combo[(energy, ssd)]
        ck = combo_kinds.get((energy, ssd), {'PDD': 0, 'Profile': 0})
        summary_rows.append({
            'Tab':          f"{energy}_{ssd}SSD"[:31],
            'Energy':       energy,
            'SSD':          ssd,
            'PDD rows':     ck.get('PDD', 0),
            'Profile rows': ck.get('Profile', 0),
            'Field Sizes':  _uniq_str(merged, 'FS'),
            'Axes':         _uniq_str(merged, 'Axis'),
            'Depths':       _uniq_str(merged, 'Depth'),
            'Detector(s)':  _uniq_str(merged, 'Detector'),
        })
    if of_all is not None:
        summary_rows.append({
            'Tab':          'OutputFactors',
            'Energy':       _uniq_str(of_all, 'Energy'),
            'SSD':          _uniq_str(of_all, 'SSD'),
            'PDD rows':     '',
            'Profile rows': '',
            'Field Sizes':  f"X: {_uniq_str(of_all, 'FS_X')}",
            'Axes':         '',
            'Depths':       _uniq_str(of_all, 'Depth'),
            'Detector(s)':  _uniq_str(of_all, 'Detector'),
        })
    df_summary = pd.DataFrame(summary_rows)

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        # Summary first so it's the landing tab.
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
        for (energy, ssd) in ordered:
            merged = merged_by_combo[(energy, ssd)]
            tab = f"{energy}_{ssd}SSD"[:31]   # Excel sheet-name limit
            merged.to_excel(writer, sheet_name=tab, index=False)
        if of_all is not None:
            of_all.to_excel(writer, sheet_name='OutputFactors', index=False)

    status_cb(f"\nWrote {len(ordered)} energy/SSD tab(s)"
              f"{' + OutputFactors' if of_frames else ''}")
    status_cb(f"Saved: {out_path}")
    return out_path


# ── GUI ──────────────────────────────────────────────────────────────────────

root_win = tk.Tk()
root_win.title("SN Data Exporter v0.1")

root_var = tk.StringVar(master=root_win, value=DEFAULT_ROOT)
sn_var   = tk.StringVar(master=root_win)

main = ttk.Frame(root_win, padding=10)
main.grid(row=0, column=0, sticky="nsew")


def log(msg):
    output.insert(tk.END, str(msg) + "\n")
    output.see(tk.END)
    output.update_idletasks()


def choose_root():
    d = filedialog.askdirectory(initialdir=root_var.get() or DEFAULT_ROOT)
    if d:
        root_var.set(d)


def scan():
    r = root_var.get().strip()
    if not os.path.isdir(r):
        messagebox.showerror("Error", "Pick a valid root folder.")
        return
    output.delete("1.0", tk.END)
    log("Scanning…")
    try:
        pp, of = _scan_root(r)
        sns = _discover_sns(pp, of)
    except Exception as e:
        messagebox.showerror("Error", str(e)); return
    sn_combo['values'] = sns
    if sns:
        sn_var.set(sns[0])
    log(f"Found {len(pp)} PDD/Profile file(s), {len(of)} OF file(s).")
    log(f"Available SNs: {', '.join(sns) if sns else '(none)'}")


def do_export():
    r  = root_var.get().strip()
    sn = sn_var.get().strip()
    if not os.path.isdir(r):
        messagebox.showerror("Error", "Pick a valid root folder."); return
    if not sn:
        messagebox.showerror("Error", "Pick an SN (scan first)."); return
    output.delete("1.0", tk.END)
    log(f"Exporting '{sn}'…")
    try:
        path = export_sn(r, sn, log)
    except Exception as e:
        messagebox.showerror("Error", str(e)); log(f"[ERROR] {e}"); return
    messagebox.showinfo("Done", f"Exported to:\n{path}")


# Row 0: root folder
ttk.Label(main, text="Root folder:").grid(row=0, column=0, sticky="w")
ttk.Entry(main, textvariable=root_var, width=70).grid(row=0, column=1, padx=5)
ttk.Button(main, text="Browse…", command=choose_root).grid(row=0, column=2)
ttk.Button(main, text="Scan",    command=scan).grid(row=0, column=3, padx=4)

# Row 1: SN selection
ttk.Label(main, text="SN:").grid(row=1, column=0, sticky="w", pady=(6, 0))
sn_combo = ttk.Combobox(main, textvariable=sn_var, width=20, state='readonly')
sn_combo.grid(row=1, column=1, sticky="w", padx=5, pady=(6, 0))
ttk.Button(main, text="Export", command=do_export).grid(row=1, column=2, pady=(6, 0))

# Row 2: output log
output = tk.Text(main, height=20, width=100, wrap='none', font=('Courier', 9))
output.grid(row=2, column=0, columnspan=4, sticky="w", pady=(10, 0))


root_win.mainloop()
