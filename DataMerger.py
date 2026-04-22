"""
Data Merger
-----------
Take a source xlsx (from IbaDataReader or PTWdatareader) with
"Depth Scans" and "Profile Scans"/"Profiles" sheets and distribute
its data into the appropriate destination files in a combined-data folder.

Destination naming convention (already in use):
    {Energy}_{SSD}SSD_PDDData.xlsx        ← depth scans
    {Energy}_{SSD}SSD_ProfileData.xlsx    ← profile scans

The user provides a machine name (sheet name to write). If the destination
file already has that sheet, the tool asks to replace or skip.
"""

import os
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import numpy as np
import pandas as pd


DEFAULT_DEST = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\Combined Consortium Data\Processed Combined Data"

# Source sheet names tolerated across IBA/PTW outputs
DEPTH_SHEET_CANDIDATES   = ['Depth Scans']
PROFILE_SHEET_CANDIDATES = ['Profile Scans', 'Profiles']


# ── Helpers ───────────────────────────────────────────────────────────────────

def _read_first_matching_sheet(xlsx_path, candidates):
    """Read the first sheet name from `candidates` that exists in the file; None if none."""
    try:
        xl = pd.ExcelFile(xlsx_path)
    except Exception as e:
        raise RuntimeError(f"Cannot open source file: {e}")
    for name in candidates:
        if name in xl.sheet_names:
            return pd.read_excel(xlsx_path, sheet_name=name)
    return None


def _format_ssd(ssd):
    """Return SSD as an int string (e.g. 90, 100). Accepts int, float, or string."""
    try:
        return f"{int(round(float(ssd)))}"
    except Exception:
        return str(ssd)


def _dest_filename(energy, ssd, kind):
    """Return the destination filename (no path)."""
    kind_tag = "PDDData" if kind == 'depth' else "ProfileData"
    return f"{energy}_{_format_ssd(ssd)}SSD_{kind_tag}.xlsx"


def _merge_into_destination(dest_path, sheet_name, new_df, conflict_policy,
                            status_cb=print):
    """Write new_df to sheet `sheet_name` in dest_path, handling conflicts.
    conflict_policy: 'ask', 'replace_all', 'skip_all'
    Returns one of: 'written', 'replaced', 'skipped', 'created_file'
    """
    created_file = False
    if os.path.exists(dest_path):
        try:
            existing = pd.read_excel(dest_path, sheet_name=None)
        except Exception as e:
            raise RuntimeError(f"Cannot read existing destination file: {e}")
    else:
        existing = {}
        created_file = True

    if sheet_name in existing:
        if conflict_policy == 'skip_all':
            status_cb(f"    [SKIP] '{sheet_name}' already in {os.path.basename(dest_path)}")
            return 'skipped'
        if conflict_policy == 'ask':
            answer = messagebox.askyesnocancel(
                title="Sheet exists",
                message=(f"Sheet '{sheet_name}' already exists in\n"
                         f"{os.path.basename(dest_path)}\n\n"
                         f"Yes  = Replace\n"
                         f"No   = Skip\n"
                         f"Cancel = Abort entire merge")
            )
            if answer is None:
                raise KeyboardInterrupt("User aborted merge.")
            if answer is False:
                status_cb(f"    [SKIP] user chose skip for {os.path.basename(dest_path)}")
                return 'skipped'
        # replace path
        existing[sheet_name] = new_df
        with pd.ExcelWriter(dest_path, engine='openpyxl') as w:
            for name, df in existing.items():
                df.to_excel(w, sheet_name=name, index=False)
        status_cb(f"    [REPLACED] '{sheet_name}' in {os.path.basename(dest_path)}")
        return 'replaced'

    # Sheet doesn't exist yet — just add it.
    existing[sheet_name] = new_df
    with pd.ExcelWriter(dest_path, engine='openpyxl') as w:
        for name, df in existing.items():
            df.to_excel(w, sheet_name=name, index=False)
    tag = 'created_file' if created_file else 'written'
    status_cb(f"    [{'CREATED FILE' if created_file else 'ADDED'}] '{sheet_name}' "
              f"in {os.path.basename(dest_path)}")
    return tag


# ── Core merge routine ───────────────────────────────────────────────────────

def run_merge(source_path, dest_folder, machine_name, conflict_policy, status_cb):
    """Main merge logic. status_cb is a function that accepts one string."""
    if not os.path.isfile(source_path):
        raise RuntimeError("Source file not found.")
    if not os.path.isdir(dest_folder):
        raise RuntimeError("Destination folder not found.")
    machine_name = machine_name.strip()
    if not machine_name:
        raise RuntimeError("Machine name is required.")

    status_cb(f"Source:  {source_path}")
    status_cb(f"Dest:    {dest_folder}")
    status_cb(f"Machine: {machine_name}")
    status_cb("")

    df_depth   = _read_first_matching_sheet(source_path, DEPTH_SHEET_CANDIDATES)
    df_profile = _read_first_matching_sheet(source_path, PROFILE_SHEET_CANDIDATES)

    if df_depth is None and df_profile is None:
        raise RuntimeError("Source has no recognised data sheets "
                           f"({DEPTH_SHEET_CANDIDATES + PROFILE_SHEET_CANDIDATES}).")

    written_count = replaced_count = skipped_count = created_count = 0

    for kind, df_src in (('depth', df_depth), ('profile', df_profile)):
        if df_src is None or df_src.empty:
            continue
        if 'Energy' not in df_src.columns or 'SSD' not in df_src.columns:
            status_cb(f"  WARNING: {kind} sheet missing 'Energy' or 'SSD' column — skipping.")
            continue

        # Iterate unique (Energy, SSD) groups
        groups = df_src[['Energy', 'SSD']].dropna().drop_duplicates()
        for _, row in groups.iterrows():
            energy = str(row['Energy']).strip()
            ssd    = row['SSD']
            subset = df_src[(df_src['Energy'] == row['Energy']) &
                            (df_src['SSD']    == row['SSD'])].copy()
            if subset.empty:
                continue
            fname = _dest_filename(energy, ssd, kind)
            dest  = os.path.join(dest_folder, fname)
            status_cb(f"\n  {kind.upper()} — {energy} {_format_ssd(ssd)} SSD → {fname}")
            status_cb(f"    rows: {len(subset)}")
            result = _merge_into_destination(dest, machine_name, subset,
                                             conflict_policy, status_cb)
            if   result == 'written':       written_count  += 1
            elif result == 'replaced':      replaced_count += 1
            elif result == 'skipped':       skipped_count  += 1
            elif result == 'created_file':  created_count  += 1

    status_cb("")
    status_cb("=" * 60)
    status_cb(f"  Added: {written_count}   Replaced: {replaced_count}   "
              f"Skipped: {skipped_count}   Files created: {created_count}")
    status_cb("=" * 60)


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Data Merger v0.1")

source_var   = tk.StringVar(master=root)
dest_var     = tk.StringVar(master=root, value=DEFAULT_DEST)
machine_var  = tk.StringVar(master=root)
conflict_var = tk.StringVar(master=root, value='ask')   # ask | replace_all | skip_all

main = ttk.Frame(root, padding=10)
main.grid(row=0, column=0, sticky="nsew")


def choose_source():
    p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls"), ("All", "*.*")])
    if p:
        source_var.set(p)


def choose_dest():
    d = filedialog.askdirectory(initialdir=dest_var.get() or DEFAULT_DEST)
    if d:
        dest_var.set(d)


def do_run():
    output.delete("1.0", tk.END)
    def _log(msg):
        output.insert(tk.END, str(msg) + "\n")
        output.see(tk.END)
        output.update_idletasks()
    try:
        run_merge(source_var.get(), dest_var.get(), machine_var.get(),
                  conflict_var.get(), _log)
    except KeyboardInterrupt as e:
        _log(f"\n[ABORTED] {e}")
    except Exception as e:
        _log(f"\n[ERROR] {e}")
        messagebox.showerror("Error", str(e))


# Row 0: source
ttk.Label(main, text="Source file:").grid(row=0, column=0, sticky="w")
ttk.Entry(main, textvariable=source_var, width=70).grid(row=0, column=1, padx=5)
ttk.Button(main, text="Browse…", command=choose_source).grid(row=0, column=2)

# Row 1: destination folder
ttk.Label(main, text="Destination folder:").grid(row=1, column=0, sticky="w", pady=(6, 0))
ttk.Entry(main, textvariable=dest_var, width=70).grid(row=1, column=1, padx=5, pady=(6, 0))
ttk.Button(main, text="Browse…", command=choose_dest).grid(row=1, column=2, pady=(6, 0))

# Row 2: machine name
ttk.Label(main, text="Machine name (sheet):").grid(row=2, column=0, sticky="w", pady=(6, 0))
ttk.Entry(main, textvariable=machine_var, width=30).grid(row=2, column=1, padx=5, sticky="w", pady=(6, 0))

# Row 3: conflict policy
ttk.Label(main, text="If sheet exists:").grid(row=3, column=0, sticky="w", pady=(6, 0))
conflict_frame = ttk.Frame(main)
conflict_frame.grid(row=3, column=1, sticky="w", padx=5, pady=(6, 0))
ttk.Radiobutton(conflict_frame, text="Ask each time",
                variable=conflict_var, value='ask').pack(side='left', padx=2)
ttk.Radiobutton(conflict_frame, text="Replace all",
                variable=conflict_var, value='replace_all').pack(side='left', padx=8)
ttk.Radiobutton(conflict_frame, text="Skip all",
                variable=conflict_var, value='skip_all').pack(side='left', padx=2)

# Row 4: run button
ttk.Button(main, text="Run merge", command=do_run).grid(row=4, column=0, columnspan=3, pady=10)

# Row 5: output
output_frame = ttk.Frame(main)
output_frame.grid(row=5, column=0, columnspan=3, sticky="w")
output = tk.Text(output_frame, height=20, width=100, wrap='none', font=('Courier', 9))
output.grid(row=0, column=0, sticky="nsew")
vsb = ttk.Scrollbar(output_frame, orient="vertical",   command=output.yview)
hsb = ttk.Scrollbar(output_frame, orient="horizontal", command=output.xview)
output.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
vsb.grid(row=0, column=1, sticky="ns")
hsb.grid(row=1, column=0, sticky="ew")


root.mainloop()
