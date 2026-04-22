"""
Profile SSD Converter
---------------------
Geometrically scale profile positions from one SSD to another.

For a ray from the source that passes through position `Pos` at depth `d` at
SSD_source, the same ray at SSD_target passes through:
    Pos_new = Pos * (SSD_target + d) / (SSD_source + d)

Dose values are kept unchanged (relative profiles normalized to CAX are
insensitive to inverse-square since it cancels across the profile).

Data format expected: same as the IBA / PTW output — a sheet (typically
"Profile Scans") with columns FS, Axis, Depth, Pos, Dose, SSD, Energy, [Detector].
"""

import os
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import numpy as np
import pandas as pd


# Sheet names tolerated for profile data.
PROFILE_SHEET_CANDIDATES = ['Profile Scans', 'Profiles']


def _find_profile_sheet(xlsx_path):
    xl = pd.ExcelFile(xlsx_path)
    for name in PROFILE_SHEET_CANDIDATES:
        if name in xl.sheet_names:
            return name
    return None


def _profile_groups(df):
    """Return a DataFrame of unique profile identifiers (one row per profile)."""
    keys = [c for c in ('Energy', 'SSD', 'FS', 'Axis', 'Depth') if c in df.columns]
    if not keys:
        return pd.DataFrame()
    return df[keys].drop_duplicates().reset_index(drop=True)


def _convert_positions(df, target_ssd):
    """Return a new DataFrame with Pos scaled to target_ssd and SSD column updated."""
    out = df.copy()
    # new_pos = old_pos * (ssd_target + depth) / (ssd_source + depth)
    src = out['SSD'].astype(float)
    dep = out['Depth'].astype(float)
    factor = (float(target_ssd) + dep) / (src + dep)
    out['Pos'] = out['Pos'].astype(float) * factor
    out['SSD'] = float(target_ssd)
    return out


# ── GUI ──────────────────────────────────────────────────────────────────────

root = tk.Tk()
root.title("Profile SSD Converter v0.1")

source_var = tk.StringVar(master=root)
target_ssd_var = tk.StringVar(master=root, value="100")
output_var = tk.StringVar(master=root)

_df_profile   = None        # full profile DataFrame
_profile_keys = None        # DataFrame of unique profile identifiers

main = ttk.Frame(root, padding=10)
main.grid(row=0, column=0, sticky="nsew")


def _default_output_name(src_path, target_ssd):
    d = os.path.dirname(src_path)
    stem, ext = os.path.splitext(os.path.basename(src_path))
    return os.path.join(d, f"{stem}_SSD{int(round(float(target_ssd)))}{ext}")


def choose_source():
    global _df_profile, _profile_keys
    p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls"), ("All", "*.*")])
    if not p:
        return
    source_var.set(p)
    try:
        sheet = _find_profile_sheet(p)
        if sheet is None:
            messagebox.showerror("Error",
                                 "No profile sheet found (looked for "
                                 f"{PROFILE_SHEET_CANDIDATES}).")
            return
        _df_profile = pd.read_excel(p, sheet_name=sheet)
    except Exception as e:
        messagebox.showerror("Error", f"Could not read source: {e}")
        return

    _profile_keys = _profile_groups(_df_profile)
    profile_list.delete(0, tk.END)
    for _, row in _profile_keys.iterrows():
        parts = [f"{k}={row[k]}" for k in _profile_keys.columns]
        profile_list.insert(tk.END, " | ".join(parts))
    # Select all by default
    if len(_profile_keys):
        profile_list.selection_set(0, tk.END)
    # Default output path
    try:
        output_var.set(_default_output_name(p, target_ssd_var.get()))
    except Exception:
        pass
    log(f"Loaded {len(_df_profile)} rows, {len(_profile_keys)} profile(s) from '{sheet}'.")


def choose_output():
    p = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                     filetypes=[("Excel", "*.xlsx")])
    if p:
        output_var.set(p)


def _update_default_output(*_):
    src = source_var.get()
    if src and target_ssd_var.get():
        try:
            output_var.set(_default_output_name(src, target_ssd_var.get()))
        except Exception:
            pass


def select_all():
    profile_list.selection_set(0, tk.END)


def select_none():
    profile_list.selection_clear(0, tk.END)


def log(msg):
    output_text.insert(tk.END, str(msg) + "\n")
    output_text.see(tk.END)
    output_text.update_idletasks()


def convert():
    if _df_profile is None or _profile_keys is None:
        messagebox.showerror("Error", "Load a source file first.")
        return
    sel = profile_list.curselection()
    if not sel:
        messagebox.showerror("Error", "Select at least one profile.")
        return
    try:
        target_ssd = float(target_ssd_var.get())
    except ValueError:
        messagebox.showerror("Error", "Target SSD must be numeric.")
        return

    selected_keys = _profile_keys.iloc[list(sel)]
    # Build a mask for matching rows in the source DataFrame
    keys = list(selected_keys.columns)
    merged = _df_profile.merge(selected_keys, on=keys, how='inner')
    if merged.empty:
        messagebox.showerror("Error", "No rows matched the selected profiles.")
        return

    converted = _convert_positions(merged, target_ssd)

    out_path = output_var.get().strip()
    if not out_path:
        out_path = _default_output_name(source_var.get(), target_ssd)

    try:
        with pd.ExcelWriter(out_path, engine='openpyxl') as w:
            converted.to_excel(w, sheet_name='Profile Scans', index=False)
    except Exception as e:
        messagebox.showerror("Error", f"Could not write output: {e}")
        return

    log(f"Converted {len(converted)} rows from {len(selected_keys)} profile(s).")
    log(f"Target SSD: {target_ssd} cm")
    log(f"Scaling: Pos_new = Pos_old * (SSD_target + Depth) / (SSD_source + Depth)")
    log(f"Wrote: {out_path}")


# Row 0: source file
ttk.Label(main, text="Source file:").grid(row=0, column=0, sticky="w")
ttk.Entry(main, textvariable=source_var, width=70).grid(row=0, column=1, padx=5)
ttk.Button(main, text="Browse…", command=choose_source).grid(row=0, column=2)

# Row 1: target SSD
ttk.Label(main, text="Target SSD [cm]:").grid(row=1, column=0, sticky="w", pady=(6, 0))
ttk.Entry(main, textvariable=target_ssd_var, width=12).grid(row=1, column=1, sticky="w",
                                                            padx=5, pady=(6, 0))
target_ssd_var.trace_add('write', _update_default_output)

# Row 2: output file
ttk.Label(main, text="Output file:").grid(row=2, column=0, sticky="w", pady=(6, 0))
ttk.Entry(main, textvariable=output_var, width=70).grid(row=2, column=1, padx=5, pady=(6, 0))
ttk.Button(main, text="Browse…", command=choose_output).grid(row=2, column=2, pady=(6, 0))

# Row 3: profile list
ttk.Label(main, text="Profiles to convert:").grid(row=3, column=0, sticky="nw", pady=(10, 0))
profile_frame = ttk.Frame(main)
profile_frame.grid(row=3, column=1, columnspan=2, sticky="w", pady=(10, 0))
profile_list = tk.Listbox(profile_frame, selectmode=tk.EXTENDED, height=12, width=80,
                          exportselection=False, font=('Courier', 9))
profile_list.grid(row=0, column=0, sticky="nsew")
vsb = ttk.Scrollbar(profile_frame, orient="vertical", command=profile_list.yview)
vsb.grid(row=0, column=1, sticky="ns")
profile_list.configure(yscrollcommand=vsb.set)

# Row 4: select all/none + convert
btn_frame = ttk.Frame(main)
btn_frame.grid(row=4, column=1, sticky="w", pady=6)
ttk.Button(btn_frame, text="Select all",  command=select_all ).pack(side='left', padx=3)
ttk.Button(btn_frame, text="Select none", command=select_none).pack(side='left', padx=3)
ttk.Button(btn_frame, text="Convert",     command=convert    ).pack(side='left', padx=20)

# Row 5: output log
output_text = tk.Text(main, height=10, width=100, wrap='none', font=('Courier', 9))
output_text.grid(row=5, column=0, columnspan=3, sticky="w", pady=(10, 0))


root.mainloop()
