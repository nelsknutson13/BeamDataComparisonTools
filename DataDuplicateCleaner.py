# ProfileDuplicateReviewGUI.py
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import pandas as pd
import numpy as np
import hashlib

REQUIRED = ["FS", "Depth", "Axis", "Pos", "Dose"]

def fingerprint_curve(df_grp, dec=6):
    arr = df_grp[["Pos","Dose"]].to_numpy(dtype=float)
    order = np.argsort(arr[:,0])
    arr = np.round(arr[order], dec)
    return hashlib.md5(arr.tobytes()).hexdigest()

def build_candidates(df):
    """
    Detect duplicate profile candidates per (FS,Depth,Axis). We consider a 'candidate'
    as a distinct curve identified by its fingerprint. We first drop exact duplicate rows.
    Returns (meta_candidates, sheet_stats)

    meta_candidates: list of dicts:
      {'key': (FS,Depth,Axis),
       'candidates': [
           {'fp': str, 'nrows': int, 'pos_min': float, 'pos_max': float, 'dose_min': float, 'dose_max': float, 'idx': index array}
       ]}
    """
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # logic-only view for grouping
    logic = df[REQUIRED].copy()
    for c in ("FS","Depth","Pos","Dose"):
        logic[c] = pd.to_numeric(logic[c], errors="coerce")
    logic["Axis"] = logic["Axis"].astype(str).str.strip()

    essential = logic.notna().all(axis=1)
    logic = logic.loc[essential]
    base  = df.loc[essential]  # original cols aligned

    # remove exact-duplicate rows from consideration to avoid double counting
    base_nodup = base.drop_duplicates()
    logic = logic.loc[base_nodup.index]

    meta_candidates = []
    total_profiles = 0
    dup_groups = 0

    for key, g in logic.groupby(["FS","Depth","Axis"]):
        # Hash whole meta-group
        # If the sheet contains multiple scans with same meta, they are usually appended,
        # so we need to split into distinct curves. Easiest correct way:
        # group by an md5 fingerprint of the entire (Pos,Dose) array across *the rows that share that exact array*.
        # We'll form the fingerprint per distinct set of rows with identical (Pos,Dose) pairs.
        # Implementation: group by the fingerprint of the whole group's (Pos,Dose) *values repeated*:
        # A direct groupby isn't trivial; compute fp for the entire group once — but that collapses all.
        # Instead, compute fp per row-block by grouping exact (Pos,Dose) pairs sequence:
        #
        # Practical workaround:
        #   - Sort by Pos
        #   - Identify consecutive duplicates of (Pos,Dose) across the group (not helpful)
        # Better approach:
        #   - Compute a fingerprint for the entire *group*.
        #   - Also compute a fingerprint for each contiguous block of indices with monotonic Pos reset.
        #   - If duplicates exist, we expect Pos to "restart" from min again. Split on large backward steps.
        s = g.sort_values("Pos")
        # If duplicates were appended, the simple sorted view will merge them.
        # Safer approach: operate on the *original order* so we can see resets.
        g_orig = logic.loc[g.index].sort_index()
        pos = g_orig["Pos"].to_numpy()
        # Find resets: where Pos decreases significantly (start of a new candidate)
        resets = [0]
        for i in range(1, len(pos)):
            if pos[i] < pos[i-1]:  # new scan appended
                resets.append(i)
        resets.append(len(pos))

        cand_list = []
        for a, b in zip(resets[:-1], resets[1:]):
            sub_idx = g_orig.iloc[a:b].index
            sub = g_orig.loc[sub_idx]
            if sub.empty:
                continue
            fp = fingerprint_curve(sub, dec=6)
            cand_list.append({
                "fp": fp,
                "nrows": int(len(sub)),
                "pos_min": float(np.nanmin(sub["Pos"])),
                "pos_max": float(np.nanmax(sub["Pos"])),
                "dose_min": float(np.nanmin(sub["Dose"])),
                "dose_max": float(np.nanmax(sub["Dose"])),
                "idx": sub_idx.to_numpy()
            })

        # Merge identical candidates (in case identical scans appear multiple times)
        merged = {}
        for cnd in cand_list:
            if cnd["fp"] in merged:
                merged[cnd["fp"]]["nrows"] += cnd["nrows"]
                merged[cnd["fp"]]["idx"] = np.concatenate([merged[cnd["fp"]]["idx"], cnd["idx"]])
                merged[cnd["fp"]]["pos_min"] = min(merged[cnd["fp"]]["pos_min"], cnd["pos_min"])
                merged[cnd["fp"]]["pos_max"] = max(merged[cnd["fp"]]["pos_max"], cnd["pos_max"])
                merged[cnd["fp"]]["dose_min"] = min(merged[cnd["fp"]]["dose_min"], cnd["dose_min"])
                merged[cnd["fp"]]["dose_max"] = max(merged[cnd["fp"]]["dose_max"], cnd["dose_max"])
            else:
                merged[cnd["fp"]] = cnd

        candidates = list(merged.values())
        total_profiles += 1
        if len(candidates) > 1:
            dup_groups += 1

        meta_candidates.append({"key": key, "candidates": candidates})

    stats = {
        "profiles_total": total_profiles,
        "dup_groups": dup_groups,
        "rows_considered": int(len(logic))
    }
    return meta_candidates, stats, base_nodup  # base_nodup keeps original columns

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Profile Duplicate Review")
        self.geometry("1000x720")
        self.file_path: Path | None = None

        self.df: pd.DataFrame | None = None
        self.sheet_names: list[str] = []
        self.candidates = []   # from build_candidates
        self.sheet_stats = {}
        self.base_nodup: pd.DataFrame | None = None

        self.keep_choice = {}  # (FS,Depth,Axis) -> fingerprint to keep

        # Top controls
        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")
        ttk.Label(top, text="Excel file:").pack(side="left")
        self.file_entry = ttk.Entry(top, width=70)
        self.file_entry.pack(side="left", padx=6)
        ttk.Button(top, text="Browse…", command=self.pick_file).pack(side="left", padx=4)

        ttk.Label(top, text="Sheet:").pack(side="left", padx=(16,4))
        self.sheet_combo = ttk.Combobox(top, width=28, state="readonly")
        self.sheet_combo.pack(side="left")
        self.sheet_combo.bind("<<ComboboxSelected>>", self.on_sheet)

        ttk.Button(top, text="Scan for duplicates", command=self.scan_sheet).pack(side="left", padx=8)

        # Summary
        self.summary = ttk.Label(self, text="—", padding=(10,4))
        self.summary.pack(fill="x")

        # Scrollable duplicate groups
        box = ttk.Frame(self)
        box.pack(fill="both", expand=True, padx=10, pady=(0,10))

        self.canvas = tk.Canvas(box, highlightthickness=0)
        self.inner = ttk.Frame(self.canvas)
        self.vsb = ttk.Scrollbar(box, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((0,0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

        # Actions
        actions = ttk.Frame(self, padding=10)
        actions.pack(fill="x")
        ttk.Button(actions, text="Auto-select (keep largest per group)", command=self.autoselect).pack(side="left")
        ttk.Button(actions, text="Create cleaned file", command=self.write_cleaned).pack(side="left", padx=10)

        self.status = ttk.Label(self, text="", padding=(10,4))
        self.status.pack(fill="x")

    def pick_file(self):
        p = filedialog.askopenfilename(filetypes=[("Excel files","*.xlsx")])
        if not p:
            return
        self.file_path = Path(p)
        self.file_entry.delete(0, tk.END)
        self.file_entry.insert(0, str(p))
        try:
            xls = pd.ExcelFile(self.file_path, engine="openpyxl")
            self.sheet_names = xls.sheet_names
            self.sheet_combo["values"] = self.sheet_names
            if self.sheet_names:
                self.sheet_combo.set(self.sheet_names[0])
                self.on_sheet()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read Excel: {e}")

    def on_sheet(self, *_):
        if not self.file_path:
            return
        try:
            sn = self.sheet_combo.get()
            self.df = pd.read_excel(self.file_path, sheet_name=sn, engine="openpyxl")
            self.candidates = []
            self.keep_choice.clear()
            self.sheet_stats = {}
            for w in self.inner.winfo_children():
                w.destroy()
            self.summary.config(text="Sheet loaded. Click 'Scan for duplicates'.")
            self.status.config(text="")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read sheet: {e}")

    def scan_sheet(self):
        if self.df is None:
            messagebox.showwarning("No sheet", "Pick a file and sheet first.")
            return
        try:
            cands, stats, base_nodup = build_candidates(self.df)
            self.candidates = cands
            self.sheet_stats = stats
            self.base_nodup = base_nodup

            # Build UI for duplicate groups
            for w in self.inner.winfo_children():
                w.destroy()
            groups_shown = 0
            for entry in self.candidates:
                key = entry["key"]
                cands = entry["candidates"]
                # Show only groups with >1 candidate (true duplicates)
                if len(cands) <= 1:
                    # Still register choice: keep the only one
                    if len(cands) == 1:
                        self.keep_choice[key] = cands[0]["fp"]
                    continue

                groups_shown += 1
                fs, depth, axis = key
                frame = ttk.LabelFrame(self.inner, text=f"FS={fs:g} cm  Depth={depth:g} cm  Axis={axis}", padding=8)
                frame.pack(fill="x", expand=False, pady=6)

                # Radiobuttons for which fingerprint to keep
                var = tk.StringVar(value="")
                # default: largest by nrows
                default_fp = max(cands, key=lambda x: x["nrows"])["fp"]
                var.set(default_fp)
                self.keep_choice[key] = default_fp

                # header
                hdr = ttk.Frame(frame)
                hdr.pack(fill="x")
                for i, col in enumerate(["Keep", "Rows", "Pos range", "Dose range", "Fingerprint (first 10)"]):
                    ttk.Label(hdr, text=col).grid(row=0, column=i, sticky="w", padx=(4,16))

                for r, cnd in enumerate(cands, start=1):
                    row = ttk.Frame(frame)
                    row.pack(fill="x")
                    ttk.Radiobutton(row, variable=var, value=cnd["fp"],
                                    command=lambda v=var, k=key: self._set_choice(k, v.get())
                                    ).grid(row=0, column=0, sticky="w", padx=4)
                    ttk.Label(row, text=str(cnd["nrows"])).grid(row=0, column=1, sticky="w", padx=(4,16))
                    ttk.Label(row, text=f"{cnd['pos_min']:.3f} … {cnd['pos_max']:.3f}").grid(row=0, column=2, sticky="w", padx=(4,16))
                    ttk.Label(row, text=f"{cnd['dose_min']:.4f} … {cnd['dose_max']:.4f}").grid(row=0, column=3, sticky="w", padx=(4,16))
                    ttk.Label(row, text=cnd["fp"][:10]+"…").grid(row=0, column=4, sticky="w", padx=(4,16))

            self.summary.config(
                text=f"Profiles total: {self.sheet_stats['profiles_total']}   |   "
                     f"Duplicate groups: {self.sheet_stats['dup_groups']}   |   "
                     f"Rows considered: {self.sheet_stats['rows_considered']}"
            )
            if groups_shown == 0:
                self.status.config(text="No duplicate groups found. All profiles unique.")

        except Exception as e:
            messagebox.showerror("Error", f"Scan failed: {e}")

    def _set_choice(self, key, fp):
        self.keep_choice[key] = fp

    def autoselect(self):
        """Keep largest by row count in each duplicate group."""
        if not self.candidates:
            return
        for entry in self.candidates:
            key = entry["key"]
            cands = entry["candidates"]
            if not cands:
                continue
            fp = max(cands, key=lambda x: x["nrows"])["fp"]
            self.keep_choice[key] = fp
        self.status.config(text="Auto-selected largest candidate per duplicate group.")

    def write_cleaned(self):
        if self.df is None or self.base_nodup is None:
            messagebox.showwarning("Nothing to write", "Scan a sheet first.")
            return
        try:
            # Build mask of rows to keep:
            # - All rows that belong to a (FS,Depth,Axis) with a chosen fingerprint fp
            # - All rows in groups with only one candidate (choice prefilled in scan)
            logic = self.base_nodup[REQUIRED].copy()
            for c in ("FS","Depth","Pos","Dose"):
                logic[c] = pd.to_numeric(logic[c], errors="coerce")
            logic["Axis"] = logic["Axis"].astype(str).str.strip()

            keep_idx = np.array([], dtype=int)
            kept_records = []

            # Group by meta and split into contiguous candidates again
            for key, g in logic.groupby(["FS","Depth","Axis"]):
                g_orig = logic.loc[g.index].sort_index()
                pos = g_orig["Pos"].to_numpy()
                # find resets
                resets = [0]
                for i in range(1, len(pos)):
                    if pos[i] < pos[i-1]:
                        resets.append(i)
                resets.append(len(pos))

                # decide which fp to keep
                fp_keep = self.keep_choice.get(key, None)

                cand_rows = []
                for a, b in zip(resets[:-1], resets[1:]):
                    sub_idx = g_orig.iloc[a:b].index
                    sub = g_orig.loc[sub_idx]
                    if sub.empty:
                        continue
                    fp = fingerprint_curve(sub, dec=6)
                    cand_rows.append((fp, sub_idx))

                if fp_keep is None and len(cand_rows) == 1:
                    fp_keep = cand_rows[0][0]  # only candidate

                # collect rows for kept fp
                for fp, sub_idx in cand_rows:
                    if fp == fp_keep:
                        keep_idx = np.concatenate([keep_idx, sub_idx.values])
                        kept_records.append({
                            "FS": key[0], "Depth": key[1], "Axis": key[2],
                            "kept_fp": fp, "rows_kept": len(sub_idx)
                        })

            # Build cleaned DataFrame (preserve ALL columns)
            cleaned = self.base_nodup.loc[np.sort(keep_idx)].copy()

            # Summary sheet
            sum_rows = []
            kept_df = pd.DataFrame(kept_records)
            # Count duplicate groups and total profiles kept
            groups_kept = kept_df.groupby(["FS","Depth","Axis"]).size().reset_index(name="kept_candidates")
            dup_groups = int((groups_kept["kept_candidates"] > 1).sum())  # should be 0 if exactly one kept
            sum_rows.append({
                "profiles_total": self.sheet_stats.get("profiles_total", 0),
                "dup_groups_found": self.sheet_stats.get("dup_groups", 0),
                "dup_groups_after_choice": dup_groups,
                "rows_in_considered": self.sheet_stats.get("rows_considered", 0),
                "rows_out_kept": len(cleaned)
            })
            summary_df = pd.DataFrame(sum_rows)

            # Write output workbook:
            in_path = self.file_path
            out_path = in_path.with_name(in_path.stem + "_cleaned.xlsx")

            # Write all sheets pass-through except current; replace current with cleaned; include a Summary tab
            xls = pd.ExcelFile(in_path, engine="openpyxl")
            with pd.ExcelWriter(out_path, engine="xlsxwriter") as w:
                for sn in xls.sheet_names:
                    if sn == self.sheet_combo.get():
                        cleaned.to_excel(w, sheet_name=sn, index=False)
                    else:
                        pd.read_excel(in_path, sheet_name=sn, engine="openpyxl").to_excel(w, sheet_name=sn, index=False)
                summary_df.to_excel(w, sheet_name="Summary", index=False)
                kept_df.to_excel(w, sheet_name="KeptProfiles", index=False)

            messagebox.showinfo("Done", f"Cleaned file written:\n{out_path}")
            self.status.config(text=f"Wrote cleaned file with {len(cleaned)} rows kept for sheet '{self.sheet_combo.get()}'")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to write cleaned workbook:\n{e}")

if __name__ == "__main__":
    App().mainloop()
