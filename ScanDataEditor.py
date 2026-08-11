"""
Scan Data Editor
----------------
Load a machine sheet from a beam data Excel file, view every individual
scan, identify duplicates, and delete unwanted scans in place.

A "scan" is uniquely identified by (Energy, SSD, FS, Axis, Depth).
Within that group, different Detector values OR repeated Pos values are
each treated as a separate "copy" and flagged as duplicates.

Workflow:
  1. Pick a source file → sheet dropdown populates automatically.
  2. Pick a sheet (machine tab, e.g. SN1003) and click Load.
  3. Every scan appears as a row. Click any column header to sort.
     Use the filter dropdowns to narrow the view.
  4. Check (☑) scans you want to DELETE — duplicate extras are
     pre-checked automatically.
  5. Click "Delete marked" → original is backed up as _orig.xlsx and
     saved in place with marked scans removed (no blank rows).
"""

import hashlib
import threading
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


CHECKED   = '☑'
UNCHECKED = '☐'

TAG_UNIQUE    = 'unique'
TAG_DUP_IDENT = 'dup_ident'
TAG_DUP_CONFL = 'dup_confl'

# Group key: a "scan position" — Detector is handled inside the group
PROFILE_KEY_COLS = ['Energy', 'SSD', 'FS', 'Axis', 'Depth']
DEPTH_KEY_COLS   = ['Energy', 'SSD', 'FS']

# (treeview column id, scan dict key, numeric sort?)
TREE_COLS = [
    ('check',    None,       False),
    ('#',        None,       True),
    ('Energy',   'energy',   False),
    ('SSD',      'ssd',      True),
    ('FS',       'fs',       True),
    ('Axis',     'axis',     False),
    ('Depth',    'depth',    True),
    ('Detector', 'detector', False),
    ('Points',   'points',   True),
    ('Status',   'status',   False),
]

# Filter bar: (label, scan dict key)
FILTER_COLS = [
    ('Energy',   'energy'),
    ('SSD',      'ssd'),
    ('FS',       'fs'),
    ('Axis',     'axis'),
    ('Detector', 'detector'),
    ('Status',   'status'),
]


# ── Core logic ────────────────────────────────────────────────────────────────

def _fingerprint(sub_df):
    arr = sub_df[['Pos', 'Dose']].apply(pd.to_numeric, errors='coerce').dropna().to_numpy()
    if len(arr) == 0:
        return ''
    arr = np.round(arr[np.argsort(arr[:, 0])], 6)
    return hashlib.md5(arr.tobytes()).hexdigest()


def _scalar(v):
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except Exception:
        return v


def _copies_from_group(det_val, det_group):
    """Split one (key + detector) sub-group into copies using Pos occurrence.
    Returns list of (det_val, sub_df, fingerprint).
    """
    pos_num  = pd.to_numeric(det_group['Pos'], errors='coerce')
    n_total  = len(det_group)
    n_unique = int(pos_num.nunique(dropna=True))

    if n_unique == 0 or n_total == n_unique or (n_unique > 0 and n_total % n_unique != 0):
        return [(det_val, det_group, _fingerprint(det_group))]

    n_copies   = n_total // n_unique
    occurrence = pos_num.groupby(pos_num).cumcount()
    result = []
    for copy_num in range(n_copies):
        copy_df = det_group[occurrence == copy_num]
        if not copy_df.empty:
            result.append((det_val, copy_df, _fingerprint(copy_df)))
    return result


def parse_scans(df):
    """Return list of scan dicts.

    Groups by (Energy, SSD, FS, Axis, Depth).  Within each group, each
    distinct Detector — and each repeated-Pos set — becomes its own copy.
    Multiple copies → duplicate status.

    Dict keys: energy, ssd, fs, axis, depth, detector, points,
               fingerprint, row_indices, status, seg_num, n_segs, delete
    """
    has_axis     = 'Axis'     in df.columns
    has_detector = 'Detector' in df.columns
    candidates   = PROFILE_KEY_COLS if has_axis else DEPTH_KEY_COLS
    key_cols     = [c for c in candidates if c in df.columns]

    df = df.copy()
    for c in ('SSD', 'FS', 'Depth'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    if has_axis:
        df['Axis'] = df['Axis'].astype(str).str.strip()
    if has_detector:
        df['Detector'] = df['Detector'].astype(str).str.strip()

    scans = []
    for key, group in df.groupby(key_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        kd = dict(zip(key_cols, key))

        energy_val = str(kd.get('Energy', '')).strip()
        ssd_val    = _scalar(kd.get('SSD',   ''))
        fs_val     = _scalar(kd.get('FS',    ''))
        axis_val   = kd.get('Axis', '')
        depth_val  = _scalar(kd.get('Depth', '')) if 'Depth' in kd else ''

        # Sub-group by Detector (each detector = separate copy candidate)
        if has_detector:
            det_groups = [(dv, dg) for dv, dg in group.groupby('Detector', dropna=False)]
        else:
            det_groups = [('', group)]

        # Expand each detector sub-group for pos-occurrence duplicates
        all_copies = []
        for det_val, det_group in det_groups:
            all_copies.extend(_copies_from_group(str(det_val).strip(), det_group))

        n = len(all_copies)
        fps = [fp for _, _, fp in all_copies]
        unique_fps = set(fps)

        if n == 1:
            status = 'Unique'
        elif len(unique_fps) == 1:
            status = 'Dup-Identical'
        else:
            status = 'Dup-Conflicting'

        max_pts   = max(len(sub) for _, sub, _ in all_copies)
        seen_fps: set = set()
        for i, (det_val, sub_df, fp) in enumerate(all_copies):
            if status == 'Unique':
                delete = False
            elif status == 'Dup-Identical':
                delete = (i > 0)          # keep first, delete rest
            else:
                # keep largest unique; delete smaller / repeated
                delete = not (len(sub_df) == max_pts and fp not in seen_fps)
            seen_fps.add(fp)

            scans.append({
                'energy':      energy_val,
                'ssd':         ssd_val,
                'fs':          fs_val,
                'axis':        axis_val,
                'depth':       depth_val,
                'detector':    det_val,
                'points':      len(sub_df),
                'fingerprint': fp,
                'row_indices': sub_df.index.tolist(),
                'status':      status,
                'seg_num':     i + 1,
                'n_segs':      n,
                'delete':      delete,
            })

    return scans


def write_sheet_inplace(source_path, sheet_name, cleaned_df):
    """Write cleaned_df back to source_path in place, backing up original."""
    src        = Path(source_path)
    backup_dir = src.parent / 'Backups'
    backup_dir.mkdir(exist_ok=True)
    backup = backup_dir / (src.stem + '_orig.xlsx')
    tmp    = src.with_name(src.stem + '_cleaned.xlsx')

    # Write cleaned version to temp file first (safe), then close the handle
    with pd.ExcelFile(source_path, engine='openpyxl') as xl:
        sheets = {sn: (cleaned_df if sn == sheet_name else xl.parse(sn))
                  for sn in xl.sheet_names}
    with pd.ExcelWriter(str(tmp), engine='openpyxl') as writer:
        for sn, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sn, index=False)

    # Swap: original → backup, tmp → original
    if backup.exists():
        backup.unlink()
    src.rename(backup)
    tmp.rename(src)
    return backup


# ── GUI ──────────────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Scan Data Editor")
        self.geometry("1200x760")
        self.resizable(True, True)

        self._source_df   = None
        self._scans       = []
        self._visible_idx = []
        self._source_path = None
        self._sheet_name  = None
        self._sort_col    = None
        self._sort_asc    = True
        self._filter_vars: dict[str, tk.StringVar] = {}

        self._build_ui()

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        pad = dict(padx=8, pady=4)
        main = ttk.Frame(self, padding=8)
        main.pack(fill='both', expand=True)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(5, weight=1)

        # Row 0: file
        ttk.Label(main, text="Source file:").grid(row=0, column=0, sticky='w', **pad)
        self.source_var = tk.StringVar()
        ttk.Entry(main, textvariable=self.source_var, width=90).grid(
            row=0, column=1, sticky='ew', **pad)
        ttk.Button(main, text="Browse…", command=self._pick_file).grid(row=0, column=2, **pad)

        # Row 1: sheet + load
        ttk.Label(main, text="Sheet (machine):").grid(row=1, column=0, sticky='w', **pad)
        self.sheet_var = tk.StringVar()
        self.sheet_cb  = ttk.Combobox(main, textvariable=self.sheet_var,
                                       state='readonly', width=32)
        self.sheet_cb.grid(row=1, column=1, sticky='w', **pad)
        ttk.Button(main, text="Load", command=self._load).grid(row=1, column=2, **pad)

        # Row 2: action buttons
        btn = ttk.Frame(main)
        btn.grid(row=2, column=0, columnspan=3, sticky='w', **pad)
        ttk.Button(btn, text="Mark all",
                   command=self._mark_all).pack(side='left', padx=4)
        ttk.Button(btn, text="Clear marks",
                   command=self._clear_marks).pack(side='left', padx=4)
        ttk.Button(btn, text="Auto-mark duplicates",
                   command=self._auto_mark).pack(side='left', padx=12)
        self.delete_btn = ttk.Button(btn, text="Delete marked →",
                                     command=self._delete)
        self.delete_btn.pack(side='right', padx=4)

        # Row 3: filter bar
        fbar = ttk.LabelFrame(main, text="Filter")
        fbar.grid(row=3, column=0, columnspan=3, sticky='ew', padx=8, pady=(2, 0))
        for label, key in FILTER_COLS:
            ttk.Label(fbar, text=f"{label}:").pack(side='left', padx=(8, 2))
            var = tk.StringVar(value='All')
            self._filter_vars[key] = var
            cb = ttk.Combobox(fbar, textvariable=var, state='readonly', width=12)
            cb['values'] = ['All']
            cb.pack(side='left', padx=(0, 4))
            var.trace_add('write', lambda *_: self._apply_view())
            setattr(self, f'_fcb_{key}', cb)
        ttk.Button(fbar, text="Clear filters",
                   command=self._clear_filters).pack(side='right', padx=8)

        # Row 4: summary
        self.summary_lbl = ttk.Label(main, text="")
        self.summary_lbl.grid(row=4, column=0, columnspan=3, sticky='w', padx=8, pady=(2, 0))

        # Row 5: treeview
        col_ids = [c for c, _, _ in TREE_COLS]
        tf = ttk.Frame(main)
        tf.grid(row=5, column=0, columnspan=3, sticky='nsew', padx=8, pady=4)
        tf.rowconfigure(0, weight=1)
        tf.columnconfigure(0, weight=1)

        self.tree = ttk.Treeview(tf, columns=col_ids, show='headings',
                                 height=20, selectmode='none')
        widths  = {'check': 32, '#': 45, 'Energy': 65, 'SSD': 50, 'FS': 55,
                   'Axis': 50, 'Depth': 65, 'Detector': 120, 'Points': 60, 'Status': 200}
        centers = {'check', '#', 'SSD', 'FS', 'Depth', 'Points', 'Axis'}
        for col_id, scan_key, _ in TREE_COLS:
            hkw = {'text': '' if col_id == 'check' else col_id}
            if scan_key:
                hkw['command'] = lambda k=scan_key: self._on_header(k)
            self.tree.heading(col_id, **hkw)
            self.tree.column(col_id, width=widths.get(col_id, 80),
                             anchor='center' if col_id in centers else 'w',
                             stretch=False)

        vsb = ttk.Scrollbar(tf, orient='vertical',   command=self.tree.yview)
        hsb = ttk.Scrollbar(tf, orient='horizontal', command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')
        self.tree.bind('<Button-1>', self._on_click)

        self.tree.tag_configure(TAG_UNIQUE,    background='')
        self.tree.tag_configure(TAG_DUP_IDENT, background='#FFF3CD')
        self.tree.tag_configure(TAG_DUP_CONFL, background='#FDDDE6')

        # Row 6: log
        lf = ttk.Frame(main)
        lf.grid(row=6, column=0, columnspan=3, sticky='ew', padx=8, pady=(4, 0))
        lf.columnconfigure(0, weight=1)
        self.log = tk.Text(lf, height=4, state='disabled',
                           font=('Courier', 9), wrap='none')
        lsb = ttk.Scrollbar(lf, orient='vertical', command=self.log.yview)
        self.log.configure(yscrollcommand=lsb.set)
        self.log.grid(row=0, column=0, sticky='ew')
        lsb.grid(row=0, column=1, sticky='ns')

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _log(self, msg):
        self.log.configure(state='normal')
        self.log.insert('end', msg + '\n')
        self.log.see('end')
        self.log.configure(state='disabled')

    def _update_summary(self):
        n_unique = sum(1 for s in self._scans if s['status'] == 'Unique')
        n_ident  = sum(1 for s in self._scans if s['status'] == 'Dup-Identical')
        n_confl  = sum(1 for s in self._scans if s['status'] == 'Dup-Conflicting')
        n_delete = sum(1 for s in self._scans if s['delete'])
        shown    = len(self._visible_idx)
        total    = len(self._scans)
        filtered = f"  (showing {shown} of {total})" if shown != total else ""
        self.summary_lbl.config(
            text=f"{total} scans total{filtered}  |  "
                 f"{n_unique} unique  |  "
                 f"{n_ident} identical dup  |  "
                 f"{n_confl} conflicting dup  |  "
                 f"{n_delete} marked for deletion"
        )

    def _check_display(self, scan):
        return CHECKED if scan['delete'] else UNCHECKED

    def _row_values(self, scan, num):
        status = scan['status']
        if scan['n_segs'] > 1:
            status += f"  (copy {scan['seg_num']} of {scan['n_segs']})"
        return (
            self._check_display(scan),
            num,
            scan['energy'], scan['ssd'], scan['fs'], scan['axis'],
            scan['depth'] if scan['depth'] != '' else '',
            scan['detector'], scan['points'], status,
        )

    # ── File / load ───────────────────────────────────────────────────────────

    def _pick_file(self):
        p = filedialog.askopenfilename(
            title='Select beam data file',
            filetypes=[('Excel files', '*.xlsx'), ('All files', '*.*')])
        if not p:
            return
        self.source_var.set(p)
        try:
            xl = pd.ExcelFile(p)
            self.sheet_cb['values'] = xl.sheet_names
            if xl.sheet_names:
                self.sheet_var.set(xl.sheet_names[0])
        except Exception as e:
            messagebox.showerror('Error', f'Cannot open file:\n{e}')

    def _load(self):
        path  = self.source_var.get().strip()
        sheet = self.sheet_var.get().strip()
        if not path or not sheet:
            messagebox.showwarning('Missing input', 'Choose a file and sheet first.')
            return
        self._log(f"Loading  {Path(path).name}  →  {sheet}…")
        try:
            df = pd.read_excel(path, sheet_name=sheet)
            self._source_df   = df
            self._source_path = path
            self._sheet_name  = sheet
            self._scans = parse_scans(df)
            self._sort_col = None
            self._sort_asc = True
            self._refresh_filter_values()
            self._clear_filters()
            self._update_summary()
            self._log(f"Loaded {len(self._scans)} scans  ({len(df)} rows total).")
        except Exception as e:
            self._log(f"ERROR: {e}")
            messagebox.showerror('Error', str(e))

    # ── Filter ────────────────────────────────────────────────────────────────

    def _refresh_filter_values(self):
        for _, key in FILTER_COLS:
            cb = getattr(self, f'_fcb_{key}', None)
            if cb is None:
                continue
            vals = sorted({str(s[key]) for s in self._scans if str(s[key]) != ''},
                          key=lambda x: (float(x),) if _is_numeric(x) else (float('inf'), x))
            cb['values'] = ['All'] + vals

    def _clear_filters(self):
        for var in self._filter_vars.values():
            var.set('All')

    def _matches_filters(self, scan):
        for key, var in self._filter_vars.items():
            val = var.get()
            if val == 'All':
                continue
            if str(scan.get(key, '')) != val:
                return False
        return True

    # ── Sort / view ───────────────────────────────────────────────────────────

    def _on_header(self, scan_key):
        if self._sort_col == scan_key:
            self._sort_asc = not self._sort_asc
        else:
            self._sort_col = scan_key
            self._sort_asc = True
        self._update_header_labels()
        self._apply_view()

    def _update_header_labels(self):
        for col_id, scan_key, _ in TREE_COLS:
            if col_id == 'check' or scan_key is None:
                continue
            arrow = (' ↑' if self._sort_asc else ' ↓') if scan_key == self._sort_col else ''
            self.tree.heading(col_id, text=col_id + arrow)

    def _sort_key_fn(self, scan_key):
        numeric = any(scan_key == k for _, k, num in TREE_COLS if num)
        if numeric:
            def key(idx):
                try:
                    return float(self._scans[idx].get(scan_key, '') or 0)
                except Exception:
                    return float('inf')
        else:
            def key(idx):
                return str(self._scans[idx].get(scan_key, '') or '')
        return key

    def _apply_view(self):
        visible = [i for i, s in enumerate(self._scans) if self._matches_filters(s)]
        if self._sort_col:
            visible.sort(key=self._sort_key_fn(self._sort_col), reverse=not self._sort_asc)
        self._visible_idx = visible
        self._repopulate_tree()
        self._update_summary()

    def _repopulate_tree(self):
        self.tree.delete(*self.tree.get_children())
        tag_map = {
            'Unique':          TAG_UNIQUE,
            'Dup-Identical':   TAG_DUP_IDENT,
            'Dup-Conflicting': TAG_DUP_CONFL,
        }
        for num, idx in enumerate(self._visible_idx, start=1):
            scan = self._scans[idx]
            self.tree.insert('', 'end', tags=(tag_map.get(scan['status'], TAG_UNIQUE),),
                             values=self._row_values(scan, num))

    # ── Click toggle ──────────────────────────────────────────────────────────

    def _on_click(self, event):
        if self.tree.identify_region(event.x, event.y) != 'cell':
            return
        if self.tree.identify_column(event.x) != '#1':
            return
        iid = self.tree.identify_row(event.y)
        if not iid:
            return
        pos  = list(self.tree.get_children()).index(iid)
        scan = self._scans[self._visible_idx[pos]]
        scan['delete'] = not scan['delete']
        self.tree.set(iid, 'check', self._check_display(scan))
        self._update_summary()

    # ── Mark buttons (operate on visible rows) ────────────────────────────────

    def _mark_all(self):
        for pos, iid in enumerate(self.tree.get_children()):
            scan = self._scans[self._visible_idx[pos]]
            scan['delete'] = True
            self.tree.set(iid, 'check', CHECKED)
        self._update_summary()

    def _clear_marks(self):
        for pos, iid in enumerate(self.tree.get_children()):
            scan = self._scans[self._visible_idx[pos]]
            scan['delete'] = False
            self.tree.set(iid, 'check', UNCHECKED)
        self._update_summary()

    def _auto_mark(self):
        """Re-apply default duplicate logic across ALL scans."""
        groups = defaultdict(list)
        for i, scan in enumerate(self._scans):
            groups[(scan['energy'], scan['ssd'], scan['fs'],
                    scan['axis'], scan['depth'])].append(i)

        for indices in groups.values():
            group = [self._scans[i] for i in indices]
            if len(group) == 1:
                group[0]['delete'] = False
                continue
            fps      = [s['fingerprint'] for s in group]
            max_pts  = max(s['points'] for s in group)
            seen_fps: set = set()
            all_identical = len(set(fps)) == 1
            for scan in group:
                if all_identical:
                    scan['delete'] = (scan['seg_num'] > 1)
                else:
                    keep = (scan['points'] == max_pts and
                            scan['fingerprint'] not in seen_fps)
                    scan['delete'] = not keep
                seen_fps.add(scan['fingerprint'])

        for pos, iid in enumerate(self.tree.get_children()):
            scan = self._scans[self._visible_idx[pos]]
            self.tree.set(iid, 'check', self._check_display(scan))
        self._update_summary()

    # ── Delete ────────────────────────────────────────────────────────────────

    def _delete(self):
        if not self._scans:
            messagebox.showinfo('Nothing loaded', 'Load a sheet first.')
            return

        to_delete = [s for s in self._scans if s['delete']]
        if not to_delete:
            messagebox.showinfo('Nothing marked', 'Check (☑) the scans you want to delete.')
            return

        to_keep       = [s for s in self._scans if not s['delete']]
        n_rows_before = sum(len(s['row_indices']) for s in self._scans)
        n_rows_after  = sum(len(s['row_indices']) for s in to_keep)
        src           = Path(self._source_path)

        if not messagebox.askyesno('Confirm delete',
                f"Delete {len(to_delete)} scan(s)  "
                f"({n_rows_before - n_rows_after} rows).\n"
                f"Keep {len(to_keep)} scan(s)  ({n_rows_after} rows).\n\n"
                f"Original backed up as:  {src.stem}_orig.xlsx\n\n"
                f"Continue?"):
            return

        self._log(f"Deleting {len(to_delete)} scans from {src.name}…")
        self.delete_btn.config(state='disabled')

        source_df   = self._source_df
        source_path = self._source_path
        sheet_name  = self._sheet_name

        def worker():
            try:
                kept_idx = []
                for s in to_keep:
                    kept_idx.extend(s['row_indices'])
                cleaned_df = source_df.loc[sorted(kept_idx)].reset_index(drop=True)
                backup = write_sheet_inplace(source_path, sheet_name, cleaned_df)
                msg = (f"Deleted {len(to_delete)} scan(s).\n"
                       f"Backup: Backups/{backup.name}\n"
                       f"Saved:  {Path(source_path).name}")
                self.after(0, lambda m=msg: self._log(m))
                self.after(0, lambda m=msg: messagebox.showinfo('Done', m))
                self.after(0, self._load)
            except Exception as e:
                err = str(e)
                self.after(0, lambda m=err: self._log(f"ERROR: {m}"))
                self.after(0, lambda m=err: messagebox.showerror('Error', m))
            finally:
                self.after(0, lambda: self.delete_btn.config(state='normal'))

        threading.Thread(target=worker, daemon=True).start()


def _is_numeric(s):
    try:
        float(s)
        return True
    except Exception:
        return False


if __name__ == '__main__':
    App().mainloop()
