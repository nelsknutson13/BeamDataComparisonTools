import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, subprocess
from scipy.stats import t as t_dist

from mcc_utils import (
    available_facets, NORM_DEPTH_CM, EM_REL_STD,
    parse_mcc_path, _group_scans, _expand_by_hv_pair,
    _group_label, _get_facet, _split_common_varying,
    _align_pair, calc_pion, _dmax_for_energy,
    _u_pion_from_beta,
)

# ── Module-level state ─────────────────────────────────────────────────────────
last_df   = None
last_meta = None
curves_summary = []

last_groups    = None
last_group_keys = []
last_label_map  = {}

group_facets = ["Energy", "Detector", "SSD_cm", "FieldSize_cm"]

_fig1 = None
_fig2 = None

# ── Legend helper ──────────────────────────────────────────────────────────────
def _attach_click_legend(fig, lines, labels):
    ax = fig.axes[0]
    short = [lbl if len(lbl) <= 55 else lbl[:52] + "…" for lbl in labels]
    leg = ax.legend(lines, short, fontsize=9, loc='upper right', framealpha=0.85)
    lined = {}
    for legline, orig in zip(leg.get_lines(), lines):
        legline.set_picker(True)
        lined[legline] = orig
    def on_pick(event):
        orig = lined.get(event.artist)
        if orig is None:
            return
        vis = not orig.get_visible()
        orig.set_visible(vis)
        event.artist.set_alpha(1.0 if vis else 0.3)
        fig.canvas.draw()
    fig.canvas.mpl_connect("pick_event", on_pick)

# ── Grouping dialog ────────────────────────────────────────────────────────────
def open_grouping_dialog():
    win = tk.Toplevel(root)
    win.title("Group by…")
    win.transient(root); win.grab_set()
    ttk.Label(win, text="Select facets:").grid(row=0, column=0, sticky="w", padx=10, pady=(10, 6))
    vars_map = {}
    for i, facet in enumerate(available_facets, start=1):
        v = tk.BooleanVar(value=(facet in group_facets))
        ttk.Checkbutton(win, text=facet, variable=v).grid(row=i, column=0, sticky="w", padx=14)
        vars_map[facet] = v
    def apply_and_close():
        chosen = [f for f in available_facets if vars_map[f].get()]
        set_group_facets(chosen or ["Energy", "Detector", "SSD_cm", "FieldSize_cm"])
        win.destroy()
    btns = ttk.Frame(win); btns.grid(row=len(available_facets)+1, column=0, pady=10, padx=10, sticky="e")
    ttk.Button(btns, text="OK",     command=apply_and_close).pack(side="left", padx=6)
    ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="left")

def set_group_facets(chosen):
    global group_facets, last_groups, last_group_keys, last_label_map
    path = file_entry.get().strip()
    if not path:
        return
    scans = parse_mcc_path(path)
    has_profile = any(len(s) >= 11 and isinstance(s[10], str) and "PROFILE" in s[10].upper() for s in scans)
    if has_profile and "ProfileDepth_cm" not in chosen:
        chosen = chosen + ["ProfileDepth_cm"]
    group_facets = chosen
    last_groups = _expand_by_hv_pair(_group_scans(scans, group_facets), group_facets)
    last_group_keys = list(last_groups.keys())
    last_label_map  = {_group_label(k, group_facets): k for k in last_group_keys}
    group_listbox.delete(0, tk.END)
    for k in last_group_keys:
        group_listbox.insert(tk.END, _group_label(k, group_facets))
    group_listbox.select_set(0, tk.END)

def browse_file():
    p = filedialog.askopenfilename(filetypes=[("PTW MCC files", "*.mcc"), ("All files", "*.*")])
    if not p: return
    file_entry.delete(0, tk.END); file_entry.insert(0, p)
    set_group_facets(group_facets)

def browse_folder():
    p = filedialog.askdirectory()
    if not p: return
    file_entry.delete(0, tk.END); file_entry.insert(0, p)
    set_group_facets(group_facets)

# ── Plotting ───────────────────────────────────────────────────────────────────
def plot_pion():
    global last_df, last_meta, last_groups, last_group_keys, last_label_map
    global curves_summary, _fig1, _fig2

    filepath = file_entry.get()
    if not filepath:
        print("No file selected"); return

    try:
        em_std = float(em_rel_std_var.get()) / 100.0
    except (ValueError, AttributeError):
        em_std = EM_REL_STD

    if last_groups is None:
        set_group_facets(group_facets)
    groups = last_groups

    sel_idxs = group_listbox.curselection()
    if not sel_idxs:
        chosen_keys = set(last_group_keys)
    else:
        labels = [group_listbox.get(i) for i in sel_idxs]
        chosen_keys = {last_label_map[lbl] for lbl in labels if lbl in last_label_map}
    if not chosen_keys:
        print("No groups selected."); return

    print(f"Plotting {len(chosen_keys)} selected groups out of {len(groups)}.")

    chosen_key_list = [k for k in last_group_keys if k in chosen_keys]
    common_str, varying_labels = _split_common_varying(chosen_key_list, group_facets)

    hv_faceted = "HVPair" in group_facets
    all_hv_tags = set()
    if not hv_faceted:
        for key in chosen_key_list:
            uvs = sorted(groups[key].keys())
            for v, u in [(v, u) for v in uvs for u in uvs if abs(u - 2*v) <= 1e-6]:
                all_hv_tags.add(f"{int(u)}/{int(v)} V")
    hv_uniform = (not hv_faceted) and len(all_hv_tags) == 1

    if _fig1 is not None and plt.fignum_exists(_fig1.number): plt.close(_fig1)
    if _fig2 is not None and plt.fignum_exists(_fig2.number): plt.close(_fig2)
    rows, curves_summary = [], []

    fig1 = plt.figure(figsize=(10.5, 5.5))
    ax1  = fig1.add_subplot(111)
    main_lines_fig1, labels_fig1 = [], []

    for key in last_group_keys:
        if key not in chosen_keys:
            continue
        energy_label = _get_facet(key, "Energy",    group_facets)
        det_name     = _get_facet(key, "Detector",  group_facets)
        ssd_cm       = _get_facet(key, "SSD_cm",    group_facets)
        by_v = groups[key]

        uniq_vs = sorted(by_v.keys())
        if hv_faceted:
            hv_str = _get_facet(key, "HVPair", group_facets) or ""
            try:
                vh_s, vl_s = hv_str.replace(" V", "").split("/")
                hv_pairs = [(float(vl_s), float(vh_s))]
            except Exception:
                hv_pairs = []
        else:
            hv_pairs = [(v, u) for v in uniq_vs for u in uniq_vs if abs(u - 2*v) <= 1e-6]
        if not hv_pairs:
            print(f"Skipping {key} — no HV pairs found: {list(by_v.keys())}")
            continue

        for (v_low, v_high) in hv_pairs:
            lows, highs = by_v.get(v_low, []), by_v.get(v_high, [])
            n_pairs = min(len(lows), len(highs))
            if n_pairs == 0:
                continue

            pair_ds, pair_mL, pair_mH = [], [], []
            for i in range(n_pairs):
                d_c, mL, mH = _align_pair(lows[i][:,[0,1]], highs[i][:,[0,1]])
                idx = np.argsort(d_c)
                d_c, mL, mH = d_c[idx], mL[idx], mH[idx]
                if d_c.size < 2: continue
                pair_ds.append(d_c); pair_mL.append(mL); pair_mH.append(mH)
            if not pair_ds:
                continue

            lo = max(d[0] for d in pair_ds)
            hi = min(d[-1] for d in pair_ds)
            base_idx  = int(np.argmin([np.median(np.diff(d)) if d.size>1 else np.inf for d in pair_ds]))
            base_grid = pair_ds[base_idx]
            base_grid = base_grid[(base_grid >= lo) & (base_grid <= hi)]
            if base_grid.size < 2: continue

            if only_beyond_dmax_var.get():
                dmax_cm = _dmax_for_energy(energy_label)
                if dmax_cm is not None and np.isfinite(dmax_cm):
                    base_grid = base_grid[base_grid >= float(dmax_cm)]
                    if base_grid.size < 2: continue
                else:
                    print(f"[dmax] No mapping for '{energy_label}'.")

            aligned_lows, aligned_highs, pion_list = [], [], []
            for i in range(len(pair_ds)):
                yL = np.interp(base_grid, pair_ds[i], pair_mL[i])
                yH = np.interp(base_grid, pair_ds[i], pair_mH[i])
                aligned_lows.append(np.column_stack([base_grid, yL]))
                aligned_highs.append(np.column_stack([base_grid, yH]))
                pion_list.append(calc_pion(base_grid, yL, yH, v_low, v_high))

            depths   = base_grid
            pion_arr = np.vstack(pion_list)

            base_lbl = varying_labels.get(key, _group_label(key, group_facets))
            hv_tag   = f"{int(v_high)}/{int(v_low)} V"
            if hv_faceted or hv_uniform:
                label = base_lbl or hv_tag
            else:
                label = f"{base_lbl}, {hv_tag}" if base_lbl else hv_tag

            beta = v_high / v_low

            if pion_arr.shape[0] == 1:
                y = pion_arr[0]
                if clamp_pion_var.get(): y = np.maximum(y, 1.0)
                mL = aligned_lows[0][:,1]; mH = aligned_highs[0][:,1]
                with np.errstate(divide='ignore', invalid='ignore'):
                    R0 = np.divide(mH, mL, out=np.full_like(mH, np.nan), where=(mL!=0))
                u_beta_pi = _u_pion_from_beta(R0, beta)
                u_R = R0 * np.sqrt(2.0) * em_std
                den = np.where(np.abs(beta - R0) < 1e-12, np.nan, np.abs(beta - R0))
                u_pion_from_R = np.abs((beta - 1.0) / (den**2)) * u_R
                yerr = 2.0 * np.sqrt(u_beta_pi**2 + u_pion_from_R**2)
                cont = ax1.errorbar(depths, y, yerr=yerr, fmt='-o', capsize=3, label=label)
                scan_type = str(key[-1]).upper()
                curves_summary.append({"x": depths, "mean": y, "yerr": yerr,
                                        "label": label, "energy": energy_label, "scan_type": scan_type})
            else:
                mean = pion_arr.mean(axis=0)
                if clamp_pion_var.get(): mean = np.maximum(mean, 1.0)
                std = pion_arr.std(axis=0, ddof=1)
                mL_mean = np.nanmean(np.vstack([a[:,1] for a in aligned_lows]),  axis=0)
                mH_mean = np.nanmean(np.vstack([a[:,1] for a in aligned_highs]), axis=0)
                with np.errstate(divide='ignore', invalid='ignore'):
                    R0 = np.divide(mH_mean, mL_mean, out=np.full_like(mH_mean, np.nan), where=(mL_mean!=0))
                u_beta_pi = _u_pion_from_beta(R0, beta)
                u_R = R0 * np.sqrt(2.0) * em_std
                den = np.where(np.abs(beta - R0) < 1e-12, np.nan, np.abs(beta - R0))
                u_pion_from_R = np.abs((beta - 1.0) / (den**2)) * u_R
                yerr_spec = 2.0 * np.sqrt(u_beta_pi**2 + u_pion_from_R**2)
                n_reps = pion_arr.shape[0]
                sem  = std / np.sqrt(n_reps)
                t95  = t_dist.ppf(0.975, df=n_reps - 1) if n_reps > 1 else 2.0
                yerr = np.maximum(t95 * sem, yerr_spec)
                cont = ax1.errorbar(depths, mean, yerr=yerr, fmt='-o', capsize=3, label=label)
                scan_type = str(key[-1]).upper() if len(key) > len(group_facets) else "PDD"
                curves_summary.append({"x": depths, "mean": mean, "yerr": yerr,
                                        "label": label, "energy": energy_label, "scan_type": scan_type})

            mainline = cont.lines[0] if hasattr(cont, "lines") else cont[0]
            main_lines_fig1.append(mainline); labels_fig1.append(label)

            for pi in range(pion_arr.shape[0]):
                d_aln  = aligned_lows[pi][:,0]
                yL_aln = aligned_lows[pi][:,1]
                yH_aln = aligned_highs[pi][:,1]
                pion_k = np.maximum(pion_arr[pi], 1.0) if clamp_pion_var.get() else pion_arr[pi]
                for j, depth in enumerate(d_aln):
                    rows.append({"Energy": energy_label, "Detector": det_name, "SSD_cm": ssd_cm,
                                 "HV_low_V": float(v_low), "HV_high_V": float(v_high), "Pair": pi+1,
                                 "Depth_cm": float(depth), "Voltage_V": float(v_low),
                                 "Charge_C": float(yL_aln[j]), "Pion": float(pion_k[j])})
                    rows.append({"Energy": energy_label, "Detector": det_name, "SSD_cm": ssd_cm,
                                 "HV_low_V": float(v_low), "HV_high_V": float(v_high), "Pair": pi+1,
                                 "Depth_cm": float(depth), "Voltage_V": float(v_high),
                                 "Charge_C": float(yH_aln[j]), "Pion": float(pion_k[j])})

    if not rows:
        print("No plottable groups found."); return

    ax1.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    first_key = next(iter(chosen_keys))
    scan_type = str(first_key[-1]).upper() if len(first_key) > len(group_facets) else ""
    subtitle_parts = ([common_str] if common_str else []) + ([next(iter(all_hv_tags))] if hv_uniform else [])
    subtitle = " — ".join(subtitle_parts)

    if "PROFILE" in scan_type:
        ax1.set_xlabel("Off-axis Position [cm]")
        ax1.set_title(f"Pion vs Off-axis Position\n{subtitle}" if subtitle else "Pion vs Off-axis Position", fontsize=10)
    else:
        ax1.set_xlabel("Depth [cm]")
        ax1.set_title(f"Pion vs Depth\n{subtitle}" if subtitle else "Pion vs Depth", fontsize=10)
    ax1.set_ylabel("Pion"); ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 30); ax1.set_ylim(0.98, 1.06)
    _fig1 = fig1
    _attach_click_legend(fig1, main_lines_fig1, labels_fig1)
    fig1.tight_layout()

    # Figure 2 — CF
    fig2 = plt.figure(figsize=(10.5, 5.5))
    ax2  = fig2.add_subplot(111)
    main_lines_fig2, labels_fig2 = [], []

    for item in curves_summary:
        x, y, err, lbl = item["x"], item["mean"], item["yerr"], item["label"]
        en = item.get("energy")
        st = str(item.get("scan_type", "PDD")).upper()
        if y is None or y.size < 2: continue

        if "PROFILE" in st:
            x0 = float(min(max(0.0, np.nanmin(x)), np.nanmax(x)))
            denom = np.interp(x0, x, y)
        elif cf_norm_var.get() == "peak":
            denom = y[int(np.argmax(y))]
        else:
            d = NORM_DEPTH_CM.get(en)
            if d is None:
                print(f"[CF norm] No depth for '{en}'. Skipping."); continue
            denom = np.interp(float(min(max(d, float(np.nanmin(x))), float(np.nanmax(x)))), x, y)

        if not np.isfinite(denom) or denom == 0:
            print(f"[CF norm] Bad norm point for '{lbl}'."); continue

        yn    = y / denom
        err_n = (err / denom) if err is not None else None
        item["cf_mean"] = yn
        item["cf_yerr"] = err_n

        if err_n is not None:
            cont = ax2.errorbar(x, yn, yerr=err_n, fmt='-o', capsize=3, label=lbl)
            mainline = cont.lines[0] if hasattr(cont, "lines") else cont[0]
        else:
            mainline, = ax2.plot(x, yn, marker='o', label=lbl)
        main_lines_fig2.append(mainline); labels_fig2.append(lbl)

    st0 = str(curves_summary[0].get("scan_type", "PDD")).upper() if curves_summary else "PDD"
    if "PROFILE" in st0:
        ax2.set_xlabel("Position [cm]")
        ax2.set_title(f"Profile CF (normalized at CAX)\n{subtitle}" if subtitle else "Profile CF (normalized at CAX)", fontsize=10)
    else:
        ax2.set_xlabel("Depth [cm]")
        ax2.set_title(f"PDD CF\n{subtitle}" if subtitle else "PDD CF", fontsize=10)
    ax2.set_ylabel("CF (normalized)"); ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 30); ax2.set_ylim(0.95, 1.02)
    _fig2 = fig2
    _attach_click_legend(fig2, main_lines_fig2, labels_fig2)
    fig2.tight_layout()

    plt.show()

    last_df = pd.DataFrame(rows, columns=[
        "Energy", "Detector", "SSD_cm", "HV_low_V", "HV_high_V",
        "Pair", "Depth_cm", "Voltage_V", "Charge_C", "Pion"
    ])
    last_meta = {
        "GroupFacets": " / ".join(group_facets),
        "Num_groups":  len(groups),
        "HV_pairs":    ", ".join(sorted({f"{int(r['HV_high_V'])}/{int(r['HV_low_V'])}" for r in rows})),
    }
    facet_names = {"Energy":"Energies","Detector":"Detectors","SSD_cm":"SSDs_cm","FieldSize_cm":"FieldSizes_cm"}
    for facet in group_facets:
        idx  = group_facets.index(facet)
        vals = {k[idx] for k in groups}
        if facet in ("SSD_cm","FieldSize_cm"):
            vals_fmt = sorted({"?" if (v is None or not np.isfinite(v)) else f"{v:g} cm" for v in vals})
        else:
            vals_fmt = sorted({str(v) if v not in (None,"") else "?" for v in vals})
        last_meta[facet_names.get(facet, facet+"s")] = ", ".join(vals_fmt)

# ── CF Report ──────────────────────────────────────────────────────────────────
def show_cf_report():
    global curves_summary
    if not curves_summary:
        print("Plot Pion first."); return
    try:
        depth = float(report_depth_var.get())
    except ValueError:
        print("Invalid depth."); return

    report_rows = []
    for item in curves_summary:
        x   = item["x"]
        y   = item.get("cf_mean", item["mean"])
        err = item.get("cf_yerr", item["yerr"])
        lbl = item["label"]
        en  = item.get("energy", "")
        if depth < float(x.min()) or depth > float(x.max()):
            cf_val = lo_val = hi_val = None
        else:
            cf_val = float(np.interp(depth, x, y))
            if err is not None:
                e = float(np.interp(depth, x, err if err.ndim == 1 else err[1]))
                lo_val, hi_val = cf_val - e, cf_val + e
            else:
                lo_val = hi_val = None
        report_rows.append({"Label": lbl, "Energy": en or "", "CF": cf_val,
                             "Lower (−2σ)": lo_val, "Upper (+2σ)": hi_val})

    win = tk.Toplevel(root)
    win.title(f"CF Report at {depth:.2f} cm")
    win.transient(root)
    cols = ["Label", "Energy", "CF", "Lower (−2σ)", "Upper (+2σ)"]
    tree = ttk.Treeview(win, columns=cols, show="headings", height=min(len(report_rows)+1, 18))
    tree.heading("Label",       text="Group");         tree.column("Label",       width=320, anchor="w")
    tree.heading("Energy",      text="Energy");        tree.column("Energy",      width=80,  anchor="center")
    tree.heading("CF",          text=f"CF @ {depth:.2f} cm"); tree.column("CF", width=110, anchor="center")
    tree.heading("Lower (−2σ)", text="Lower (−2σ)");  tree.column("Lower (−2σ)", width=110, anchor="center")
    tree.heading("Upper (+2σ)", text="Upper (+2σ)");  tree.column("Upper (+2σ)", width=110, anchor="center")
    for r in report_rows:
        tree.insert("", "end", values=(
            r["Label"], r["Energy"],
            f"{r['CF']:.5f}"          if r['CF']          is not None else "out of range",
            f"{r['Lower (−2σ)']:.5f}" if r['Lower (−2σ)'] is not None else "—",
            f"{r['Upper (+2σ)']:.5f}" if r['Upper (+2σ)'] is not None else "—",
        ))
    sb = ttk.Scrollbar(win, orient="vertical", command=tree.yview)
    tree.configure(yscrollcommand=sb.set)
    tree.grid(row=0, column=0, sticky="nsew", padx=(10,0), pady=(10,0))
    sb.grid(row=0, column=1, sticky="ns", padx=(0,10), pady=(10,0))
    win.rowconfigure(0, weight=1); win.columnconfigure(0, weight=1)

    def copy_to_clipboard():
        hdr = "\t".join(["Group","Energy",f"CF @ {depth:.2f} cm","Lower (-2σ)","Upper (+2σ)"])
        lines = [hdr]
        for r in report_rows:
            lines.append("\t".join([
                r["Label"], r["Energy"],
                f"{r['CF']:.5f}"          if r['CF']          is not None else "out of range",
                f"{r['Lower (−2σ)']:.5f}" if r['Lower (−2σ)'] is not None else "",
                f"{r['Upper (+2σ)']:.5f}" if r['Upper (+2σ)'] is not None else "",
            ]))
        win.clipboard_clear(); win.clipboard_append("\n".join(lines))
        print(f"CF report ({len(report_rows)} rows) copied.")

    bf = ttk.Frame(win); bf.grid(row=1, column=0, columnspan=2, pady=8, sticky="e", padx=10)
    ttk.Button(bf, text="Copy to Clipboard", command=copy_to_clipboard).pack(side="right", padx=5)
    ttk.Button(bf, text="Close", command=win.destroy).pack(side="right")

# ── Export ─────────────────────────────────────────────────────────────────────
def export_excel():
    if last_df is None:
        print("Run Plot Pion first."); return
    save_path = filedialog.asksaveasfilename(
        defaultextension=".xlsx", initialfile="pion_results.xlsx",
        filetypes=[("Excel", "*.xlsx")]
    )
    if not save_path: return
    with pd.ExcelWriter(save_path, engine="xlsxwriter") as xw:
        last_df.to_excel(xw, index=False, sheet_name="Pion")
        pd.DataFrame([last_meta]).to_excel(xw, index=False, sheet_name="Meta")
    try:
        if sys.platform.startswith("win"):   os.startfile(save_path)
        elif sys.platform == "darwin":       subprocess.run(["open",     save_path], check=False)
        else:                                subprocess.run(["xdg-open", save_path], check=False)
        print(f"Saved and opened: {save_path}")
    except Exception as e:
        print(f"Saved: {save_path} (auto-open failed: {e})")

# ── GUI ────────────────────────────────────────────────────────────────────────
root = tk.Tk()
root.title("Pion Calculator")

def on_close():
    root.destroy(); root.quit()
root.protocol("WM_DELETE_WINDOW", on_close)

root.columnconfigure(0, weight=1)

top = ttk.Frame(root, padding=10)
top.grid(row=0, column=0, sticky="ew")
top.columnconfigure(1, weight=1)
ttk.Label(top, text="MCC File:").grid(row=0, column=0, sticky="w")
file_entry = ttk.Entry(top, width=50)
file_entry.grid(row=0, column=1, padx=5, sticky="ew")
ttk.Button(top, text="Browse",        command=browse_file).grid(row=0, column=2, padx=5)
ttk.Button(top, text="Browse Folder", command=browse_folder).grid(row=0, column=3, padx=5)

sel_frame = ttk.Frame(root, padding=(10,0,10,10))
sel_frame.grid(row=1, column=0, sticky="nsew")
root.rowconfigure(1, weight=1)
sel_frame.rowconfigure(1, weight=1); sel_frame.columnconfigure(0, weight=1)
ttk.Label(sel_frame, text="Select Groups to Plot:").grid(row=0, column=0, sticky="w")

lb_frame = ttk.Frame(sel_frame)
lb_frame.grid(row=1, column=0, sticky="nsew", pady=(4,0))
lb_frame.rowconfigure(0, weight=1); lb_frame.columnconfigure(0, weight=1)
group_listbox = tk.Listbox(lb_frame, selectmode="extended", height=8)
lb_scroll = ttk.Scrollbar(lb_frame, orient="vertical", command=group_listbox.yview)
group_listbox.configure(yscrollcommand=lb_scroll.set)
group_listbox.grid(row=0, column=0, sticky="nsew")
lb_scroll.grid(row=0, column=1, sticky="ns")

btns = ttk.Frame(sel_frame); btns.grid(row=2, column=0, sticky="w", pady=6)
ttk.Button(btns, text="Select All", command=lambda: group_listbox.select_set(0, tk.END)).pack(side="left", padx=(0,6))
ttk.Button(btns, text="Clear",      command=lambda: group_listbox.selection_clear(0, tk.END)).pack(side="left")

action_frame = ttk.Frame(root, padding=10)
action_frame.grid(row=2, column=0, sticky="ew")
ttk.Button(action_frame, text="Plot Pion",    command=plot_pion).pack(side="left", padx=5)
ttk.Button(action_frame, text="Export Excel", command=export_excel).pack(side="left", padx=5)
ttk.Button(action_frame, text="Grouping…",    command=open_grouping_dialog).pack(side="left", padx=5)
only_beyond_dmax_var = tk.BooleanVar(value=False)
ttk.Checkbutton(action_frame, text="Beyond dmax only", variable=only_beyond_dmax_var).pack(side="left", padx=12)
clamp_pion_var = tk.BooleanVar(value=False)
ttk.Checkbutton(action_frame, text="Set Pion <1 = 1", variable=clamp_pion_var).pack(side="left", padx=4)

report_frame = ttk.Frame(root, padding=(10,0,10,10))
report_frame.grid(row=3, column=0, sticky="ew")
ttk.Label(report_frame, text="Report depth (cm):").pack(side="left")
report_depth_var = tk.StringVar(value="10.0")
ttk.Entry(report_frame, textvariable=report_depth_var, width=7).pack(side="left", padx=5)
ttk.Button(report_frame, text="CF Report", command=show_cf_report).pack(side="left", padx=5)
ttk.Separator(report_frame, orient="vertical").pack(side="left", fill="y", padx=10)
ttk.Label(report_frame, text="CF norm:").pack(side="left")
cf_norm_var = tk.StringVar(value="fixed")
ttk.Radiobutton(report_frame, text="Fixed depth", variable=cf_norm_var, value="fixed").pack(side="left", padx=4)
ttk.Radiobutton(report_frame, text="Peak pion",   variable=cf_norm_var, value="peak").pack(side="left", padx=4)
ttk.Separator(report_frame, orient="vertical").pack(side="left", fill="y", padx=10)
ttk.Label(report_frame, text="EM rel std (%):").pack(side="left")
em_rel_std_var = tk.StringVar(value="0.05")
ttk.Entry(report_frame, textvariable=em_rel_std_var, width=6).pack(side="left", padx=5)

root.mainloop()
