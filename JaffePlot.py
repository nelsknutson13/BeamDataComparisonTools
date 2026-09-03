import numpy as np
import tkinter as tk
from tkinter import filedialog, ttk
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, subprocess

from mcc_utils import (
    available_facets, EM_REL_STD, HV_REL_STD,
    parse_mcc_path, _group_scans,
    _group_label, _get_facet, _split_common_varying,
    align_all_voltages, jaffe_fit, jaffe_pion,
)

# HVPair is a two-voltage concept — not meaningful for Jaffe (all voltages are individual points)
JAFFE_FACETS = [f for f in available_facets if f != "HVPair"]

# Reference dose per pulse (Gy/pulse) and dmax (cm) from BGM calibration, 100 SSD, 10.5 FS
# 1 MU = 1 cGy → divide BGM MU/pulse by 100
BGM_DPP = {
    "6MV":      (2.776e-04, 1.5),
    "6MV FFF":  (7.630e-04, 1.4),
    "8MV FFF":  (1.634e-03, 1.9),
    "10MV":     (2.776e-04, 2.4),
    "10MV FFF": (2.163e-03, 2.2),
    "15MV":     (5.553e-04, 3.0),
}

# ── Module-level state ─────────────────────────────────────────────────────────
last_groups     = None
last_group_keys = []
last_label_map  = {}
last_results    = {}   # key -> {depths, M_inf, alpha, R_sq, pion, label}

group_facets = ["Energy", "Detector", "SSD_cm", "FieldSize_cm"]

_fig_jaffe = None
_fig_pion  = None
_fig_rsq   = None
_fig_dpp   = None
_fig_boag  = None

# ── Legend helper (shared style) ───────────────────────────────────────────────
def _attach_click_legend(fig, lines, labels):
    ax = fig.axes[0]
    short = [lbl if len(lbl) <= 55 else lbl[:52] + "…" for lbl in labels]
    leg = ax.legend(lines, short, fontsize=9, loc='best', framealpha=0.85)
    lined = {}
    for legline, orig in zip(leg.get_lines(), lines):
        legline.set_picker(True)
        lined[legline] = orig
    def on_pick(event):
        orig = lined.get(event.artist)
        if orig is None: return
        vis = not orig.get_visible()
        orig.set_visible(vis)
        event.artist.set_alpha(1.0 if vis else 0.3)
        fig.canvas.draw()
    fig.canvas.mpl_connect("pick_event", on_pick)

# ── Grouping ───────────────────────────────────────────────────────────────────
def open_grouping_dialog():
    win = tk.Toplevel(root)
    win.title("Group by…")
    win.transient(root); win.grab_set()
    ttk.Label(win, text="Select facets:").grid(row=0, column=0, sticky="w", padx=10, pady=(10,6))
    vars_map = {}
    for i, facet in enumerate(JAFFE_FACETS, start=1):
        v = tk.BooleanVar(value=(facet in group_facets))
        ttk.Checkbutton(win, text=facet, variable=v).grid(row=i, column=0, sticky="w", padx=14)
        vars_map[facet] = v
    def apply_and_close():
        chosen = [f for f in JAFFE_FACETS if vars_map[f].get()]
        set_group_facets(chosen or ["Energy", "Detector", "SSD_cm", "FieldSize_cm"])
        win.destroy()
    btns = ttk.Frame(win); btns.grid(row=len(available_facets)+1, column=0, pady=10, padx=10, sticky="e")
    ttk.Button(btns, text="OK",     command=apply_and_close).pack(side="left", padx=6)
    ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="left")

def set_group_facets(chosen):
    global group_facets, last_groups, last_group_keys, last_label_map
    path = file_entry.get().strip()
    if not path: return
    scans = parse_mcc_path(path)
    has_profile = any(len(s) >= 11 and isinstance(s[10], str) and "PROFILE" in s[10].upper() for s in scans)
    if has_profile and "ProfileDepth_cm" not in chosen:
        chosen = chosen + ["ProfileDepth_cm"]
    group_facets = chosen
    last_groups     = _group_scans(scans, group_facets)
    last_group_keys = list(last_groups.keys())
    last_label_map  = {_group_label(k, group_facets): k for k in last_group_keys}
    group_listbox.delete(0, tk.END)
    for k in last_group_keys:
        group_listbox.insert(tk.END, _group_label(k, group_facets))
    group_listbox.select_set(0, tk.END)
    # auto-populate ref DPP from BGM table when a single energy is present
    energies = {_get_facet(k, "Energy", group_facets) for k in last_group_keys}
    if len(energies) == 1:
        energy = next(iter(energies))
        if energy in BGM_DPP:
            dpp_gy, dmax = BGM_DPP[energy]
            try:
                ref_dpp_var.set(f"{dpp_gy:.4e}")
                ref_depth_var.set(f"{dmax}")
            except Exception:
                pass

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

# ── Core analysis ──────────────────────────────────────────────────────────────
def run_jaffe():
    global _fig_jaffe, _fig_pion, _fig_rsq, _fig_dpp, _fig_boag, last_results

    if last_groups is None:
        set_group_facets(group_facets)
    if last_groups is None:
        print("No data loaded."); return

    try:
        v_op = float(op_voltage_var.get())
    except ValueError:
        print("Invalid operating voltage."); return
    try:
        diag_depths = [float(x.strip()) for x in diag_depth_var.get().split(",") if x.strip()]
        if not diag_depths:
            raise ValueError
    except ValueError:
        print("Invalid diagnostic depths — enter one or more comma-separated values."); return

    ref_dpp_str = ref_dpp_var.get().strip()
    ref_dpp     = float(ref_dpp_str) if ref_dpp_str else None   # Gy/pulse at cal_ssd; None = skip DPP figure
    try:
        cal_ssd = float(cal_ssd_var.get())
    except ValueError:
        cal_ssd = 100.0
    try:
        ref_depth = float(ref_depth_var.get())
    except ValueError:
        ref_depth = diag_depths[0]

    sel_idxs = group_listbox.curselection()
    if not sel_idxs:
        chosen_keys = set(last_group_keys)
    else:
        labels = [group_listbox.get(i) for i in sel_idxs]
        chosen_keys = {last_label_map[lbl] for lbl in labels if lbl in last_label_map}
    if not chosen_keys:
        print("No groups selected."); return

    chosen_key_list = [k for k in last_group_keys if k in chosen_keys]
    common_str, varying_labels = _split_common_varying(chosen_key_list, group_facets)
    subtitle = common_str

    print(f"Running Jaffe fit on {len(chosen_keys)} group(s). Operating voltage: {v_op} V")

    last_results = {}

    for fig in (_fig_jaffe, _fig_pion, _fig_rsq, _fig_dpp, _fig_boag):
        if fig is not None and plt.fignum_exists(fig.number):
            plt.close(fig)

    fig_j = plt.figure(figsize=(8, 5))
    ax_j  = fig_j.add_subplot(111)
    fig_p = plt.figure(figsize=(10.5, 5.5))
    ax_p  = fig_p.add_subplot(111)
    fig_r = plt.figure(figsize=(10.5, 4.5))
    ax_r  = fig_r.add_subplot(111)
    fig_d = plt.figure(figsize=(10.5, 4.5)) if ref_dpp else None
    ax_d  = fig_d.add_subplot(111) if fig_d else None
    fig_b = plt.figure(figsize=(7, 5)) if ref_dpp else None
    ax_b  = fig_b.add_subplot(111) if fig_b else None

    lines_j, lbls_j = [], []
    lines_p, lbls_p = [], []
    lines_r, lbls_r = [], []
    lines_d, lbls_d = [], []
    lines_b, lbls_b = [], []
    boag_pts = []   # (D_p_Gy, b, label) — collected for α_Boag extraction

    for key in last_group_keys:
        if key not in chosen_keys:
            continue
        by_v  = last_groups[key]
        label = varying_labels.get(key, _group_label(key, group_facets))

        depths, aligned, aligned_sem = align_all_voltages(by_v)
        if depths is None or len(aligned) < 2:
            print(f"[skip] {label}: need ≥2 voltages, got {len(aligned)}."); continue

        M_inf, alpha, R_sq = jaffe_fit(depths, aligned)
        pion = jaffe_pion(M_inf, alpha, v_op)

        # b × M_inf = α_Boag × DPP  (units: volts, a characteristic recombination voltage)
        b_Minf = alpha * M_inf   # shape (n_depths,)

        dpp = None
        if ref_dpp is not None:
            # apply inverse-square correction: BGM is at cal_ssd, scan may differ
            energy_key = _get_facet(key, "Energy", group_facets) or ""
            dmax_bgm   = BGM_DPP.get(energy_key, (None, ref_depth))[1]
            try:
                ssd_scan = float(_get_facet(key, "SSD_cm", group_facets) or cal_ssd)
            except (TypeError, ValueError):
                ssd_scan = cal_ssd
            isf = ((cal_ssd + dmax_bgm) / (ssd_scan + dmax_bgm)) ** 2
            ref_dpp_corrected = ref_dpp * isf
            print(f"  [{label}] ISF: cal_SSD={cal_ssd}, scan_SSD={ssd_scan:.1f}, "
                  f"dmax={dmax_bgm} → ×{isf:.4f}, "
                  f"ref_DPP {ref_dpp*1e4:.3f}→{ref_dpp_corrected*1e4:.3f} ×10⁻⁴ Gy/pulse")
            # Use highest-voltage PDD as dose proxy (least recombination)
            # DPP(depth) = DPP_ref × M(depth,V_max) / M(ref_depth,V_max)
            v_max     = max(aligned.keys())
            m_vmax    = aligned[v_max]                        # charge at all depths, highest V
            ref_idx   = int(np.argmin(np.abs(depths - ref_depth)))
            m_ref_val = m_vmax[ref_idx]
            if np.isfinite(m_ref_val) and m_ref_val != 0:
                dpp = ref_dpp_corrected * m_vmax / m_ref_val  # Gy/pulse at each depth

        last_results[key] = {
            "depths": depths, "M_inf": M_inf, "alpha": alpha,
            "R_sq": R_sq, "pion": pion, "dpp": dpp, "label": label,
            "aligned": aligned,
        }

        # ── collect (D_p, b) for α_Boag extraction figure (R²≥0.99 only) ──────
        if dpp is not None:
            for i_d in range(len(depths)):
                b_val  = alpha[i_d]
                dp_val = dpp[i_d]
                r2_val = R_sq[i_d]
                if (np.isfinite(b_val) and np.isfinite(dp_val) and dp_val > 0
                        and np.isfinite(r2_val) and r2_val >= 0.99):
                    boag_pts.append((dp_val, b_val, label))

        # ── Jaffe diagnostic: one line per selected depth, normalized to lowest V ─
        voltages = np.array(sorted(aligned.keys()))
        inv_V    = 1.0 / voltages
        v_low    = voltages[0]   # lowest voltage = most recombination → normalize to 1.0
        multi_group = len(chosen_key_list) > 1
        for d_target in diag_depths:
            d_idx    = int(np.argmin(np.abs(depths - d_target)))
            actual_d = depths[d_idx]
            inv_M_raw = np.array([1.0 / aligned[v][d_idx] for v in voltages])
            ok        = np.isfinite(inv_M_raw) & (inv_M_raw != 0)
            norm_val  = inv_M_raw[voltages == v_low][0] if (voltages == v_low).any() else inv_M_raw[ok][0]
            if norm_val == 0 or not np.isfinite(norm_val):
                continue
            inv_M_n  = inv_M_raw / norm_val   # normalized so lowest-V point = 1.0
            # label: use DPP if available, otherwise fall back to depth
            if dpp is not None and np.isfinite(dpp[d_idx]):
                dpp_cGy = dpp[d_idx] * 100   # Gy → cGy
                dpp_str = f"{dpp_cGy:.4f} cGy/pulse"
                d_lbl = f"{label} — {dpp_str}" if multi_group else dpp_str
            else:
                d_lbl = f"{label} @ {actual_d:.1f} cm" if multi_group else f"{actual_d:.1f} cm"

            # ── uncertainty propagation ───────────────────────────────────────
            # σ(M) = max(SEM_from_repeats, EM_REL_STD × M)  per voltage
            M_raw   = np.array([aligned[v][d_idx] for v in voltages])
            sem_M   = np.array([aligned_sem[v][d_idx] for v in voltages])
            sigma_M = np.maximum(sem_M, EM_REL_STD * M_raw)
            # normalization point
            norm_v_idx  = int(np.argmin(np.abs(voltages - v_low)))
            sigma_M_nrm = sigma_M[norm_v_idx]
            M_nrm       = M_raw[norm_v_idx]
            # propagate: σ(y_i)/y_i = sqrt((σ_Mi/M_i)² + (σ_Mnorm/M_norm)²)
            sigma_y = inv_M_n * np.sqrt((sigma_M / M_raw)**2 + (sigma_M_nrm / M_nrm)**2)
            # x-bars from HV spec (tiny but correct)
            sigma_x = HV_REL_STD * inv_V

            eb = ax_j.errorbar(inv_V[ok], inv_M_n[ok],
                               yerr=sigma_y[ok], xerr=sigma_x[ok],
                               fmt='o', markersize=6, capsize=3, linewidth=1,
                               zorder=3, label=d_lbl)
            color = eb[0].get_color()
            if ok.sum() >= 2:
                coeffs = np.polyfit(inv_V[ok], inv_M_n[ok], 1)
                x_fit  = np.linspace(inv_V[ok].min() * 0.95, inv_V[ok].max() * 1.05, 100)
                ax_j.plot(x_fit, np.polyval(coeffs, x_fit), color=color,
                          linewidth=1.2, alpha=0.7, linestyle='--')
            lines_j.append(eb[0]); lbls_j.append(d_lbl)

        # ── Pion vs depth ─────────────────────────────────────────────────────
        valid = np.isfinite(pion)
        if valid.any():
            line, = ax_p.plot(depths[valid], pion[valid], '-o', markersize=3, label=label)
            lines_p.append(line); lbls_p.append(label)

        # ── R² vs depth ───────────────────────────────────────────────────────
        valid_r = np.isfinite(R_sq)
        if valid_r.any():
            line, = ax_r.plot(depths[valid_r], R_sq[valid_r], '-', linewidth=1.5, label=label)
            lines_r.append(line); lbls_r.append(label)

        # ── Dose per pulse vs depth ───────────────────────────────────────────
        if ax_d is not None and dpp is not None:
            valid_d = np.isfinite(dpp)
            if valid_d.any():
                dpp_mGy = dpp[valid_d] * 1e3   # convert Gy → mGy
                line, = ax_d.plot(depths[valid_d], dpp_mGy, '-o', markersize=3, label=label)
                lines_d.append(line); lbls_d.append(label)

    if not last_results:
        print("No groups produced valid Jaffe fits."); return

    # ── Finish Jaffe diagnostic ───────────────────────────────────────────────
    depth_str = ", ".join(f"{d:.1f}" for d in diag_depths)
    ax_j.axhline(1.0, color="gray", linewidth=0.8, linestyle="--")
    ax_j.set_xlabel("1 / Voltage  (V⁻¹)")
    ax_j.set_ylabel(f"1/Q  (normalized to 1/Q at {int(v_low)} V)")
    ax_j.set_title(f"Jaffe Plot — depths {depth_str} cm\n{subtitle}" if subtitle
                   else f"Jaffe Plot — depths {depth_str} cm", fontsize=10)
    ax_j.grid(True, alpha=0.3)
    _attach_click_legend(fig_j, lines_j, lbls_j)
    fig_j.tight_layout()

    # ── Finish Pion vs depth ──────────────────────────────────────────────────
    ax_p.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax_p.set_xlabel("Depth [cm]"); ax_p.set_ylabel("Pion")
    ax_p.set_title(f"Pion vs Depth (Jaffe, V_op = {v_op:.0f} V)\n{subtitle}" if subtitle
                   else f"Pion vs Depth (Jaffe, V_op = {v_op:.0f} V)", fontsize=10)
    ax_p.set_xlim(0, 30); ax_p.set_ylim(0.98, 1.06)
    ax_p.grid(True, alpha=0.3)
    _attach_click_legend(fig_p, lines_p, lbls_p)
    fig_p.tight_layout()

    # ── Finish R² vs depth ────────────────────────────────────────────────────
    ax_r.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax_r.axhline(0.99, color="orange", linestyle=":", linewidth=0.8)
    ax_r.set_xlabel("Depth [cm]"); ax_r.set_ylabel("R²")
    ax_r.set_title(f"Jaffe Fit R² vs Depth\n{subtitle}" if subtitle
                   else "Jaffe Fit R² vs Depth", fontsize=10)
    ax_r.set_xlim(0, 30); ax_r.set_ylim(0.9, 1.005)
    ax_r.grid(True, alpha=0.3)
    _attach_click_legend(fig_r, lines_r, lbls_r)
    fig_r.tight_layout()

    # ── Finish DPP vs depth ───────────────────────────────────────────────────
    if fig_d is not None and lines_d:
        ax_d.set_xlabel("Depth [cm]"); ax_d.set_ylabel("Dose per pulse [mGy/pulse]")
        ax_d.set_title(f"Dose per Pulse vs Depth (ref {ref_dpp*1e3:.3f} mGy/pulse @ {ref_depth:.1f} cm)"
                       + (f"\n{subtitle}" if subtitle else ""), fontsize=10)
        ax_d.set_xlim(0, 30)
        ax_d.grid(True, alpha=0.3)
        _attach_click_legend(fig_d, lines_d, lbls_d)
        fig_d.tight_layout()

    # ── Finish α_Boag extraction: b vs D_p across all depths × groups ────────
    if fig_b is not None and boag_pts:
        dp_all = np.array([p[0] for p in boag_pts]) * 1e3   # Gy → mGy
        b_all  = np.array([p[1] for p in boag_pts])
        # forced-through-origin OLS: min Σ(b_i - α·Dp_i)²  →  α = Σ(b·Dp)/Σ(Dp²)
        alpha_boag = np.dot(b_all, dp_all) / np.dot(dp_all, dp_all)
        # scatter per group
        groups_b = {}
        for dp, b, lbl in boag_pts:
            if lbl not in groups_b:
                groups_b[lbl] = ([], [])
            groups_b[lbl][0].append(dp * 1e3)
            groups_b[lbl][1].append(b)
        for lbl, (dp_list, b_list) in groups_b.items():
            line, = ax_b.plot(dp_list, b_list, 'o', markersize=6, alpha=0.8,
                              linewidth=0, label=lbl)
            lines_b.append(line); lbls_b.append(lbl)
        # overlay forced-origin fit line
        dp_fit = np.linspace(0, dp_all.max() * 1.08, 200)
        ax_b.plot(dp_fit, alpha_boag * dp_fit, 'k--', linewidth=1.5,
                  label=f"b = {alpha_boag:.4f}·Dp  (n={len(boag_pts)})")
        print(f"  α_Boag (forced-origin slope) = {alpha_boag:.5f} V·pulse/mGy  "
              f"from {len(boag_pts)} depth×group points")
        ax_b.set_xlabel("Dose per pulse [mGy/pulse]")
        ax_b.set_ylabel("Jaffé slope  b  [V / charge-unit]")
        ax_b.set_title("α_Boag Extraction: Jaffé Slope vs Dose per Pulse"
                       + (f"\n{subtitle}" if subtitle else ""), fontsize=10)
        ax_b.set_xlim(left=0); ax_b.set_ylim(bottom=0)
        ax_b.grid(True, alpha=0.3)
        _attach_click_legend(fig_b, lines_b, lbls_b)
        fig_b.tight_layout()
    elif fig_b is not None:
        plt.close(fig_b); fig_b = None

    _fig_jaffe, _fig_pion, _fig_rsq, _fig_dpp, _fig_boag = fig_j, fig_p, fig_r, fig_d, fig_b
    plt.show()

# ── Export ─────────────────────────────────────────────────────────────────────
def export_excel():
    if not last_results:
        print("Run Jaffe first."); return
    try:
        v_op = float(op_voltage_var.get())
    except ValueError:
        v_op = float("nan")

    rows = []
    for key, res in last_results.items():
        lbl    = res["label"]
        depths = res["depths"]
        for i, d in enumerate(depths):
            dpp_i = float(res["dpp"][i]) if (res["dpp"] is not None and np.isfinite(res["dpp"][i])) else None
            rows.append({
                "Label":             lbl,
                "Depth_cm":          float(d),
                "M_inf":             float(res["M_inf"][i])  if np.isfinite(res["M_inf"][i])  else None,
                "b_jaffe_slope":     float(res["alpha"][i])  if np.isfinite(res["alpha"][i])  else None,
                "R_sq":              float(res["R_sq"][i])   if np.isfinite(res["R_sq"][i])   else None,
                f"Pion_at_{int(v_op)}V": float(res["pion"][i]) if np.isfinite(res["pion"][i]) else None,
                "DPP_mGy_pulse":     dpp_i * 1e3 if dpp_i is not None else None,
            })

    df = pd.DataFrame(rows)
    save_path = filedialog.asksaveasfilename(
        defaultextension=".xlsx", initialfile="jaffe_results.xlsx",
        filetypes=[("Excel", "*.xlsx")]
    )
    if not save_path: return
    with pd.ExcelWriter(save_path, engine="xlsxwriter") as xw:
        df.to_excel(xw, index=False, sheet_name="Jaffe")
    try:
        if sys.platform.startswith("win"):   os.startfile(save_path)
        elif sys.platform == "darwin":       subprocess.run(["open",     save_path], check=False)
        else:                                subprocess.run(["xdg-open", save_path], check=False)
        print(f"Saved and opened: {save_path}")
    except Exception as e:
        print(f"Saved: {save_path} (auto-open failed: {e})")

# ── GUI ────────────────────────────────────────────────────────────────────────
root = tk.Tk()
root.title("Jaffe Plot")

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
ttk.Label(sel_frame, text="Select Groups:").grid(row=0, column=0, sticky="w")

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
ttk.Button(action_frame, text="Run Jaffe",    command=run_jaffe).pack(side="left", padx=5)
ttk.Button(action_frame, text="Export Excel", command=export_excel).pack(side="left", padx=5)
ttk.Button(action_frame, text="Grouping…",    command=open_grouping_dialog).pack(side="left", padx=5)

param_frame = ttk.Frame(root, padding=(10,0,10,10))
param_frame.grid(row=3, column=0, sticky="ew")
ttk.Label(param_frame, text="Operating voltage (V):").pack(side="left")
op_voltage_var = tk.StringVar(value="300")
ttk.Entry(param_frame, textvariable=op_voltage_var, width=7).pack(side="left", padx=5)
ttk.Separator(param_frame, orient="vertical").pack(side="left", fill="y", padx=10)
ttk.Label(param_frame, text="Depths for Jaffe plot (cm, comma-sep):").pack(side="left")
diag_depth_var = tk.StringVar(value="1.5, 5, 10, 20, 30")
ttk.Entry(param_frame, textvariable=diag_depth_var, width=22).pack(side="left", padx=5)

ref_frame = ttk.Frame(root, padding=(10,0,10,10))
ref_frame.grid(row=4, column=0, sticky="ew")
ttk.Label(ref_frame, text="Ref DPP (Gy/pulse, optional):").pack(side="left")
ref_dpp_var = tk.StringVar(value="")
ttk.Entry(ref_frame, textvariable=ref_dpp_var, width=12).pack(side="left", padx=5)
ttk.Separator(ref_frame, orient="vertical").pack(side="left", fill="y", padx=10)
ttk.Label(ref_frame, text="at ref depth (cm):").pack(side="left")
ref_depth_var = tk.StringVar(value="1.5")
ttk.Entry(ref_frame, textvariable=ref_depth_var, width=7).pack(side="left", padx=5)
ttk.Separator(ref_frame, orient="vertical").pack(side="left", fill="y", padx=10)
ttk.Label(ref_frame, text="Cal SSD (cm):").pack(side="left")
cal_ssd_var = tk.StringVar(value="100")
ttk.Entry(ref_frame, textvariable=cal_ssd_var, width=6).pack(side="left", padx=5)

root.mainloop()
