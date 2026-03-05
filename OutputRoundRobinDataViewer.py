# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 11:06:44 2025
OutputRound Robin data viewer
@author: nknutson
"""

# -*- coding: utf-8 -*-
"""
Output Round Robin data viewer (minimal, clean boxplots)
- No Institution exclusion
- No titles
- No point markers / legend
- Fixed box order: IROC, Consortium Form, Consortium Audit
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf

# ---------- Config ----------
DEFAULT_PATH = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\OutputRoundRobinData.xlsx"
EXPECTED_ENERGIES = ["6X", "10X", "15X", "6FFF", "8FFF", "10FFF"]
SYSTEM_ORDER = ["IROC", "RDS", "Consortium Audit", "Institution Reported (Consortium Form)", "Consortium Audit (Form-Corrected)"]
ORDER = {name: i for i, name in enumerate(SYSTEM_ORDER)}
OUTLIER_MIN_N = 15   # only show 1.5*IQR fliers when N >= this

MARKER_COLOR = "black"   # or "#444", "0.3", etc.

# ---------- Data helpers ----------
def read_table(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)  # single sheet assumed
    elif ext in (".csv", ".txt"):
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, sep="\t")
    else:
        raise ValueError(f"Unsupported file type: {ext}")



def to_long(df: pd.DataFrame):
    energies = [e for e in EXPECTED_ENERGIES if e in df.columns]
    if not energies:
        raise ValueError("No energy columns found (need one of: " + ", ".join(EXPECTED_ENERGIES) + ")")
    for c in ["Date", "SN", "System", "Chamber"]:
        if c not in df.columns:
            df[c] = np.nan
    for e in energies:
        df[e] = pd.to_numeric(df[e], errors="coerce")
    long = df.melt(
        id_vars=["Date", "SN", "System", "Chamber"],
        value_vars=energies,
        var_name="Energy",
        value_name="Ratio"
    ).dropna(subset=["Ratio"]).reset_index(drop=True)
    # parse dates if present (optional)
    if long["Date"].dtype == "object":
        with pd.option_context("mode.chained_assignment", None):
            long["Date"] = pd.to_datetime(long["Date"], errors="coerce")
    return long, energies


def _draw_tukey_fliers_threshold(ax, series, positions, min_n, color, markersize=10, zorder=6):
    """
    Draw Tukey (1.5*IQR) fliers ONLY for boxes with N >= min_n.
    Works across matplotlib versions regardless of boxplot(showfliers=...) behavior.
    """
    for vals, x in zip(series, positions):
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) < min_n:
            continue

        q1 = np.percentile(vals, 25)
        q3 = np.percentile(vals, 75)
        iqr = q3 - q1
        if not np.isfinite(iqr) or iqr <= 0:
            continue

        lo = q1 - 1.5 * iqr
        hi = q3 + 1.5 * iqr
        fliers = vals[(vals < lo) | (vals > hi)]
        if fliers.size == 0:
            continue

        ax.plot(
            np.full(fliers.size, x, dtype=float), fliers,
            linestyle="None",
            marker=r'$\ast$',
            color=color,
            markersize=markersize,
            zorder=zorder
        )


def _system_boxplot(long: pd.DataFrame, show_dates: bool = False, show_sn_labels: bool = False, show_energy_labels: bool = False):


    data = long.copy()
    data["System"] = data["System"].astype(str).str.strip()

    # Keep only systems in your fixed order and only those present
    systems_present = [s for s in SYSTEM_ORDER if s in set(data["System"].dropna())]
    if not systems_present:
        raise ValueError("No recognized systems found to plot.")

    # Build series in fixed left→right order
    series = [data.loc[data["System"] == sys, "Ratio"].dropna().values for sys in systems_present]
    positions = list(range(len(systems_present)))
    pos_map = {sys: i for i, sys in enumerate(systems_present)}

    # Figure/axes
    fig = plt.figure(figsize=(8.5, 5))
    ax = fig.add_subplot(111)

    # Boxplot styling to match existing grouped plots
    bp = ax.boxplot(
        series,
        positions=positions,
        widths=0.6,
        showmeans=False,
        showfliers=False,  # <- always off; we draw fliers ourselves
        patch_artist=False,
        medianprops=dict(color='black', linewidth=1.2),
    )

    # Draw outlier stars only when N >= OUTLIER_MIN_N
    _draw_tukey_fliers_threshold(ax, series, positions, OUTLIER_MIN_N, MARKER_COLOR, markersize=10, zorder=6)

    # Baseline at 1.0
    ax.axhline(1.0, linestyle="--", color="black", linewidth=1.0)

    # Overlay sample points with same marker map as grouped plot
    marker_map = {
        "IROC": "D",
        "RDS": "o",
        "Institution Reported (Consortium Form)": "s",
        "Consortium Audit": "X",
        "Consortium Audit (Form-Corrected)": "^",
    }
    for sys in systems_present:
        sub = data.loc[data["System"] == sys, ["Ratio", "Date", "SN", "Energy"]].dropna(subset=["Ratio"])
        vals = sub["Ratio"].values
        dates = sub["Date"].values
        sns = sub["SN"].values
        ens = sub["Energy"].values

        if len(vals) == 0:
            continue
        x = pos_map[sys]
        #xs = np.full(len(vals), x, dtype=float)
        xs = x + np.random.uniform(-0.01, 0.01, size=len(vals))

        ax.scatter(
            xs, vals,
            marker=marker_map.get(sys, "D"),
            s=40, alpha=0.9,
            facecolors="none", edgecolors=MARKER_COLOR, linewidths=1.0,
            zorder=5
        )
        if show_sn_labels or show_energy_labels:
            for xi, yi, sn, en in zip(xs, vals, sns, ens):
                parts = []
                if show_sn_labels and pd.notna(sn): parts.append(f"SN {sn}")
                if show_energy_labels and pd.notna(en): parts.append(str(en))
                if parts: ax.text(xi + 0.01, yi + 0.00015, " | ".join(parts), fontsize=7, alpha=0.8)
        



    # Labels (no title)
    ax.set_xticks(positions)
    ax.set_xticklabels(systems_present, rotation=0)
    ax.set_ylabel("Output [cGy/MU]")

    # Legend to match grouped plot
    handles = [
        Line2D([], [], linestyle="None",
               marker=marker_map.get(sys, "D"), markersize=10,
               markerfacecolor="none", markeredgecolor=MARKER_COLOR, color=MARKER_COLOR,
               label=sys)
        for sys in systems_present
    ]
    handles.append(Line2D([0], [0], linestyle="--", color="black", label="Institution reported (1.0 cGy/MU)"))
    handles.append(Line2D([], [], linestyle="None", marker=r'$\ast$', markersize=7,
                          color=MARKER_COLOR, label=r"Outlier (1.5$\times$ IQR, only if N $\geq$ %d)" % OUTLIER_MIN_N))
    ax.legend(handles=handles, loc="best", frameon=True)

    plt.tight_layout()
    plt.show()


# ---------- Plot ----------
def _grouped_boxplot(long: pd.DataFrame, group_col: str, energies_order,
                     show_dates: bool = False,
                     show_sn_labels: bool = False,
                     show_energy_labels: bool = False):


    data = long.copy()
    data["System"] = data["System"].astype(str).str.strip()
    data = data[data["System"].isin(ORDER.keys())]  # keep only your ordered systems

    # Fixed system order; include only those present
    present = set(data["System"].astype(str).str.strip().unique())
    systems = sorted(present, key=lambda s: ORDER.get(s, 999))

    if not systems:
        raise ValueError("No systems to plot (none of IROC / Consortium Form / Consortium Audit present).")

    # Choose grouping axis
    if group_col == "Energy":
        groups = [e for e in energies_order if e in set(data["Energy"])]
        xlabels = groups
    else:
        groups = sorted([g for g in pd.unique(data["SN"]) if not (isinstance(g, float) and np.isnan(g))])
        xlabels = [f"SN {int(g)}" for g in groups]

    # Build box data in fixed left→right system order inside each group
    series, positions, xticks, xtlabs = [], [], [], []
    pos_map = {}  # (group, system) -> x position for overlay
    n_sys, gap, base = len(systems), 1.0, 0.0
    for g in groups:
        xticks.append(base + (n_sys - 1) / 2)
        xtlabs.append(xlabels[groups.index(g)])
        for j, sys in enumerate(systems):
            sub = data[(data[group_col] == g) & (data["System"] == sys)]
            vals = sub["Ratio"].values
            series.append(vals)
            x = base + j
            positions.append(x)
            pos_map[(g, sys)] = x
        base += n_sys + gap

    # Plot
    fig = plt.figure(figsize=(11.5, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(
        series,
        positions=positions,
        widths=0.6,
        showmeans=False,
        showfliers=False,  # <- always off; we draw fliers ourselves
        patch_artist=False,
        medianprops=dict(color='black', linewidth=1.2),
    )

    # Draw outlier stars only when N >= OUTLIER_MIN_N
    _draw_tukey_fliers_threshold(ax, series, positions, OUTLIER_MIN_N, MARKER_COLOR, markersize=10, zorder=6)

    # Baseline at 1.0 (and use black so legend matches)
    ax.axhline(1.0, linestyle="--", color="black", linewidth=1.0)

    # Overlay sample points with distinct markers per system (basic markers)
    marker_map = {"IROC": "D", "RDS": "o", "Institution Reported (Consortium Form)": "s",
                  "Consortium Audit": "X", "Consortium Audit (Form-Corrected)": "^"}

    for g in groups:
        for sys in systems:
            sub = data[(data[group_col] == g) & (data["System"] == sys)][["Ratio", "Date", "SN", "Energy"]].dropna(subset=["Ratio"])

            vals = sub["Ratio"].values
            dates = sub["Date"].values
            sns = sub["SN"].values
            ens = sub["Energy"].values

            if len(vals) == 0:
                continue
            x = pos_map[(g, sys)]
            #xs = np.full(len(vals), x, dtype=float)
            xs = x + np.random.uniform(-0.01, 0.01, size=len(vals))

            ax.scatter(xs, vals,
                       marker=marker_map.get(sys, "D"), s=40, alpha=0.9,
                       facecolors="none", edgecolors=MARKER_COLOR, linewidths=1.0,
                       zorder=5)  # draw above boxplot fliers (we draw fliers at zorder=6)
            if show_dates:
                for xi, yi, di in zip(xs, vals, dates):
                    if pd.notna(di):
                        ax.text(
                            xi, yi + 0.00015,
                            pd.to_datetime(di).strftime("%Y-%m-%d"),
                            fontsize=7, alpha=0.8
                        )
            
            if show_sn_labels and group_col == "Energy":
                for xi, yi, sn in zip(xs, vals, sns):
                    if pd.notna(sn):
                        ax.text(xi + 0.01, yi + 0.00015, f"SN {sn}", fontsize=7, alpha=0.8)
            if show_energy_labels and group_col == "SN":
                for xi, yi, en in zip(xs, vals, ens):
                    if pd.notna(en):
                        ax.text(xi + 0.01, yi + 0.00015, str(en), fontsize=7, alpha=0.8)



    # Labels (no title)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtlabs)
    ax.set_ylabel("Output [cGy/MU]")

    # Legend: markers for each system + dashed baseline
    handles = [Line2D([], [], linestyle="None",
                      marker=marker_map.get(sys, "D"), markersize=10,
                      markerfacecolor="none", markeredgecolor=MARKER_COLOR, color=MARKER_COLOR,
                      label=sys) for sys in systems]

    handles.append(Line2D([0], [0], linestyle="--", color="black", label="Institution reported (1.0 cGy/MU)"))
    handles.append(Line2D([], [], linestyle="None", marker=r'$\ast$', markersize=7,
                          color=MARKER_COLOR, label=r"Outlier (1.5$\times$ IQR, only if N $\geq$ %d)" % OUTLIER_MIN_N))

    ax.legend(handles=handles, loc="best", frameon=True)

    plt.tight_layout()
    plt.show()


def analysis_B_mixed_effects(long_filt: pd.DataFrame, verbose: bool = True):
    """
    Analysis B: Mixed-effects model on Delta = Ratio - 1.00
      Delta ~ 0 + C(System)
      Random intercepts:
        - groups = SN
        - variance components for Date (optional) and Energy (optional)
    """
    df = long_filt.copy()

    # ---- Clean required fields ----
    if "Ratio" not in df.columns or "System" not in df.columns or "SN" not in df.columns:
        raise ValueError("Need columns: 'Ratio', 'System', 'SN' (and optionally 'Date', 'Energy').")

    df["Ratio"] = pd.to_numeric(df["Ratio"], errors="coerce")
    df = df.dropna(subset=["Ratio", "System", "SN"]).copy()

    # System as plain python string, trimmed
    df["System"] = df["System"].astype(str).str.strip()

    # SN: keep as-is, but make it categorical-ish & stable
    # (MixedLM groups can be numbers or strings; strings are safest)
    df["SN"] = df["SN"].astype(str).str.strip()

    # Delta from target
    df["Delta"] = df["Ratio"] - 1.0

    # ---- Optional clustering factors as categoricals ----
    vc = {}

    # DateCat: use date only, string label; ignore NaT
    if "Date" in df.columns:
        dt = pd.to_datetime(df["Date"], errors="coerce")
        df["DateCat"] = dt.dt.date.astype(object)  # python date or NaT->NaN
        # convert to string labels but keep NaN as NaN (not "NaT")
        df["DateCat"] = df["DateCat"].where(pd.notna(df["DateCat"]), np.nan)

        if df["DateCat"].notna().sum() > 0 and pd.Series(df["DateCat"]).nunique(dropna=True) > 1:
            vc["Date"] = "0 + C(DateCat)"

    # EnergyCat: string labels
    if "Energy" in df.columns:
        df["EnergyCat"] = df["Energy"].astype(str).str.strip()
        # treat empty strings as missing
        df.loc[df["EnergyCat"].isin(["", "nan", "None"]), "EnergyCat"] = np.nan

        if df["EnergyCat"].notna().sum() > 0 and df["EnergyCat"].nunique(dropna=True) > 1:
            vc["Energy"] = "0 + C(EnergyCat)"

    # ---- Fit model ----
    # No intercept => each system coefficient is that system's mean Delta
    formula = "Delta ~ 0 + C(System)"
    model = smf.mixedlm(
        formula,
        df,
        groups=df["SN"],
        vc_formula=vc if vc else None,
        re_formula="1"
    )
    res = model.fit(reml=True, method="lbfgs")

    if verbose:
        print("\n=== Analysis B: Mixed Effects on Delta = Ratio - 1.00 ===")
        print(res.summary())

        # system-specific deltas & 95% CI
        params = res.params
        bse = res.bse

        rows = []
        for k in params.index:
            if k.startswith("C(System)["):
                sys_name = k.split("[", 1)[1].rstrip("]")
                mu_delta = float(params[k])
                se = float(bse[k])
                ci_lo = mu_delta - 1.96 * se
                ci_hi = mu_delta + 1.96 * se
                rows.append({
                    "System": sys_name,
                    "Mean_Delta": mu_delta,
                    "CI95_Delta_L": ci_lo,
                    "CI95_Delta_U": ci_hi,
                    "Mean_Ratio": mu_delta + 1.0,
                    "CI95_Ratio_L": ci_lo + 1.0,
                    "CI95_Ratio_U": ci_hi + 1.0
                })

        out = pd.DataFrame(rows).sort_values("System")
        print("\nSystem-level estimated bias (vs 1.00) with 95% CI:")
        print(out.to_string(index=False))

    return res


def make_plots(df: pd.DataFrame,
               show_system: bool = True,
               show_sn: bool = True,
               show_energy: bool = True,
               sn_filter=None,
               energy_filter=None,
               system_filter=None,
               show_dates: bool = False,
               show_sn_labels: bool = False,
               show_energy_labels: bool = False):

    long, energies = to_long(df)

    if not (show_system or show_sn or show_energy):
        raise ValueError("No plots selected. Please enable at least one plot type.")

    mask = np.ones(len(long), dtype=bool)

    if sn_filter:
        mask &= long["SN"].isin(sn_filter)

    if energy_filter:
        mask &= long["Energy"].isin(energy_filter)

    if system_filter:
        mask &= long["System"].astype(str).str.strip().isin(system_filter)

    long_filt = long[mask]
    if long_filt.empty:
        raise ValueError("No data left after applying SN/Energy/System filters.")

    # Analysis B
    try:
        analysis_B_mixed_effects(long_filt, verbose=True)
    except Exception as e:
        print(f"Analysis B failed: {e}")

    # Plots
    if show_system:
        _system_boxplot(long_filt, show_dates=show_dates, show_sn_labels=show_sn_labels, show_energy_labels=show_energy_labels)
    
    if show_sn:
        _grouped_boxplot(
            long_filt,
            group_col="SN",
            energies_order=energies,
            show_dates=show_dates,
            show_energy_labels=show_energy_labels
        )

    
    if show_energy:
        _grouped_boxplot(
            long_filt,
            group_col="Energy",
            energies_order=energies,
            show_dates=show_dates,
            show_sn_labels=show_sn_labels,
            show_energy_labels=show_energy_labels
        )




# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Output Round Robin — Boxplots")
        self.geometry("720x520")  # taller to fit 3 listboxes
        pad = {"padx": 8, "pady": 6}

        ttk.Label(self, text="Data file (XLSX/CSV):").grid(row=0, column=0, sticky="w", **pad)
        self.path_var = tk.StringVar(value=DEFAULT_PATH)
        ttk.Entry(self, textvariable=self.path_var, width=70).grid(row=0, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.browse).grid(row=0, column=2, **pad)

        # ---- Plot selection checkboxes ----
        self.do_system = tk.BooleanVar(value=True)
        self.do_sn = tk.BooleanVar(value=True)
        self.do_energy = tk.BooleanVar(value=True)
        self.show_dates = tk.BooleanVar(value=False)
        self.show_sn_labels = tk.BooleanVar(value=False)
        self.show_energy_labels = tk.BooleanVar(value=False)

        ttk.Checkbutton(
            self, text="System boxplot (all data)",
            variable=self.do_system
        ).grid(row=1, column=0, columnspan=3, sticky="w", **pad)

        ttk.Checkbutton(
            self, text="Grouped boxplot by SN",
            variable=self.do_sn
        ).grid(row=2, column=0, columnspan=3, sticky="w", **pad)

        ttk.Checkbutton(
            self, text="Grouped boxplot by Energy",
            variable=self.do_energy
        ).grid(row=3, column=0, columnspan=3, sticky="w", **pad)
        
        ttk.Checkbutton(
            self, text="Show Dates",
            variable=self.show_dates
        ).grid(row=4, column=0, columnspan=3, sticky="w", **pad)
        ttk.Checkbutton(self, text="Show SN labels (System + Energy plots)", variable=self.show_sn_labels).grid(row=5, column=0, columnspan=3, sticky="w", **pad)
        ttk.Checkbutton(self, text="Show Energy labels (System plot only)", variable=self.show_energy_labels).grid(row=6, column=0, columnspan=3, sticky="w", **pad)

        # ---- System selection ----
        ttk.Label(self, text="Select System(s) for global filter:").grid(
            row=7, column=0, columnspan=3, sticky="w", **pad
        )
        self.system_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.system_listbox.grid(row=8, column=0, columnspan=3, sticky="nsew", **pad)

        # ---- SN selection ----
        ttk.Label(self, text="Select SN(s) for SN-grouped plot:").grid(
            row=9, column=0, columnspan=3, sticky="w", **pad
        )
        self.sn_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.sn_listbox.grid(row=10, column=0, columnspan=3, sticky="nsew", **pad)

        # ---- Energy selection ----
        ttk.Label(self, text="Select Energy(ies) for Energy-grouped plot:").grid(
            row=11, column=0, columnspan=3, sticky="w", **pad
        )
        self.energy_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.energy_listbox.grid(row=12, column=0, columnspan=3, sticky="nsew", **pad)

        # backing lists for listbox indices
        self.system_values = []
        self.sn_values = []
        self.energy_values = []

        # Plot button
        ttk.Button(self, text="Plot", command=self.plot).grid(row=13, column=0, columnspan=3, **pad)

        # populate lists from default file if possible
        self.populate_lists_from_file()

    def populate_lists_from_file(self):
        """Read the current file and populate System/SN/Energy listboxes if possible."""
        path = self.path_var.get().strip()
        if not path:
            return
        if not os.path.exists(path):
            return

        try:
            df = read_table(path)
            long, energies = to_long(df)
        except Exception:
            # If it fails (missing file, bad format, etc.), just skip; Plot will show error later
            return

        # ---- System list ----
        self.system_listbox.delete(0, tk.END)
        self.system_values = []

        sys_series = long["System"].astype(str).str.strip()
        sys_vals_raw = pd.unique(sys_series)
        sys_vals_raw = [s for s in sys_vals_raw if s not in (None, "", "nan")]

        sys_set = set(sys_vals_raw)
        # order known systems according to SYSTEM_ORDER, then any extras alphabetically
        ordered_known = [s for s in SYSTEM_ORDER if s in sys_set]
        extras = sorted(sys_set - set(SYSTEM_ORDER))
        systems = ordered_known + extras

        self.system_values = systems
        for s in systems:
            self.system_listbox.insert(tk.END, s)
        if systems:
            self.system_listbox.selection_set(0, tk.END)

        # ---- SN list ----
        self.sn_listbox.delete(0, tk.END)
        self.sn_values = []

        sn_vals = sorted([
            g for g in pd.unique(long["SN"])
            if not (isinstance(g, float) and np.isnan(g))
        ])
        self.sn_values = sn_vals
        for sn in sn_vals:
            try:
                label = f"SN {int(sn)}"
            except (ValueError, TypeError):
                label = str(sn)
            self.sn_listbox.insert(tk.END, label)
        if sn_vals:
            self.sn_listbox.selection_set(0, tk.END)

        # ---- Energy list ----
        self.energy_listbox.delete(0, tk.END)
        self.energy_values = energies[:]  # copy
        for e in self.energy_values:
            self.energy_listbox.insert(tk.END, e)
        if self.energy_values:
            self.energy_listbox.selection_set(0, tk.END)

    def browse(self):
        p = filedialog.askopenfilename(filetypes=[
            ("Excel/CSV", "*.xlsx *.xls *.csv *.txt"),
            ("Excel", "*.xlsx *.xls"), ("CSV", "*.csv *.txt"),
            ("All files", "*.*"),
        ])
        if p:
            self.path_var.set(p)
            self.populate_lists_from_file()

    def plot(self):
        path = self.path_var.get().strip()
        if not path:
            messagebox.showerror("Error", "Pick a data file.")
            return
        try:
            df = read_table(path)

            # ---- System filter (global) ----
            system_filter = None
            if self.system_values:
                idxs = self.system_listbox.curselection()
                if idxs:
                    system_filter = [self.system_values[i] for i in idxs]

            # ---- SN filter (global) ----
            sn_filter = None
            if self.sn_values:
                idxs = self.sn_listbox.curselection()
                if idxs:
                    sn_filter = [self.sn_values[i] for i in idxs]

            # ---- Energy filter (global) ----
            energy_filter = None
            if self.energy_values:
                idxs = self.energy_listbox.curselection()
                if idxs:
                    energy_filter = [self.energy_values[i] for i in idxs]

            make_plots(df,
               show_system=self.do_system.get(),
               show_sn=self.do_sn.get(),
               show_energy=self.do_energy.get(),
               sn_filter=sn_filter,
               energy_filter=energy_filter,
               system_filter=system_filter,
               show_dates=self.show_dates.get(),
               show_sn_labels=self.show_sn_labels.get(),
               show_energy_labels=self.show_energy_labels.get())
        except Exception as e:
            messagebox.showerror("Plot error", str(e))


if __name__ == "__main__":
    App().mainloop()
