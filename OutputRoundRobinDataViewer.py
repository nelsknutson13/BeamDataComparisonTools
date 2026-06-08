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
from matplotlib.ticker import MultipleLocator
import statsmodels.formula.api as smf

# ---------- Config ----------
DEFAULT_PATH = r"C:\Users\nknutson\OneDrive - Washington University in St. Louis\NGDS QA Consortium\OutputRoundRobinData.xlsx"
EXPECTED_ENERGIES = ["6X", "10X", "15X", "6FFF", "8FFF", "10FFF"]
SYSTEM_ORDER = ["Institution", "IROC", "RDS", "Consortium Audit", "Institution Reported (Consortium Form)", "Consortium Audit (Form-Corrected)"]
ORDER = {name: i for i, name in enumerate(SYSTEM_ORDER)}
# Marker per system (used in both boxplots)
MARKER_MAP = {
    "Institution": "P",
    "IROC": "D",
    "RDS": "o",
    "Institution Reported (Consortium Form)": "s",
    "Consortium Audit": "X",
    "Consortium Audit (Form-Corrected)": "^",
}
OUTLIER_MIN_N = 15   # only show 1.5*IQR fliers when N >= this

MARKER_COLOR = "black"   # or "#444", "0.3", etc.

DEFAULT_TOL = 0.05   # default acceptance tolerance (±5%), drawn as dashed lines around the 1.0 baseline

# Default y-axis range and tick spacing, all in PERCENT around the 1.0 baseline.
DEFAULT_YRANGE_PCT = 7.0    # y-limits = 1 ± 7%
DEFAULT_MAJOR_PCT  = 1.0    # major tick every 1%
DEFAULT_MINOR_PCT  = 0.5    # minor tick every 0.5%

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


def renormalize_to_reference(long: pd.DataFrame, reference_system: str):
    """Re-reference Ratio values so `reference_system` sits at 1.0 within each
    Institution-anchored session. Returns (new_long, n_dropped).

    Sessions: within each (SN, Energy) sorted by date, every 'Institution' row
    starts a new session; subsequent non-institution rows attach to the most
    recent preceding Institution row. Each session is divided by its own
    reference value (no cross-session averaging). If the reference system
    appears multiple times in a session, those are averaged (with a warning).

    'Institution' is already the 1.0 baseline → no-op.
    """
    ref = (reference_system or "Institution").strip()
    if ref == "Institution":
        return long, 0

    data = long.copy()
    data["System"] = data["System"].astype(str).str.strip()
    data["_date"] = pd.to_datetime(data["Date"], errors="coerce")

    kept = []
    dropped = orphan = 0
    multi_ref_warned = False

    for (sn, en), grp in data.groupby(["SN", "Energy"], dropna=False):
        grp = grp.sort_values("_date", kind="mergesort")
        # Assign session ids: each Institution row increments the counter.
        sids, sid = [], 0
        for sysname in grp["System"]:
            if sysname == "Institution":
                sid += 1
            sids.append(sid)
        grp = grp.assign(_session=sids)

        # Rows before the first Institution (session 0) are orphans.
        orphan += int((grp["_session"] == 0).sum())
        grp = grp[grp["_session"] > 0]

        for _sess, sgrp in grp.groupby("_session"):
            ref_vals = sgrp.loc[sgrp["System"] == ref, "Ratio"].dropna()
            if ref_vals.empty:
                dropped += len(sgrp)
                continue
            if len(ref_vals) > 1 and not multi_ref_warned:
                print(f"[renormalize] reference '{ref}' has multiple measurements in a "
                      f"session (e.g. SN {sn}, {en}) — averaging them for the divisor.")
                multi_ref_warned = True
            ref_val = float(ref_vals.mean())
            if ref_val == 0 or not np.isfinite(ref_val):
                dropped += len(sgrp)
                continue
            s = sgrp.copy()
            s["Ratio"] = s["Ratio"] / ref_val
            kept.append(s)

    if orphan:
        print(f"[renormalize] dropped {orphan} orphan row(s) appearing before any "
              f"Institution measurement.")
    if dropped:
        print(f"[renormalize] dropped {dropped} row(s) in sessions with no '{ref}' measurement.")

    if not kept:
        return data.iloc[0:0].drop(columns=["_date", "_session"], errors="ignore"), dropped + orphan
    out = pd.concat(kept, ignore_index=True).drop(columns=["_date", "_session"], errors="ignore")
    return out, dropped + orphan


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


def _draw_tolerance(ax, tol: float = DEFAULT_TOL):
    """Dashed lines at 1±tol marking the acceptance tolerance band."""
    for y in (1.0 - tol, 1.0 + tol):
        ax.axhline(y, linestyle="--", color="0.5", linewidth=0.9, zorder=1)


def _tolerance_handle(tol: float = DEFAULT_TOL):
    return Line2D([0], [0], linestyle="--", color="0.5",
                  label=f"Tolerance ±{tol * 100:g}%")


def _apply_yaxis(ax, range_pct=DEFAULT_YRANGE_PCT,
                 major_pct=DEFAULT_MAJOR_PCT, minor_pct=DEFAULT_MINOR_PCT):
    """Set y-limits to 1±range_pct and major/minor tick spacing — all in percent
    of the 1.0 baseline. A non-positive range or spacing skips that piece."""
    if range_pct and range_pct > 0:
        r = range_pct / 100.0
        ax.set_ylim(1.0 - r, 1.0 + r)
    if major_pct and major_pct > 0:
        ax.yaxis.set_major_locator(MultipleLocator(major_pct / 100.0))
    if minor_pct and minor_pct > 0:
        ax.yaxis.set_minor_locator(MultipleLocator(minor_pct / 100.0))


def _system_boxplot(long: pd.DataFrame, show_dates: bool = False, show_sn_labels: bool = False, show_energy_labels: bool = False, ref_label: str = "Institution", show_tolerance: bool = True, tol: float = DEFAULT_TOL,
                    y_range: float = DEFAULT_YRANGE_PCT, major_tick: float = DEFAULT_MAJOR_PCT, minor_tick: float = DEFAULT_MINOR_PCT,
                    show_points: bool = True):


    data = long.copy()
    data["System"] = data["System"].astype(str).str.strip()

    # Keep only systems in your fixed order and only those present; exclude the
    # reference system (it's the degenerate 1.0 baseline, shown as the dashed line).
    present_set = set(data["System"].dropna())
    systems_present = [s for s in SYSTEM_ORDER if s in present_set and s != ref_label]
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
    if show_tolerance:
        _draw_tolerance(ax, tol)

    # Overlay sample points with shared marker map
    marker_map = MARKER_MAP
    for sys in (systems_present if show_points else []):
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
    ax.set_ylabel("Output [cGy/MU]" if ref_label == "Institution"
                  else f"Output (relative to {ref_label})")
    _apply_yaxis(ax, y_range, major_tick, minor_tick)

    # Legend to match grouped plot
    handles = [
        Line2D([], [], linestyle="None",
               marker=marker_map.get(sys, "D"), markersize=10,
               markerfacecolor="none", markeredgecolor=MARKER_COLOR, color=MARKER_COLOR,
               label=sys)
        for sys in systems_present
    ]
    handles.append(Line2D([0], [0], linestyle="--", color="black", label=ref_label))
    handles.append(Line2D([], [], linestyle="None", marker=r'$\ast$', markersize=7,
                          color=MARKER_COLOR, label="Outlier"))
    if show_tolerance:
        handles.append(_tolerance_handle(tol))
    ax.legend(handles=handles, loc="best", frameon=True)

    plt.tight_layout()
    plt.show()


# ---------- Plot ----------
def _grouped_boxplot(long: pd.DataFrame, group_col: str, energies_order,
                     show_dates: bool = False,
                     show_sn_labels: bool = False,
                     show_energy_labels: bool = False,
                     ref_label: str = "Institution",
                     show_tolerance: bool = True,
                     tol: float = DEFAULT_TOL,
                     y_range: float = DEFAULT_YRANGE_PCT,
                     major_tick: float = DEFAULT_MAJOR_PCT,
                     minor_tick: float = DEFAULT_MINOR_PCT,
                     show_points: bool = True):


    data = long.copy()
    data["System"] = data["System"].astype(str).str.strip()
    data = data[data["System"].isin(ORDER.keys())]  # keep only your ordered systems

    # Fixed system order; include only those present, excluding the reference
    # system (degenerate 1.0 baseline, shown as the dashed line).
    present = set(data["System"].astype(str).str.strip().unique())
    systems = [s for s in sorted(present, key=lambda s: ORDER.get(s, 999)) if s != ref_label]

    if not systems:
        raise ValueError("No systems to plot after excluding the reference system.")

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
    if show_tolerance:
        _draw_tolerance(ax, tol)

    # Overlay sample points with shared marker map
    marker_map = MARKER_MAP

    for g in (groups if show_points else []):
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
    ax.set_ylabel("Output [cGy/MU]" if ref_label == "Institution"
                  else f"Output (relative to {ref_label})")
    _apply_yaxis(ax, y_range, major_tick, minor_tick)

    # Legend: markers for each system + dashed baseline
    handles = [Line2D([], [], linestyle="None",
                      marker=marker_map.get(sys, "D"), markersize=10,
                      markerfacecolor="none", markeredgecolor=MARKER_COLOR, color=MARKER_COLOR,
                      label=sys) for sys in systems]

    handles.append(Line2D([0], [0], linestyle="--", color="black", label=ref_label))
    handles.append(Line2D([], [], linestyle="None", marker=r'$\ast$', markersize=7,
                          color=MARKER_COLOR, label="Outlier"))
    if show_tolerance:
        handles.append(_tolerance_handle(tol))

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
        # statsmodels rounds to 3 decimals, so tiny variances print as 0.000.
        # On the variance-component rows only (they contain 'Var'), show
        # '<0.001' instead — and replace '0.000 ' (incl. one trailing space) so
        # the column stays the same width. Coefficient rows are left untouched.
        summary_lines = res.summary().as_text().split("\n")
        for i, ln in enumerate(summary_lines):
            if "Var" in ln and "0.000 " in ln:
                summary_lines[i] = ln.replace("0.000 ", "<0.001", 1)
        print("\n".join(summary_lines))

        # The summary table rounds variances to 3 decimals, so small (but
        # nonzero) random-effect variances print as 0.000. Re-report them at
        # full precision and as standard deviations in % (sqrt of variance,
        # ×100 since Delta is in fraction units) — far easier to interpret.
        def _fmt_var(label, var):
            var = float(var)
            sd_pct = np.sqrt(var) * 100 if var > 0 else 0.0
            print(f"  {label:<22} var={var:.3e}   SD={sd_pct:.4f}%")

        print("\nVariance components (full precision):")
        _fmt_var("Residual (Scale)", res.scale)
        try:
            _fmt_var("Group (SN)", res.cov_re.iloc[0, 0])
        except Exception as e:
            print(f"  Group (SN) var unavailable: {e}")
        try:
            for name, v in zip(vc.keys(), np.atleast_1d(res.vcomp)):
                _fmt_var(name, v)
        except Exception as e:
            print(f"  Variance components unavailable: {e}")

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
               show_energy_labels: bool = False,
               normalize_to: str = "Institution",
               show_tolerance: bool = True,
               tolerance: float = DEFAULT_TOL,
               y_range: float = DEFAULT_YRANGE_PCT,
               major_tick: float = DEFAULT_MAJOR_PCT,
               minor_tick: float = DEFAULT_MINOR_PCT,
               show_points: bool = True):

    long, energies = to_long(df)

    # Re-reference per Institution-anchored session BEFORE filtering, so the
    # session anchors (Institution rows) are still present even if the user
    # later filters out a system.
    long, _ = renormalize_to_reference(long, normalize_to)
    if long.empty:
        raise ValueError(f"No data left after normalizing to '{normalize_to}'.")
    ref_label = (normalize_to or "Institution").strip()

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
        _system_boxplot(long_filt, show_dates=show_dates, show_sn_labels=show_sn_labels,
                        show_energy_labels=show_energy_labels, ref_label=ref_label,
                        show_tolerance=show_tolerance, tol=tolerance,
                        y_range=y_range, major_tick=major_tick, minor_tick=minor_tick,
                        show_points=show_points)

    if show_sn:
        _grouped_boxplot(
            long_filt,
            group_col="SN",
            energies_order=energies,
            show_dates=show_dates,
            show_energy_labels=show_energy_labels,
            ref_label=ref_label,
            show_tolerance=show_tolerance,
            tol=tolerance,
            y_range=y_range, major_tick=major_tick, minor_tick=minor_tick,
            show_points=show_points
        )


    if show_energy:
        _grouped_boxplot(
            long_filt,
            group_col="Energy",
            energies_order=energies,
            show_dates=show_dates,
            show_sn_labels=show_sn_labels,
            show_energy_labels=show_energy_labels,
            ref_label=ref_label,
            show_tolerance=show_tolerance,
            tol=tolerance,
            y_range=y_range, major_tick=major_tick, minor_tick=minor_tick,
            show_points=show_points
        )




# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Output Round Robin — Boxplots")
        self.geometry("720x650")  # taller to fit 3 listboxes + normalize dropdown + tolerance + axis controls
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
        self.show_tolerance = tk.BooleanVar(value=True)
        self.tol_var = tk.StringVar(value="5")
        self.y_range_var = tk.StringVar(value="7")
        self.major_tick_var = tk.StringVar(value="2")
        self.minor_tick_var = tk.StringVar(value="1")
        self.show_points = tk.BooleanVar(value=True)

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
        ttk.Checkbutton(self, text="Show tolerance  ±", variable=self.show_tolerance).grid(row=7, column=0, sticky="e", **pad)
        ttk.Entry(self, textvariable=self.tol_var, width=6).grid(row=7, column=1, sticky="w", **pad)
        ttk.Label(self, text="%").grid(row=7, column=1, sticky="w", padx=(60, 0))
        ttk.Checkbutton(self, text="Show all points (off = outliers only)", variable=self.show_points).grid(row=7, column=2, sticky="w", **pad)

        # ---- System selection ----
        ttk.Label(self, text="Select System(s) for global filter:").grid(
            row=8, column=0, columnspan=3, sticky="w", **pad
        )
        self.system_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.system_listbox.grid(row=9, column=0, columnspan=3, sticky="nsew", **pad)

        # ---- SN selection ----
        ttk.Label(self, text="Select SN(s) for SN-grouped plot:").grid(
            row=10, column=0, columnspan=3, sticky="w", **pad
        )
        self.sn_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.sn_listbox.grid(row=11, column=0, columnspan=3, sticky="nsew", **pad)

        # ---- Energy selection ----
        ttk.Label(self, text="Select Energy(ies) for Energy-grouped plot:").grid(
            row=12, column=0, columnspan=3, sticky="w", **pad
        )
        self.energy_listbox = tk.Listbox(self, selectmode="extended", height=4, exportselection=False)
        self.energy_listbox.grid(row=13, column=0, columnspan=3, sticky="nsew", **pad)

        # backing lists for listbox indices
        self.system_values = []
        self.sn_values = []
        self.energy_values = []

        # ---- Normalize-to selection ----
        ttk.Label(self, text="Normalize to:").grid(row=14, column=0, sticky="w", **pad)
        self.normalize_var = tk.StringVar(value="Institution")
        self.normalize_combo = ttk.Combobox(self, textvariable=self.normalize_var,
                                             state="readonly", width=40)
        self.normalize_combo["values"] = ["Institution"]
        self.normalize_combo.grid(row=14, column=1, sticky="w", **pad)

        # ---- Y-axis range / tick spacing (all in percent) ----
        axis_frame = ttk.Frame(self)
        axis_frame.grid(row=15, column=0, columnspan=3, sticky="w", **pad)
        ttk.Label(axis_frame, text="Y-axis ±(%):").pack(side="left")
        ttk.Entry(axis_frame, textvariable=self.y_range_var, width=5).pack(side="left", padx=(2, 12))
        ttk.Label(axis_frame, text="Major tick (%):").pack(side="left")
        ttk.Entry(axis_frame, textvariable=self.major_tick_var, width=5).pack(side="left", padx=(2, 12))
        ttk.Label(axis_frame, text="Minor tick (%):").pack(side="left")
        ttk.Entry(axis_frame, textvariable=self.minor_tick_var, width=5).pack(side="left", padx=(2, 0))

        # Plot button
        ttk.Button(self, text="Plot", command=self.plot).grid(row=16, column=0, columnspan=3, **pad)

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

        # ---- Normalize-to dropdown options (Institution first, then others) ----
        norm_opts = (["Institution"] if "Institution" in sys_set else []) + \
                    [s for s in systems if s != "Institution"]
        if not norm_opts:
            norm_opts = ["Institution"]
        self.normalize_combo["values"] = norm_opts
        if self.normalize_var.get() not in norm_opts:
            self.normalize_var.set(norm_opts[0])

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

            # ---- Tolerance (% -> fraction) ----
            tolerance = DEFAULT_TOL
            if self.show_tolerance.get():
                try:
                    tolerance = float(self.tol_var.get()) / 100.0
                except ValueError:
                    messagebox.showerror("Error", "Tolerance must be a number (percent).")
                    return

            # ---- Y-axis range / tick spacing (percent) ----
            try:
                y_range    = float(self.y_range_var.get())
                major_tick = float(self.major_tick_var.get())
                minor_tick = float(self.minor_tick_var.get())
            except ValueError:
                messagebox.showerror("Error", "Y-axis range and tick spacing must be numbers (percent).")
                return

            make_plots(df,
               show_system=self.do_system.get(),
               show_sn=self.do_sn.get(),
               show_energy=self.do_energy.get(),
               sn_filter=sn_filter,
               energy_filter=energy_filter,
               system_filter=system_filter,
               show_dates=self.show_dates.get(),
               show_sn_labels=self.show_sn_labels.get(),
               show_energy_labels=self.show_energy_labels.get(),
               normalize_to=self.normalize_var.get(),
               show_tolerance=self.show_tolerance.get(),
               tolerance=tolerance,
               y_range=y_range,
               major_tick=major_tick,
               minor_tick=minor_tick,
               show_points=self.show_points.get())
        except Exception as e:
            messagebox.showerror("Plot error", str(e))


if __name__ == "__main__":
    App().mainloop()
