"""
Batch PDD Composite Analysis
-----------------------------
Loops over all energy folders, runs composite analysis on
Measurement vs MC sheets, saves a figure per energy and a
summary CSV of pass-rate statistics.

Parameters are set in the CONFIG section below.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend (no GUI needed)
from matplotlib import pyplot as plt
from scipy import interpolate as interp
from scipy.ndimage import gaussian_filter1d

from comp import dosedif, dta as dtafunc

# ─────────────────────────────────────────────
#  CONFIG  — edit these as needed
# ─────────────────────────────────────────────
BASE_PATH = (
    r"C:\Users\nknutson\Box\Knutson\Research\Projects underway"
    r"\Radformation research\RADMC Photon VALIDATION\TrueBeam\ProcessedData"
)

SHEET1_NAME = "Measurement"   # reference / measured
SHEET2_NAME = "MC"            # dataset to convolve

ANALYSIS      = 'plot'        # 'comp', 'dif', 'dist', or 'plot'
DD_CRITERIA   = 2.0           # dose difference threshold [%]
DTA_CRITERIA  = 0.2           # DTA threshold [cm]  (2 mm)
NORM          = 2             # 1 = normalize to dmax, 2 = normalize to fixed depth per energy

# Fixed normalization depth [cm] per energy folder name
NORM_DEPTH = {
    '6X':   1.3,
    '6FFF': 1.3,
    '10X':  2.2,
    '10FFF':2.2,
    '15X':  3.5,
}
DEPTH_SHIFT1  = 0.0           # depth shift for sheet 1 [cm]
DEPTH_SHIFT2  = 0.0           # depth shift for sheet 2 [cm]
CUTOFF_DEPTH  = 0.0           # discard points shallower than this [cm]; 0 = no cutoff
CONV_FWHM_CM  = 0.48          # PTW 31021: 2 × 2.4 mm cavity radius = 4.8 mm = 0.48 cm
MARKER_SIZE   = 2             # plot marker size
DPI           = 400           # figure output resolution

RESULTS_DIR = os.path.join(BASE_PATH, "Results")
# ─────────────────────────────────────────────


# ── helpers copied from PDDCompare.py ────────

def apply_detector_convolution(profile_y, profile_dose, fwhm_cm):
    """Gaussian blur along depth axis to simulate detector volume averaging."""
    if fwhm_cm <= 0:
        return profile_dose
    sigma_cm  = fwhm_cm / 2.355
    step_cm   = np.mean(np.diff(np.asarray(profile_y, dtype=float)))
    sigma_idx = sigma_cm / step_cm
    dose_arr  = np.asarray(profile_dose, dtype=float)
    return gaussian_filter1d(dose_arr, sigma=sigma_idx)


def downsample_to_native(x_interp, vals_interp, x_native):
    """Thin interpolated output back to roughly native data spacing."""
    x_interp    = np.asarray(x_interp,    dtype=float)
    vals_interp = np.asarray(vals_interp, dtype=float)
    x_native    = np.asarray(x_native,    dtype=float)

    if x_interp.size < 2 or x_native.size < 2:
        return x_interp, vals_interp

    xu = np.unique(x_native[np.isfinite(x_native)])
    xu.sort()
    diffs = np.diff(xu)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return x_interp, vals_interp

    native_step = float(np.median(diffs))
    xi = np.unique(x_interp[np.isfinite(x_interp)])
    xi.sort()
    interp_step = float(np.median(np.diff(xi))) if xi.size >= 2 else 0.0
    if interp_step <= 0:
        return x_interp, vals_interp

    factor = max(1, int(round(native_step / interp_step)))
    return x_interp[::factor], vals_interp[::factor]


# ── analysis for one file ─────────────────────

def run_one_file(xlsx_path, energy_label):
    """
    Run composite analysis on one PDD Excel file.
    Returns a list of result dicts (one per FS).
    """
    print(f"\n{'='*60}")
    print(f"  {energy_label}  —  {os.path.basename(xlsx_path)}")
    print(f"{'='*60}")

    df1 = pd.read_excel(xlsx_path, sheet_name=SHEET1_NAME, header=0)
    df2 = pd.read_excel(xlsx_path, sheet_name=SHEET2_NAME, header=0)
    df1 = df1.sort_values(by=['FS', 'Axis', 'Pos'])
    df2 = df2.sort_values(by=['FS', 'Axis', 'Pos'])

    # FS values common to both sheets (Z-axis rows only)
    fs1 = set(df1.loc[df1['Axis'] == 'Z', 'FS'].unique())
    fs2 = set(df2.loc[df2['Axis'] == 'Z', 'FS'].unique())
    fsl = sorted(fs1 & fs2)

    if not fsl:
        print("  No common field sizes found — skipping.")
        return []

    # ── figure setup ─────────────────────────
    if ANALYSIS == 'comp':
        fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 11),
            gridspec_kw={'height_ratios': [1.5, 1, 1]})
        plt.rcParams.update({'font.size': 20})
        ax1.set_ylabel('DTA [mm]')
        ax2.set_ylabel('Dose Difference [%]')
    elif ANALYSIS in ('dif', 'dist'):
        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(15, 9),
            gridspec_kw={'height_ratios': [1.5, 1]})
        plt.rcParams.update({'font.size': 18})
        ax1.set_ylabel('Dose Difference [%]' if ANALYSIS == 'dif' else 'DTA [mm]')
    else:  # 'plot'
        fig, ax0 = plt.subplots(1, 1, figsize=(15, 8))
        plt.rcParams.update({'font.size': 16})

    ax0.plot([], '+r', ms=10, label=SHEET1_NAME)
    ax0.plot([], '.k', ms=10, label=SHEET2_NAME)
    ax0.legend()
    ax0.set_ylabel('Percentage Depth Dose [%]')
    ax0.set_xlabel('Depth [cm]')

    results = []
    total_points_cumulative = 0
    total_fails_cumulative  = 0

    for fs_val in fsl:
        # ── extract data ───────────────────
        mask1 = (df1['FS'] == fs_val) & (df1['Axis'] == 'Z')
        mask2 = (df2['FS'] == fs_val) & (df2['Axis'] == 'Z')
        y1 = df1.loc[mask1, 'Pos'].reset_index(drop=True)
        d1 = df1.loc[mask1, 'Dose'].reset_index(drop=True)
        y2 = df2.loc[mask2, 'Pos'].reset_index(drop=True)
        d2 = df2.loc[mask2, 'Dose'].reset_index(drop=True)

        if len(y1) < 2 or len(y2) < 2:
            print(f"  FS {fs_val}: not enough points — skipping.")
            continue

        # ── depth shifts ───────────────────
        y1 = y1 + DEPTH_SHIFT1
        y2 = y2 + DEPTH_SHIFT2

        # ── cutoff ─────────────────────────
        if CUTOFF_DEPTH > 0:
            m1 = y1 >= CUTOFF_DEPTH
            m2 = y2 >= CUTOFF_DEPTH
            y1, d1 = y1[m1].reset_index(drop=True), d1[m1].reset_index(drop=True)
            y2, d2 = y2[m2].reset_index(drop=True), d2[m2].reset_index(drop=True)
            if len(y1) < 2 or len(y2) < 2:
                print(f"  FS {fs_val}: cutoff removed too many points — skipping.")
                continue

        # ── normalization ──────────────────
        if NORM == 1:
            d1 = d1 / d1.max() * 100
            d2 = d2 / d2.max() * 100
        elif NORM == 2:
            fixed_depth = NORM_DEPTH.get(energy_label, None)
            if fixed_depth is None:
                print(f"  WARNING: no fixed depth for '{energy_label}', falling back to max norm.")
                d1 = d1 / d1.max() * 100
                d2 = d2 / d2.max() * 100
            else:
                d1_fixed = float(interp.pchip(y1, d1)(fixed_depth))
                d2_fixed = float(interp.pchip(y2, d2)(fixed_depth))
                d1 = d1 / d1_fixed * 100
                d2 = d2 / d2_fixed * 100

        # ── detector convolution (MC = curve 2) ──
        d2 = apply_detector_convolution(y2, d2, CONV_FWHM_CM)

        # ── plot curves ───────────────────
        ax0.plot(y1, d1, '+r', ms=MARKER_SIZE)
        ax0.plot(y2, d2, '.k', ms=MARKER_SIZE)

        # ── analysis metrics & subplot plots ──
        xlim = (max(float(y1.min()), float(y2.min())),
                min(float(y1.max()), float(y2.max())))
        ax0.set_xlim(*xlim)

        if ANALYSIS == 'plot':
            ax0.set_ylim(0, max(1.05 * np.max(d1), 1.05 * np.max(d2)))

        elif ANALYSIS == 'comp':
            dtax, dtav = dtafunc(y1, d1, y2, d2, DTA_CRITERIA)
            difx, difv = dosedif(y1, d1, y2, d2, 0)
            dtax, dtav = downsample_to_native(dtax, dtav, y1)
            difx, difv = downsample_to_native(difx, difv, y1)
            dtafail  = np.where(np.abs(dtav) > DTA_CRITERIA)[0]
            ddiffail = np.where(np.abs(difv) > DD_CRITERIA)[0]
            fail_set = set(dtafail) & set(ddiffail)
            totalfail = len(fail_set)
            total     = len(dtav)
            mean_dd  = float(np.mean(np.abs(difv))) if total > 0 else 0.0
            mean_dta = float(np.mean(np.abs(dtav))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDoseDiff_pct': round(mean_dd, 3),
                'MeanDTA_mm': round(mean_dta * 10, 3), 'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDD={mean_dd:.2f}%  MeanDTA={mean_dta*10:.2f}mm  Fail={totalfail}/{total}")
            for ax in (ax1, ax2):
                ax.set_xlim(*xlim)
            ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
            ax2.plot(difx, difv,      '.g', ms=0.5)
            for n in fail_set:
                ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)
                ax2.plot(difx[n], difv[n],      '.r', ms=0.5)

        elif ANALYSIS == 'dif':
            difx, difv = dosedif(y1, d1, y2, d2, 0)
            difx, difv = downsample_to_native(difx, difv, y1)
            ddiffail = np.where(np.abs(difv) > DD_CRITERIA)[0]
            totalfail = len(ddiffail)
            total     = len(difv)
            mean_dd  = float(np.mean(np.abs(difv))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDoseDiff_pct': round(mean_dd, 3),
                'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDD={mean_dd:.2f}%  Fail={totalfail}/{total}")
            ax1.set_xlim(*xlim)
            ax1.plot(difx, difv, '.g', ms=0.5)
            for n in ddiffail:
                ax1.plot(difx[n], difv[n], '.r', ms=0.5)

        elif ANALYSIS == 'dist':
            dtax, dtav = dtafunc(y1, d1, y2, d2, DTA_CRITERIA)
            dtax, dtav = downsample_to_native(dtax, dtav, y1)
            dtafail  = np.where(np.abs(dtav) > DTA_CRITERIA)[0]
            totalfail = len(dtafail)
            total     = len(dtav)
            mean_dta = float(np.mean(np.abs(dtav))) if total > 0 else 0.0
            passrate = (total - totalfail) / total * 100 if total > 0 else 0.0
            total_points_cumulative += total
            total_fails_cumulative  += totalfail
            results.append({'Energy': energy_label, 'FS': fs_val,
                'PassRate_pct': round(passrate, 2), 'MeanDTA_mm': round(mean_dta * 10, 3),
                'TotalPoints': total, 'FailPoints': totalfail})
            print(f"  FS {fs_val:.1f} cm : Pass={passrate:.2f}%  "
                  f"MeanDTA={mean_dta*10:.2f}mm  Fail={totalfail}/{total}")
            ax1.set_xlim(*xlim)
            ax1.plot(dtax, dtav * 10, '.g', ms=0.5)
            for n in dtafail:
                ax1.plot(dtax[n], dtav[n] * 10, '.r', ms=0.5)

    # ── figure summary label & title ──────────
    if total_points_cumulative > 0:
        overall_pass = 100 * (1 - total_fails_cumulative / total_points_cumulative)
    else:
        overall_pass = 0.0

    if ANALYSIS == 'comp':
        ax2.set_xlabel(f'Points outside {DD_CRITERIA:.1f}% & {DTA_CRITERIA*10:.1f}mm  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'Composite {DD_CRITERIA:.0f}%/{DTA_CRITERIA*10:.0f}mm'
    elif ANALYSIS == 'dif':
        ax1.set_xlabel(f'Points outside {DD_CRITERIA:.1f}%  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'Dose Difference {DD_CRITERIA:.0f}%'
    elif ANALYSIS == 'dist':
        ax1.set_xlabel(f'Points outside {DTA_CRITERIA*10:.1f}mm  '
                       f'{total_fails_cumulative}/{total_points_cumulative}  '
                       f'Pass Rate: {overall_pass:.2f}%')
        title_tag = f'DTA {DTA_CRITERIA*10:.0f}mm'
    else:
        title_tag = 'Plots Only'

    fig.suptitle(f'{energy_label} — {SHEET1_NAME} vs {SHEET2_NAME}  {title_tag}', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.subplots_adjust(hspace=0.45)

    # ── save figure ───────────────────────────
    stem     = os.path.splitext(os.path.basename(xlsx_path))[0]
    safe_tag = title_tag.replace(' ', '_').replace('/', '-').replace('%', 'pct')
    out_png  = os.path.join(RESULTS_DIR, f"{energy_label}_{stem}_{safe_tag}.png")
    fig.savefig(out_png, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"  Figure saved: {out_png}")

    return results


# ── main ─────────────────────────────────────

def main():
    from datetime import datetime
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_csv = os.path.join(RESULTS_DIR, f"pdd_comparison_summary_{timestamp}.csv")
    all_results = []

    energy_dirs = sorted([
        d for d in os.listdir(BASE_PATH)
        if os.path.isdir(os.path.join(BASE_PATH, d))
        and os.path.join(BASE_PATH, d) != RESULTS_DIR
    ])

    for energy_label in energy_dirs:
        folder = os.path.join(BASE_PATH, energy_label)
        pdd_files = glob.glob(os.path.join(folder, '*PDD*.xlsx'))
        if not pdd_files:
            print(f"No PDD xlsx found in {folder} — skipping.")
            continue

        xlsx_path = pdd_files[0]
        results = run_one_file(xlsx_path, energy_label)
        all_results.extend(results)

    # ── save summary CSV ─────────────────────
    if all_results:
        df_out = pd.DataFrame(all_results, columns=[
            'Energy', 'FS', 'PassRate_pct', 'MeanDoseDiff_pct',
            'MeanDTA_mm', 'TotalPoints', 'FailPoints'
        ])

        # ── build pivot table (pass rate % by energy × FS) ───────────
        energy_order = [e for e in ['6X', '6FFF', '10X', '10FFF', '15X']
                        if e in df_out['Energy'].unique()]
        fss = sorted(df_out['FS'].unique())

        def weighted_pass(g):
            t = g['TotalPoints'].sum()
            f = g['FailPoints'].sum()
            return (t - f) / t * 100 if t > 0 else 0.0

        pivot = df_out.pivot_table(
            index='Energy', columns='FS',
            values='PassRate_pct', aggfunc='first'
        ).reindex(index=energy_order, columns=fss)

        all_fs   = df_out.groupby('Energy').apply(weighted_pass).reindex(energy_order)
        all_en   = df_out.groupby('FS').apply(weighted_pass).reindex(fss)
        all_data = weighted_pass(df_out)

        def fmt(v):
            return f"{v:.1f}%" if not np.isnan(v) else "—"

        fs_col_names = [f"{int(fs)} cm" if float(fs).is_integer() else f"{fs} cm"
                        for fs in fss]
        table_rows = []
        for e in energy_order:
            row = {'Energy': e}
            for fs, col in zip(fss, fs_col_names):
                row[col] = fmt(pivot.loc[e, fs])
            row['All Field Sizes'] = fmt(all_fs[e])
            table_rows.append(row)

        footer = {'Energy': 'All Energies'}
        for fs, col in zip(fss, fs_col_names):
            footer[col] = fmt(all_en[fs])
        footer['All Field Sizes'] = f"All Data: {fmt(all_data)}"
        table_rows.append(footer)

        df_table = pd.DataFrame(table_rows)

        try:
            df_table.to_csv(summary_csv, index=False)
            print(f"\nSummary CSV saved: {summary_csv}")
        except PermissionError:
            from datetime import datetime
            fallback = summary_csv.replace('.csv', f'_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
            df_table.to_csv(fallback, index=False)
            print(f"\nCSV locked (close it in Excel first next time).")
            print(f"Saved to fallback: {fallback}")
        print(df_table.to_string(index=False))
    else:
        print("\nNo results to save.")


if __name__ == '__main__':
    main()
