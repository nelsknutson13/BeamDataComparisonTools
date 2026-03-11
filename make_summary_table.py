"""
Publication-Ready PDD Pass Rate Summary Table
----------------------------------------------
Reads the pivot-format pdd_comparison_summary.csv produced by
batch_pdd_compare.py and renders a formatted table PNG.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─────────────────────────────────────────────
#  CONFIG
# ─────────────────────────────────────────────
SUMMARY_CSV = (
    r"C:\Users\nknutson\Box\Knutson\Research\Projects underway"
    r"\Radformation research\RADMC Photon VALIDATION\TrueBeam"
    r"\ProcessedData\Results\pdd_comparison_summary.csv"
)

DD_CRITERIA  = 2    # % — used only for caption
DTA_CRITERIA = 2    # mm — used only for caption

OUT_PNG = os.path.join(os.path.dirname(SUMMARY_CSV), "pdd_passrate_table.png")
DPI     = 300

# Colours
COL_HEADER_BG  = '#2c5f8a'   # dark blue — column headers
COL_HEADER_FG  = 'white'
ROW_HEADER_BG  = '#2c5f8a'   # dark blue — row (energy) headers
ROW_HEADER_FG  = 'white'
MARGINAL_BG    = '#d0e4f5'   # light blue — all marginal cells
CELL_BG        = 'white'
CELL_BG_ALT    = '#f7f7f7'   # subtle alternating row shade
CELL_FG        = '#1a1a1a'   # near-black for all pass rate values
# ─────────────────────────────────────────────


def build_table():
    df = pd.read_csv(SUMMARY_CSV)

    # Last row = "All Energies" footer, last col = "All Field Sizes" marginal
    energy_rows = df[df['Energy'] != 'All Energies'].copy()
    footer_row  = df[df['Energy'] == 'All Energies'].iloc[0]

    fs_cols      = [c for c in df.columns if c not in ('Energy', 'All Field Sizes')]
    row_labels   = list(energy_rows['Energy']) + ['All\nEnergies']
    col_labels   = fs_cols + ['All Field\nSizes']

    # Build cell value arrays
    cell_vals = []
    for _, row in energy_rows.iterrows():
        cell_vals.append([str(row[c]) for c in fs_cols] + [str(row['All Field Sizes'])])

    # Footer: replace "All Data: XX%" corner with two-line version
    corner_raw = str(footer_row['All Field Sizes'])   # e.g. "All Data: 99.0%"
    corner_str = corner_raw.replace(': ', ':\n')
    cell_vals.append([str(footer_row[c]) for c in fs_cols] + [corner_str])

    n_rows = len(row_labels)
    n_cols = len(col_labels)

    # ── figure ────────────────────────────────────────────────────────
    fig_w = max(10, 1.5 * n_cols + 2)
    fig_h = max(4,  0.7 * n_rows + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis('off')

    cell_h = 0.12
    cell_w = 1.0 / (n_cols + 1)

    y_header = 1.0 - 0.05
    y_starts = [y_header - cell_h * (i + 1) for i in range(n_rows)]

    def draw_cell(ax, x, y, w, h, text, bg, fg, fontsize=10, bold=False):
        rect = mpatches.FancyBboxPatch(
            (x, y - h), w, h,
            boxstyle='square,pad=0',
            facecolor=bg, edgecolor='#cccccc', linewidth=0.5,
            transform=ax.transAxes, clip_on=False
        )
        ax.add_patch(rect)
        ax.text(
            x + w / 2, y - h / 2, text,
            ha='center', va='center',
            fontsize=fontsize, color=fg, weight='bold' if bold else 'normal',
            transform=ax.transAxes, clip_on=False, linespacing=1.3
        )

    # — column headers —
    draw_cell(ax, 0, y_header, cell_w, cell_h, 'Energy',
              COL_HEADER_BG, COL_HEADER_FG, bold=True)
    for j, cl in enumerate(col_labels):
        draw_cell(ax, cell_w * (j + 1), y_header, cell_w, cell_h,
                  cl, COL_HEADER_BG, COL_HEADER_FG, bold=True)

    # — energy rows —
    for i, (rl, row) in enumerate(zip(row_labels[:-1], cell_vals[:-1])):
        y   = y_starts[i]
        alt = CELL_BG_ALT if i % 2 else CELL_BG
        draw_cell(ax, 0, y, cell_w, cell_h, rl,
                  ROW_HEADER_BG, ROW_HEADER_FG, bold=True)
        for j, val_str in enumerate(row):
            is_marginal = (j == n_cols - 1)
            bg = MARGINAL_BG if is_marginal else alt
            draw_cell(ax, cell_w * (j + 1), y, cell_w, cell_h,
                      val_str, bg, CELL_FG, bold=is_marginal)

    # — footer row (All Energies) —
    y = y_starts[-1]
    draw_cell(ax, 0, y, cell_w, cell_h, row_labels[-1],
              MARGINAL_BG, '#1a3a5c', bold=True)
    for j, val_str in enumerate(cell_vals[-1]):
        draw_cell(ax, cell_w * (j + 1), y, cell_w, cell_h,
                  val_str, MARGINAL_BG, CELL_FG, bold=True)

    # ── title and caption ─────────────────────────────────────────────
    ax.text(0.5, y_header + 0.01,
            f'PDD Composite Pass Rate: {DD_CRITERIA}% / {DTA_CRITERIA} mm',
            ha='center', va='bottom', fontsize=13, weight='bold',
            transform=ax.transAxes)

    caption = (
        f'Composite criterion: dose difference ≤ {DD_CRITERIA}%  or  '
        f'DTA ≤ {DTA_CRITERIA} mm.  '
        '"All Field Sizes" and "All Energies" are weighted by point count.'
    )
    bottom_y = y_starts[-1] - cell_h - 0.01
    ax.text(0.5, bottom_y, caption,
            ha='center', va='top', fontsize=8, color='#555555',
            transform=ax.transAxes, style='italic')

    plt.savefig(OUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Table saved: {OUT_PNG}")
    os.startfile(OUT_PNG)


if __name__ == '__main__':
    build_table()
