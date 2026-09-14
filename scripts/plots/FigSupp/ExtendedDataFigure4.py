"""
Extended Data Figure 4: overall shape of the per-pocket docking-score distributions across
the 3 screening libraries (HLL, REAL 9.56M, REAL 9.92B), stepping back from Figure 2c's
per-gene grid to show all 276 pockets pooled together. Reuses figure_2_calculations.py's
already-computed per-pocket percentiles (output/plots/figure_2/figure_2c_docking_percentiles.csv)
rather than recomputing anything.

(a) Median docking score, 3 distributions (HLL / REAL 9.56M / REAL 9.92B), 276 pockets each.
(b) Percentile-1 (1st-percentile) docking score, same 3 distributions.

Each distribution: a violin (colored per library, matching figure_2_plot.py's LIBRARIES, with a
black outline and median line).

Usage:
    python ExtendedDataFigure4.py
"""
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import matplotlib.colors as mcolors
import pandas as pd
import stylia
from stylia.config import get_fg_color

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataFigure4")
os.makedirs(plots_dir, exist_ok=True)

DOCKING_PERCENTILES_CSV = os.path.join(output_dir, "plots", "figure_2", "figure_2c_docking_percentiles.csv")

# (CSV library value, display label, NamedColors name) - same 3 libraries/colors as
# figure_2_plot.py's LIBRARIES, so colors stay consistent with Figure 2.
LIBRARY_DISPLAY = [
    ("HL", "HLL", "crimson"),
    ("REAL 10M", "REAL 9.56M", "turquoise"),
    ("REAL 10B", "REAL 9.92B", "amber"),
]

VIOLIN_WIDTH = 0.5

# Panel letter sits above each panel's own full rendered extent (axes + tick labels +
# title), same tightbbox-based approach as ExtendedDataFigure1's add_panel_label.
PANEL_LABEL_Y_PAD = 0.015


def add_panel_label(fig, ax, letter):
    renderer = fig.canvas.get_renderer()
    bbox = ax.get_tightbbox(renderer)
    x0 = fig.transFigure.inverted().transform((bbox.x0, 0))[0]
    top_y = fig.transFigure.inverted().transform((0, bbox.y1))[1]
    fig.text(x0, top_y + PANEL_LABEL_Y_PAD, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="bottom",
              transform=fig.transFigure)


def plot_metric_panel(ax, df, metric, ylabel):
    nc = stylia.NamedColors()
    positions = range(len(LIBRARY_DISPLAY))
    box_data = [df.loc[df["library"] == lib, metric].values for lib, _, _ in LIBRARY_DISPLAY]

    # ax.violinplot sets a default alpha=0.3 on each body Collection itself (separate from any
    # per-color alpha); that stored value silently overrides the alpha of whatever facecolor/
    # edgecolor is set afterwards (black at 0.3 opacity on white renders as RGB ~178, i.e. gray -
    # user-flagged). set_alpha(None) clears it so the RGBA alpha we assign below actually sticks.
    vp = ax.violinplot(box_data, positions=positions, widths=VIOLIN_WIDTH, showmedians=True,
                        showextrema=False)
    for i, body in enumerate(vp["bodies"]):
        body.set_alpha(None)
        body.set_facecolor(mcolors.to_rgba(nc.get(LIBRARY_DISPLAY[i][2]), alpha=0.95))
        body.set_edgecolor((0, 0, 0, 1))
        body.set_linewidth(stylia.LINEWIDTH)
        body.set_zorder(1)
    vp["cmedians"].set_color("black")
    vp["cmedians"].set_linewidth(stylia.LINEWIDTH)
    vp["cmedians"].set_zorder(2)

    ax.set_xticks(list(positions))
    ax.set_xticklabels([label for _, label, _ in LIBRARY_DISPLAY])
    stylia.label(ax, xlabel="", ylabel=ylabel)
    ax.set_box_aspect(1)
    # Shared range across both panels for direct comparison (user-requested) - clips 2/276 pockets
    # per library in panel a (median, up to ~3.6) out of the visible frame; panel b (Percentile-1)
    # has no points above -2. Accepted as a deliberate trade-off for a common scale, not a
    # data-loss bug.
    ax.set_ylim(-12, -2)


def main():
    df = pd.read_csv(DOCKING_PERCENTILES_CSV)

    data_path = os.path.join(plots_dir, "ExtendedDataFigure4_data.csv")
    df[["pocket", "gene", "library", "n", "median", "p1"]].to_csv(data_path, index=False)
    print(f"Saved {len(df)} row(s) to {data_path}")

    fig, axs = stylia.create_figure(1, 2, width=0.6, height=0.3)
    panel_axes = []

    ax_a = axs.next()
    plot_metric_panel(ax_a, df, "median", "Median docking score")
    panel_axes.append(ax_a)

    ax_b = axs.next()
    plot_metric_panel(ax_b, df, "p1", "Percentile-1 docking score")
    panel_axes.append(ax_b)

    # Needs a real renderer for add_panel_label's tightbbox measurements.
    fig.canvas.draw()
    for ax, letter in zip(panel_axes, ["a", "b"]):
        add_panel_label(fig, ax, letter)

    pdf_path = os.path.join(plots_dir, "ExtendedDataFigure4.pdf")
    png_path = os.path.join(plots_dir, "ExtendedDataFigure4.png")
    stylia.save_figure(pdf_path)
    stylia.save_figure(png_path)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
