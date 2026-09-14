"""
Extended Data Figure 2: relationship between local and global structural comparisons of tRNA
synthetases, classified by Class pairs. All structure pairs are used for comparison, and the
best (lowest) result (RMSD) is kept per protein pair. Structural pairs having a coverage lower
than 10% are discarded. A whole protein pair is discarded if all their associated structure
pairs are discarded for this reason.

Reconstructs a scatter plot originally produced ad hoc (never saved as a script) while tuning
figure_1_calculations.py's GLOBAL_COVERAGE_MIN=0.10 threshold for its structural RMSD matrix -
this figure shows the per-pair values behind that matrix's upper/lower triangles directly, split
by aaRS class combination.

Usage:
    python ExtendedDataFigure2.py
"""
import itertools
import json
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pandas as pd
import stylia

from default import AARS_CLASS_LABELS

# Format: print | Style: article - matches ExtendedDataFigure1.py (publication supplementary figure)
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "..", "output", "plots", "FigSupp", "ExtendedDataFigure2")
os.makedirs(plots_dir, exist_ok=True)

COMPARISONS_DIR = os.path.join(output_dir, "structural_comparisons_full")

# Canonical gene ordering/coloring, reused from figure_1_calculations.py so genes line up with
# every other figure in this project.
with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    figure_1_mappings = json.load(f)
uniprot_to_gene = figure_1_mappings["uniprot_to_gene"]

# Same aaRS Class I/II mapping as figure_1_plot.py, now centralized in src/default.py
# (AARS_CLASS_LABELS) - gatA/gatB have no aaRS class (transamidases, not ligases) and are
# excluded from this figure entirely.
classified_uniprots = sorted(uid for uid, gene in uniprot_to_gene.items() if gene in AARS_CLASS_LABELS)

# Same coverage cutoff and rationale as figure_1_calculations.py's structural RMSD matrix: a
# plain unfiltered min is misleading (low-coverage alignments can produce spuriously low RMSD),
# so the global (cmd.super) min is restricted to combinations clearing this coverage - a pair
# with no combination clearing it is dropped (NaN), not silently backfilled.
GLOBAL_COVERAGE_MIN = 0.10

CATEGORY_ORDER = ["I vs I", "II vs II", "I vs II"]


def pair_category(gene_1, gene_2):
    # Shortened from AARS_CLASS_LABELS' "Class I"/"Class II" to fit the legend at this figure's
    # panel size (user-requested).
    classes = sorted(AARS_CLASS_LABELS[gene].replace("Class ", "") for gene in (gene_1, gene_2))
    return f"{classes[0]} vs {classes[1]}"


rows = []
dropped_pairs = []
for uid1, uid2 in itertools.combinations(classified_uniprots, 2):
    ac_lo, ac_hi = sorted((uid1, uid2))
    path = os.path.join(COMPARISONS_DIR, f"{ac_lo}_{ac_hi}_global_local.csv")
    df = pd.read_csv(path)

    gene_1, gene_2 = uniprot_to_gene[uid1], uniprot_to_gene[uid2]
    filtered = df[df["global_coverage"] >= GLOBAL_COVERAGE_MIN]
    if not len(filtered):
        dropped_pairs.append((gene_1, gene_2))
        continue

    global_rmsd = filtered["global_rmsd"].min()
    local_rmsd = df["local_rmsd"].min()
    category = pair_category(gene_1, gene_2)
    rows.append((gene_1, gene_2, category, global_rmsd, local_rmsd))

print(f"Pairs dropped for no combination clearing {GLOBAL_COVERAGE_MIN:.0%} global coverage: {dropped_pairs}")

pair_data = pd.DataFrame(rows, columns=["gene_1", "gene_2", "category", "global_rmsd", "local_rmsd"])
print(f"Category counts: {pair_data['category'].value_counts().to_dict()}")

data_path = os.path.join(plots_dir, "ExtendedDataFigure2_data.csv")
pair_data.to_csv(data_path, index=False)
print(f"Saved {len(pair_data)} row(s) to {data_path}")

# ===========================================================================
# Plot
# ===========================================================================
pal = stylia.CategoricalPalette("npg")
category_colors = dict(zip(CATEGORY_ORDER, pal.get(len(CATEGORY_ORDER))))


def plot_global_vs_local_rmsd(ax, pair_data):
    for category in CATEGORY_ORDER:
        subset = pair_data[pair_data["category"] == category]
        ax.scatter(subset["global_rmsd"], subset["local_rmsd"], color=category_colors[category],
                   label=f"{category} (n={len(subset)})")
    stylia.label(ax, xlabel="Global RMSD,\n10%-coverage-filtered min (Å)", ylabel="Local RMSD, min (Å)",
                 title="Global vs. local structural RMSD\nper protein pair, by Class I/II category")
    # box_aspect(1) forces a physically square plot box (square data space, per stylia convention).
    ax.set_box_aspect(1)
    # Placed to the right of the axes rather than inside - at this panel's smaller size (matched to
    # the other square-scatter panels, see main()) no in-axes placement (auto "best" or a pinned
    # corner) was free of data: this panel's 3 categories don't fit any gap in the data without
    # touching a point (user-flagged after two prior attempts).
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)


def main():
    # width/height=0.3 matches the per-panel size used by the other square-scatter panels in this
    # supplementary series (ExtendedDataFigure4/5/6 split the same 0.6-fraction canvas across 2-4
    # panels) - 0.6/0.6 here was sized for a single panel and rendered ~4x larger in area than the
    # rest of the series (user-flagged).
    fig, axs = stylia.create_figure(1, 1, width=0.35, height=0.35)
    ax = axs.next()
    plot_global_vs_local_rmsd(ax, pair_data)

    png_path = os.path.join(plots_dir, "ExtendedDataFigure2.png")
    stylia.save_figure(png_path)
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
