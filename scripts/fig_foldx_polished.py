#!/usr/bin/env python3
"""
Generate polished FoldX DeltaG figure for publication.
Fig 7: FoldX Thermodynamic Stability Comparison

Three panels:
  A. Violin/box plot: FoldX total energy by dataset group
  B. GroEL substrates vs background (strip + box)
  C. FoldX DeltaG distribution with substrate overlay
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

PROJECT_DIR = os.path.expanduser("~/Downloads/Antah_Asti_Prarambh")
RESULTS = f"{PROJECT_DIR}/results/phase2"
FIG_DIR = f"{RESULTS}/figures"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Colorblind-friendly palette (Wong 2011)
COLORS = {
    "groel": "#0072B2",
    "hsp60": "#D55E00",
    "matrix_bg": "#009E73",
    "mito_bg": "#56B4E9",
    "proteome_bg": "#999999"
}

print("=" * 70)
print("Generating polished FoldX figure")
print("=" * 70)

# Load data
foldx = pd.read_csv(f"{RESULTS}/foldx/foldx_stability_all.tsv", sep="\t")
foldx = foldx[foldx["status"] == "success"].copy()
print(f"FoldX data: {len(foldx)} proteins")

paired = pd.read_csv(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", sep="\t")
print(f"N-vs-C paired: {len(paired)} proteins")

# Load substrate lists
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
matrix = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv", sep="\t")
mito = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv", sep="\t")

groel_acc = set(groel["current_accession"].dropna())
hsp60_acc = set(hsp60["uniprot_id"].dropna())
matrix_acc = set(matrix["uniprot_id"].dropna())
mito_acc = set(mito["uniprot_id"].dropna())

# Assign dataset labels to FoldX data
def assign_dataset(acc):
    if acc in groel_acc: return "GroEL"
    if acc in hsp60_acc: return "HSP60"
    if acc in matrix_acc: return "Matrix bg"
    if acc in mito_acc: return "Mito bg"
    return "Proteome bg"

foldx["dataset"] = foldx["accession"].apply(assign_dataset)

# Stats
groel_vals = foldx[foldx["dataset"] == "GroEL"]["total_energy"]
hsp60_vals = foldx[foldx["dataset"] == "HSP60"]["total_energy"]
bg_vals = foldx[foldx["dataset"] == "Proteome bg"]["total_energy"]
matrix_vals = foldx[foldx["dataset"] == "Matrix bg"]["total_energy"]

stat_gh, p_gh = stats.mannwhitneyu(groel_vals, bg_vals, alternative="two-sided")
stat_hh, p_hh = stats.mannwhitneyu(hsp60_vals, bg_vals, alternative="two-sided")

print(f"\nGroEL vs Proteome: U={stat_gh:.0f}, p={p_gh:.2e}, n={len(groel_vals)}")
print(f"HSP60 vs Proteome: U={stat_hh:.0f}, p={p_hh:.2e}, n={len(hsp60_vals)}")
print(f"GroEL median: {groel_vals.median():.1f}, HSP60 median: {hsp60_vals.median():.1f}, Proteome median: {bg_vals.median():.1f}")

# ---- Figure 7: FoldX DeltaG Comparison ----
fig = plt.figure(figsize=(14, 5))
gs = GridSpec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

# Panel A: Violin plot by group (clipped to reasonable range)
ax1 = fig.add_subplot(gs[0])
plot_data = foldx[foldx["dataset"].isin(["GroEL", "HSP60", "Matrix bg", "Mito bg"])].copy()
# Clip for visualization (outliers compress the plot)
plot_data["total_energy_clipped"] = plot_data["total_energy"].clip(-300, 1500)

order = ["GroEL", "HSP60", "Matrix bg", "Mito bg"]
palette = {"GroEL": "#0072B2", "HSP60": "#D55E00", "Matrix bg": "#009E73", "Mito bg": "#56B4E9"}

sns.violinplot(data=plot_data, x="dataset", y="total_energy_clipped", order=order,
               palette=palette, cut=0, inner="quartile", ax=ax1, linewidth=1)
ax1.set_xlabel("")
ax1.set_ylabel("FoldX Total Energy (kcal/mol)")
ax1.set_title("A. FoldX Stability by Group")

# Add sample sizes
for i, ds in enumerate(order):
    n = len(plot_data[plot_data["dataset"] == ds])
    med = plot_data[plot_data["dataset"] == ds]["total_energy"].median()
    ax1.text(i, ax1.get_ylim()[1] * 0.95, f"n={n}\nmed={med:.0f}",
             ha="center", va="top", fontsize=8, style="italic")

# Significance bracket for GroEL
y_max = 1400
ax1.plot([0, 0, 3, 3], [y_max, y_max + 30, y_max + 30, y_max], "k-", lw=1)
ax1.text(1.5, y_max + 40, f"p={p_gh:.1e}" if p_gh < 0.001 else f"p={p_gh:.3f}",
         ha="center", fontsize=9)
ax1.axhline(y=0, color="gray", linestyle="--", alpha=0.4, linewidth=0.8)

# Panel B: GroEL substrates vs E. coli background (focused)
ax2 = fig.add_subplot(gs[1])
ecoli_acc = set(pd.read_csv(f"{PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.tsv", sep="\t",
                            usecols=["Entry"])["Entry"].values)
ecoli_foldx = foldx[foldx["accession"].isin(ecoli_acc) & ~foldx["accession"].isin(groel_acc)]
groel_foldx = foldx[foldx["accession"].isin(groel_acc)]

if len(ecoli_foldx) > 0 and len(groel_foldx) > 0:
    ecoli_clip = ecoli_foldx["total_energy"].clip(-300, 600)
    groel_clip = groel_foldx["total_energy"].clip(-300, 600)

    data_b = pd.DataFrame({
        "FoldX Total Energy": pd.concat([groel_clip, ecoli_clip]),
        "Group": ["GroEL\nsubstrates"] * len(groel_clip) + ["E. coli\nbackground"] * len(ecoli_clip)
    })
    sns.boxplot(data=data_b, x="Group", y="FoldX Total Energy",
                palette=["#0072B2", "#999999"], width=0.5, ax=ax2,
                fliersize=2, linewidth=1)
    ax2.set_ylabel("FoldX Total Energy (kcal/mol)")
    ax2.set_xlabel("")
    ax2.set_title("B. GroEL vs E. coli Background")

    stat_ec, p_ec = stats.mannwhitneyu(groel_foldx["total_energy"], ecoli_foldx["total_energy"])
    ax2.text(0.5, 0.95, f"Mann-Whitney p={p_ec:.1e}\nn={len(groel_foldx)} vs {len(ecoli_foldx)}",
             transform=ax2.transAxes, ha="center", va="top", fontsize=9,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))
    ax2.axhline(y=0, color="gray", linestyle="--", alpha=0.4, linewidth=0.8)

# Panel C: Distribution overlay
ax3 = fig.add_subplot(gs[2])
bins = np.linspace(-300, 1000, 60)
ax3.hist(bg_vals.clip(-300, 1000), bins=bins, alpha=0.4, color="#999999",
         label=f"Proteome (n={len(bg_vals)})", density=True, edgecolor="none")
ax3.hist(groel_vals.clip(-300, 1000), bins=bins, alpha=0.7, color="#0072B2",
         label=f"GroEL (n={len(groel_vals)})", density=True, edgecolor="none")
ax3.hist(hsp60_vals.clip(-300, 1000), bins=bins, alpha=0.6, color="#D55E00",
         label=f"HSP60 (n={len(hsp60_vals)})", density=True, edgecolor="none")

ax3.axvline(groel_vals.median(), color="#0072B2", linestyle="--", linewidth=1.5, label=f"GroEL med={groel_vals.median():.0f}")
ax3.axvline(hsp60_vals.median(), color="#D55E00", linestyle="--", linewidth=1.5, label=f"HSP60 med={hsp60_vals.median():.0f}")
ax3.set_xlabel("FoldX Total Energy (kcal/mol)")
ax3.set_ylabel("Density")
ax3.set_title("C. FoldX Energy Distributions")
ax3.legend(fontsize=7, loc="upper right")

plt.suptitle("FoldX Thermodynamic Stability of Chaperonin Substrates", fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()

out_path = os.path.join(FIG_DIR, "fig7_foldx_stability")
fig.savefig(f"{out_path}.pdf", bbox_inches="tight", dpi=300)
fig.savefig(f"{out_path}.png", bbox_inches="tight", dpi=300)
plt.close()
print(f"\nSaved: {out_path}.pdf/png")

print("\nFoldX figure complete.")
