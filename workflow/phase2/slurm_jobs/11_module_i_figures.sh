#!/bin/bash
#SBATCH --job-name=aap_mod_i
#SBATCH --output=logs/11_module_i_%j.out
#SBATCH --error=logs/11_module_i_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

# =============================================================================
# Antah Asti Prarambh — Step 11: Module I (Publication figures)
# Generates Phase 2 figures from full-scale analysis results
# Requires: Step 10 (Module H stats) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 11: Module I — Figures — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Use non-interactive matplotlib backend (no display on HPC)
export MPLBACKEND=Agg

python3 << 'PYTHON_SCRIPT'
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"
FIG_DIR = f"{RESULTS}/figures"
os.makedirs(FIG_DIR, exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 12, "axes.titlesize": 13,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

print("=" * 70)
print("Module I: Phase 2 Figure Generation")
print("=" * 70)

# ---- Load data ----
def safe_load(path, desc):
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        print(f"  Loaded {desc}: {len(df)} records")
        return df
    print(f"  [SKIP] {desc} not found: {path}")
    return None

print("\nLoading data...")
domains = safe_load(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", "domains")
dist = safe_load(f"{RESULTS}/domains/domain_distribution_full.tsv", "distribution")
paired = safe_load(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", "N-vs-C paired")
pvalues = safe_load(f"{RESULTS}/stats/corrected_pvalues_full.tsv", "p-values")

# Phase 1 reference data
groel = safe_load(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", "GroEL")
hsp60 = safe_load(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", "HSP60")

def get_accessions(df, preferred_cols=["accession", "current_accession", "uniprot_accession", "uniprot_id"]):
    """Extract accession set from a DataFrame, trying multiple column names."""
    if df is None:
        return set()
    for col in preferred_cols:
        if col in df.columns:
            return set(df[col].dropna().values)
    return set()

groel_acc = get_accessions(groel)
hsp60_acc = get_accessions(hsp60)

figures_generated = 0

# ---- Figure 1: Domain Distribution ----
if dist is not None and len(dist) > 0:
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel A: domain count distribution by dataset
        ax = axes[0]
        datasets = dist["dataset"].unique()
        x = np.arange(max(dist["n_domains"].max() + 1, 6))
        width = 0.15
        colors = sns.color_palette("Set2", len(datasets))

        for i, ds in enumerate(datasets):
            ds_data = dist[dist["dataset"] == ds]
            heights = [ds_data[ds_data["n_domains"] == nd]["percent"].sum() for nd in x]
            ax.bar(x + i * width, heights, width, label=ds, color=colors[i])

        ax.set_xlabel("Number of Domains")
        ax.set_ylabel("Percentage")
        ax.set_title("A. Domain Count Distribution")
        ax.legend(fontsize=8)
        ax.set_xticks(x + width * len(datasets) / 2)
        ax.set_xticklabels([str(int(v)) for v in x])

        # Panel B: CATH vs Chainsaw source
        ax = axes[1]
        if domains is not None and "source" in domains.columns:
            source_counts = domains["source"].value_counts()
            ax.pie(source_counts.values, labels=source_counts.index,
                   autopct="%1.1f%%", colors=["#4ECDC4", "#FF6B6B"])
            ax.set_title("B. Domain Assignment Source")

        plt.tight_layout()
        for fmt in ["pdf", "png"]:
            fig.savefig(f"{FIG_DIR}/fig1_domain_distribution_full.{fmt}")
        plt.close(fig)
        figures_generated += 1
        print(f"\nGenerated: fig1_domain_distribution_full")
    except Exception as e:
        print(f"Error generating Figure 1: {e}")

# ---- Figure 2: N-vs-C Stability with FoldX ----
if paired is not None and "foldx_deltaG" in paired.columns:
    try:
        valid = paired[paired["foldx_deltaG"].notna()]
        if len(valid) > 10:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # Panel A: FoldX DeltaG by dataset
            ax = axes[0]
            plot_data = []
            for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
                sub = valid[valid["accession"].isin(ds_accs)]
                bg = valid[~valid["accession"].isin(ds_accs)]
                for _, row in sub.iterrows():
                    plot_data.append({"Group": f"{ds_name}\nsubstrates", "FoldX DeltaG": row["foldx_deltaG"]})
                for _, row in bg.head(500).iterrows():
                    plot_data.append({"Group": "Background", "FoldX DeltaG": row["foldx_deltaG"]})

            if plot_data:
                plot_df = pd.DataFrame(plot_data)
                sns.violinplot(data=plot_df, x="Group", y="FoldX DeltaG", ax=ax,
                              palette="Set2", cut=0)
                ax.set_title("A. FoldX Stability by Group")
                ax.set_ylabel("FoldX DeltaG (kcal/mol)")

            # Panel B: N-domain length distribution
            ax = axes[1]
            multi = paired[paired["n_domains"] >= 2]
            if "n_domain_length" in multi.columns and len(multi) > 0:
                for ds_name, ds_accs, color in [("GroEL", groel_acc, "#4ECDC4"),
                                                  ("HSP60", hsp60_acc, "#FF6B6B")]:
                    sub = multi[multi["accession"].isin(ds_accs)]
                    if len(sub) > 0:
                        ax.hist(sub["n_domain_length"], bins=30, alpha=0.5,
                                label=ds_name, color=color, density=True)
                ax.set_xlabel("N-domain Length (residues)")
                ax.set_ylabel("Density")
                ax.set_title("B. First Domain Length Distribution")
                ax.legend()

            plt.tight_layout()
            for fmt in ["pdf", "png"]:
                fig.savefig(f"{FIG_DIR}/fig2_n_vs_c_stability_full.{fmt}")
            plt.close(fig)
            figures_generated += 1
            print(f"Generated: fig2_n_vs_c_stability_full")
    except Exception as e:
        print(f"Error generating Figure 2: {e}")
elif paired is not None:
    print("Skipping Figure 2: no FoldX DeltaG data in paired results")

# ---- Figure 3: FoldX DeltaG Distribution ----
if paired is not None and "foldx_deltaG" in paired.columns:
    try:
        valid = paired[paired["foldx_deltaG"].notna()]
        if len(valid) > 10:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(valid["foldx_deltaG"], bins=50, color="#2C3E50", alpha=0.7, edgecolor="white")
            ax.axvline(valid["foldx_deltaG"].median(), color="red", linestyle="--",
                      label=f"Median: {valid['foldx_deltaG'].median():.1f}")
            ax.set_xlabel("FoldX DeltaG (kcal/mol)")
            ax.set_ylabel("Count")
            ax.set_title("FoldX Stability Distribution (Full Proteome)")
            ax.legend()
            plt.tight_layout()
            for fmt in ["pdf", "png"]:
                fig.savefig(f"{FIG_DIR}/fig3_foldx_deltag_comparison.{fmt}")
            plt.close(fig)
            figures_generated += 1
            print(f"Generated: fig3_foldx_deltag_comparison")
    except Exception as e:
        print(f"Error generating Figure 3: {e}")

# ---- Figure 4: Statistical Summary ----
if pvalues is not None and len(pvalues) > 0:
    try:
        fig, ax = plt.subplots(figsize=(12, 6))
        pvalues_sorted = pvalues.sort_values("p_bh")
        colors = ["#E74C3C" if sig else "#BDC3C7" for sig in pvalues_sorted["significant"]]
        y_pos = range(len(pvalues_sorted))
        ax.barh(y_pos, -np.log10(pvalues_sorted["p_bh"].clip(lower=1e-50)), color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(pvalues_sorted["hypothesis"].values, fontsize=8)
        ax.axvline(-np.log10(0.05), color="black", linestyle="--", alpha=0.5, label="BH=0.05")
        ax.set_xlabel("-log10(p_BH)")
        ax.set_title("Phase 2 Statistical Tests")
        ax.legend()
        plt.tight_layout()
        for fmt in ["pdf", "png"]:
            fig.savefig(f"{FIG_DIR}/fig4_statistics_summary_full.{fmt}")
        plt.close(fig)
        figures_generated += 1
        print(f"Generated: fig4_statistics_summary_full")
    except Exception as e:
        print(f"Error generating Figure 4: {e}")

print(f"\n{'=' * 70}")
print(f"Module I complete. Generated {figures_generated} figures in {FIG_DIR}")
print(f"{'=' * 70}")
PYTHON_SCRIPT

echo ""
echo "Module I complete — $(date)"
echo "============================================================"
