#!/bin/bash
#SBATCH --job-name=aap_mod_h
#SBATCH --output=logs/10_module_h_%j.out
#SBATCH --error=logs/10_module_h_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

# =============================================================================
# Antah Asti Prarambh — Step 10: Module H (Comparative statistics)
# Full-scale statistical analysis with hierarchical BH correction
# Requires: Steps 8 (Module E) + 9 (Module F) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 10: Module H — Statistics — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Check inputs
for f in \
  "${RESULTS}/domains/domain_distribution_full.tsv" \
  "${RESULTS}/stability/n_vs_c_paired_full.tsv" \
  "${RESULTS}/domains/unified_domain_assignments_full.tsv"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing input: $f"
    exit 1
  fi
done

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import numpy as np
import os
import warnings
from scipy import stats
from datetime import datetime

warnings.filterwarnings("ignore")

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("WARNING: statsmodels not available, using manual BH correction")

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Module H: Phase 2 Full-Scale Comparative Statistics")
print(f"Run date: {datetime.now()}")
print("=" * 70)

# ---- Load Phase 2 data ----
domains = pd.read_csv(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", sep="\t")
dist = pd.read_csv(f"{RESULTS}/domains/domain_distribution_full.tsv", sep="\t")
paired = pd.read_csv(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", sep="\t")

print(f"\nDomain data: {len(domains)} proteins")
print(f"Distribution data: {len(dist)} records")
print(f"N-vs-C paired: {len(paired)} proteins")

# Try loading Foldseek clusters
clusters_path = f"{RESULTS}/foldseek/analysis/foldseek_clusters_full.tsv"
if os.path.exists(clusters_path):
    clusters = pd.read_csv(clusters_path, sep="\t")
    print(f"Cluster data: {len(clusters)} proteins")
else:
    clusters = pd.DataFrame()
    print("No cluster data available")

# ---- Load substrate lists ----
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
homologs = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv", sep="\t")

def get_accessions(df, preferred_cols=["accession", "current_accession", "uniprot_accession", "uniprot_id"]):
    """Extract accession set from a DataFrame, trying multiple column names."""
    for col in preferred_cols:
        if col in df.columns:
            return set(df[col].dropna().values)
    return set()

groel_acc = get_accessions(groel)
hsp60_acc = get_accessions(hsp60)

print(f"\nGroEL substrates: {len(groel_acc)}")
print(f"HSP60 substrates: {len(hsp60_acc)}")
print(f"Homolog pairs: {len(homologs)}")

# ---- Helper: BH correction ----
def bh_correct(pvalues):
    """Benjamini-Hochberg correction."""
    if HAS_STATSMODELS:
        reject, corrected, _, _ = multipletests(pvalues, method="fdr_bh")
        return corrected
    # Manual BH
    n = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    corrected = np.zeros(n)
    for rank, idx in enumerate(sorted_idx, 1):
        corrected[idx] = pvalues[idx] * n / rank
    # Enforce monotonicity
    corrected[sorted_idx[-1]] = min(corrected[sorted_idx[-1]], 1.0)
    for i in range(n - 2, -1, -1):
        corrected[sorted_idx[i]] = min(corrected[sorted_idx[i]], corrected[sorted_idx[i + 1]])
    return corrected

# ---- Helper: Cohen's d ----
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return 0.0
    pooled_std = np.sqrt(((nx - 1) * np.std(x, ddof=1)**2 + (ny - 1) * np.std(y, ddof=1)**2) / (nx + ny - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(x) - np.mean(y)) / pooled_std

# ===========================================================================
# FAMILY 1: Domain Architecture (Fisher exact tests)
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 1: Domain Architecture")
print("=" * 70)

all_tests = []

# Test: multi-domain enrichment in substrates
for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
    sub = domains[domains["accession"].isin(ds_accs)]
    bg = domains[~domains["accession"].isin(ds_accs)]

    if len(sub) == 0 or len(bg) == 0:
        continue

    # Multi-domain (>= 2) vs single-domain
    sub_multi = (sub["n_domains"] >= 2).sum()
    sub_single = (sub["n_domains"] < 2).sum()
    bg_multi = (bg["n_domains"] >= 2).sum()
    bg_single = (bg["n_domains"] < 2).sum()

    table = [[sub_multi, sub_single], [bg_multi, bg_single]]
    oddsratio, pval = stats.fisher_exact(table)
    d = cohens_d(sub["n_domains"].values, bg["n_domains"].values)

    print(f"\n{ds_name} multi-domain enrichment:")
    print(f"  Substrate: {sub_multi}/{len(sub)} ({100*sub_multi/len(sub):.1f}%) multi-domain")
    print(f"  Background: {bg_multi}/{len(bg)} ({100*bg_multi/len(bg):.1f}%) multi-domain")
    print(f"  OR={oddsratio:.3f}, p={pval:.2e}, Cohen's d={d:.3f}")

    all_tests.append({
        "family": "domain_architecture",
        "hypothesis": f"H1.1_{ds_name}_multi_domain",
        "test": "Fisher exact",
        "statistic": oddsratio,
        "p_value": pval,
        "effect_size": d,
        "n_substrate": len(sub),
        "n_background": len(bg),
    })

# ===========================================================================
# FAMILY 2: N-vs-C Stability Asymmetry
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 2: N-vs-C Stability Asymmetry")
print("=" * 70)

# Test FoldX DeltaG differences between substrates and backgrounds
if "foldx_deltaG" in paired.columns:
    for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
        sub = paired[paired["accession"].isin(ds_accs) & paired["foldx_deltaG"].notna()]
        bg = paired[~paired["accession"].isin(ds_accs) & paired["foldx_deltaG"].notna()]

        if len(sub) < 5 or len(bg) < 5:
            continue

        stat, pval = stats.mannwhitneyu(sub["foldx_deltaG"], bg["foldx_deltaG"], alternative="two-sided")
        d = cohens_d(sub["foldx_deltaG"].values, bg["foldx_deltaG"].values)

        print(f"\n{ds_name} FoldX DeltaG (substrate vs background):")
        print(f"  Substrate: n={len(sub)}, mean={sub['foldx_deltaG'].mean():.2f}")
        print(f"  Background: n={len(bg)}, mean={bg['foldx_deltaG'].mean():.2f}")
        print(f"  U={stat:.0f}, p={pval:.2e}, Cohen's d={d:.3f}")

        all_tests.append({
            "family": "stability_asymmetry",
            "hypothesis": f"H2.1_{ds_name}_foldx_deltaG",
            "test": "Mann-Whitney U",
            "statistic": stat,
            "p_value": pval,
            "effect_size": d,
            "n_substrate": len(sub),
            "n_background": len(bg),
        })

# N-domain length asymmetry
multi_paired = paired[paired["n_domains"] >= 2]
for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
    sub = multi_paired[multi_paired["accession"].isin(ds_accs)]
    bg = multi_paired[~multi_paired["accession"].isin(ds_accs)]

    if len(sub) < 5:
        continue

    # Pre-tail length comparison
    if "pre_tail_length" in sub.columns:
        stat, pval = stats.mannwhitneyu(
            sub["pre_tail_length"].dropna(), bg["pre_tail_length"].dropna(),
            alternative="two-sided"
        )
        d = cohens_d(sub["pre_tail_length"].dropna().values, bg["pre_tail_length"].dropna().values)
        print(f"\n{ds_name} pre-tail length:")
        print(f"  Substrate: mean={sub['pre_tail_length'].mean():.1f}")
        print(f"  Background: mean={bg['pre_tail_length'].mean():.1f}")
        print(f"  p={pval:.2e}, d={d:.3f}")

        all_tests.append({
            "family": "stability_asymmetry",
            "hypothesis": f"H2.2_{ds_name}_pre_tail",
            "test": "Mann-Whitney U",
            "statistic": stat,
            "p_value": pval,
            "effect_size": d,
            "n_substrate": len(sub),
            "n_background": len(bg),
        })

# ===========================================================================
# FAMILY 3: MTS Targeting
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 3: MTS Targeting Statistics")
print("=" * 70)

# Load Phase 1 targeting data (if available)
targ_path = f"{PROJECT_DIR}/results/mts/combined_targeting.tsv"
if os.path.exists(targ_path):
    targeting = pd.read_csv(targ_path, sep="\t")
    print(f"Targeting data: {len(targeting)} proteins")

    # Matrix enrichment in HSP60 substrates
    # Targeting file uses uniprot_accession
    targ_acc_col = None
    for col in ["accession", "uniprot_accession", "uniprot_id"]:
        if col in targeting.columns:
            targ_acc_col = col
            break

    if targ_acc_col:
        hsp60_targ = targeting[targeting[targ_acc_col].isin(hsp60_acc)]
        bg_targ = targeting[~targeting[targ_acc_col].isin(hsp60_acc)]

        matrix_col = None
        for col in ["mitocarta_is_matrix", "matrix_localization", "is_matrix", "localization", "sub_localization"]:
            if col in targeting.columns:
                matrix_col = col
                break

        if matrix_col:
            hsp60_matrix = hsp60_targ[hsp60_targ[matrix_col].astype(str).str.contains("True|matrix|Matrix|1", case=False, na=False)]
            bg_matrix = bg_targ[bg_targ[matrix_col].astype(str).str.contains("True|matrix|Matrix|1", case=False, na=False)]

            table = [[len(hsp60_matrix), len(hsp60_targ) - len(hsp60_matrix)],
                     [len(bg_matrix), len(bg_targ) - len(bg_matrix)]]
            oddsratio, pval = stats.fisher_exact(table)

            print(f"\nHSP60 matrix enrichment:")
            print(f"  HSP60 matrix: {len(hsp60_matrix)}/{len(hsp60_targ)}")
            print(f"  Background matrix: {len(bg_matrix)}/{len(bg_targ)}")
            print(f"  OR={oddsratio:.3f}, p={pval:.2e}")

            all_tests.append({
                "family": "mts_targeting",
                "hypothesis": "H3.1_HSP60_matrix_enrichment",
                "test": "Fisher exact",
                "statistic": oddsratio,
                "p_value": pval,
                "effect_size": oddsratio,
                "n_substrate": len(hsp60_targ),
                "n_background": len(bg_targ),
            })
else:
    print("No targeting data available (Phase 1 file not found)")

# ===========================================================================
# HIERARCHICAL BH CORRECTION
# ===========================================================================
print("\n" + "=" * 70)
print("HIERARCHICAL CORRECTION")
print("=" * 70)

if all_tests:
    results_df = pd.DataFrame(all_tests)

    # Apply BH correction within each family
    corrected_pvals = []
    for family in results_df["family"].unique():
        mask = results_df["family"] == family
        family_pvals = results_df.loc[mask, "p_value"].values
        if len(family_pvals) > 0:
            corrected = bh_correct(family_pvals)
            for idx, m in zip(results_df.loc[mask].index, range(len(corrected))):
                corrected_pvals.append({"index": idx, "p_bh": corrected[m]})

    corr_df = pd.DataFrame(corrected_pvals).set_index("index")
    results_df["p_bh"] = corr_df["p_bh"]
    results_df["significant"] = results_df["p_bh"] < 0.05

    print(f"\nTotal tests: {len(results_df)}")
    print(f"Significant (BH < 0.05): {results_df['significant'].sum()}")
    print(f"\nBy family:")
    for family in results_df["family"].unique():
        fam = results_df[results_df["family"] == family]
        sig = fam["significant"].sum()
        print(f"  {family}: {sig}/{len(fam)} significant")

    # Save results
    os.makedirs(f"{RESULTS}/stats", exist_ok=True)

    pval_path = f"{RESULTS}/stats/corrected_pvalues_full.tsv"
    results_df.to_csv(pval_path, sep="\t", index=False)
    print(f"\nSaved: {pval_path}")

    # Write summary report
    report_path = f"{RESULTS}/stats/statistics_summary_full.txt"
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("PHASE 2 FULL-SCALE STATISTICS SUMMARY\n")
        f.write(f"Generated: {datetime.now()}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Total tests: {len(results_df)}\n")
        f.write(f"Significant (BH < 0.05): {results_df['significant'].sum()}\n\n")

        for _, row in results_df.iterrows():
            sig_mark = " ***" if row["significant"] else ""
            f.write(f"{row['hypothesis']}: p={row['p_value']:.2e}, "
                    f"p_BH={row['p_bh']:.2e}, d={row['effect_size']:.3f}{sig_mark}\n")

        f.write("\n\nKey reminders:\n")
        f.write("  - pLDDT is NOT stability; FoldX DeltaG is the primary stability metric\n")
        f.write("  - N-vs-C asymmetry was universal in Phase 1 (verify with Phase 2)\n")
        f.write("  - Hierarchical BH correction across 3 goal families\n")

    print(f"Saved: {report_path}")
else:
    print("No statistical tests could be run (insufficient data)")

print("\n" + "=" * 70)
print("Module H complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "Module H complete — $(date)"
echo "============================================================"
