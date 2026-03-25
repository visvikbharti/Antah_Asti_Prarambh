#!/bin/bash
#SBATCH --job-name=aap_mod_e
#SBATCH --output=logs/08_module_e_%j.out
#SBATCH --error=logs/08_module_e_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

# =============================================================================
# Antah Asti Prarambh — Step 8: Module E (Domain architecture at full scale)
# v2 fix: corrected column names for all input files
# Integrates Phase 1 CATH + Phase 2 Chainsaw, computes domain distributions
# Requires: Steps 4 (Foldseek clusters) + 5 (Chainsaw) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 8: Module E — Domain architecture — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Check inputs exist
for f in \
  "${RESULTS}/domains/chainsaw_full_predictions.tsv" \
  "${RESULTS}/foldseek/analysis/foldseek_clusters_full.tsv"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing input: $f"
    exit 1
  fi
done

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import numpy as np
import os
from collections import Counter

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Module E: Full-Scale Domain Architecture Analysis (v2)")
print("=" * 70)

# ---- Load Phase 1 CATH assignments ----
cath_path = f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv"
cath_summary_path = f"{PROJECT_DIR}/results/domains/cath_protein_summary.tsv"

if os.path.exists(cath_path):
    cath = pd.read_csv(cath_path, sep="\t")
    print(f"Phase 1 CATH domains: {len(cath)} domain records")
    # Column is 'uniprot_accession' in Phase 1 CATH files
    cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
    cath_proteins = set(cath[cath_acc_col].unique())
    print(f"  Proteins with CATH: {len(cath_proteins)} (using column '{cath_acc_col}')")
else:
    print(f"WARNING: CATH file not found at {cath_path}")
    cath = pd.DataFrame()
    cath_proteins = set()
    cath_acc_col = "accession"

if os.path.exists(cath_summary_path):
    cath_summary = pd.read_csv(cath_summary_path, sep="\t")
    print(f"Phase 1 CATH summary: {len(cath_summary)} proteins")
    # Identify the accession column in summary too
    cath_summ_acc_col = "uniprot_accession" if "uniprot_accession" in cath_summary.columns else "accession"
    print(f"  Summary accession column: '{cath_summ_acc_col}'")
else:
    cath_summary = pd.DataFrame()
    cath_summ_acc_col = "accession"

# ---- Load Phase 2 Chainsaw predictions ----
chainsaw = pd.read_csv(f"{RESULTS}/domains/chainsaw_full_predictions.tsv", sep="\t")
print(f"\nPhase 2 Chainsaw predictions: {len(chainsaw)} proteins")

# Identify accession column
acc_col = None
for col in ["accession", "chain_id", "name", "protein"]:
    if col in chainsaw.columns:
        acc_col = col
        break
if acc_col is None:
    acc_col = chainsaw.columns[0]
print(f"  Using column '{acc_col}' for protein IDs")

# Extract accession from AlphaFold filename pattern if needed
if chainsaw[acc_col].str.contains("AF-").any():
    chainsaw["accession"] = chainsaw[acc_col].str.extract(r"AF-([A-Z0-9]+)-", expand=False)
else:
    chainsaw["accession"] = chainsaw[acc_col]

chainsaw_proteins = set(chainsaw["accession"].dropna().unique())
print(f"  Chainsaw proteins: {len(chainsaw_proteins)}")

# Quick stats on Chainsaw predictions
ndom_col = "ndom" if "ndom" in chainsaw.columns else "n_domains"
if ndom_col in chainsaw.columns:
    ndom_counts = chainsaw[ndom_col].value_counts().sort_index()
    print(f"  Domain count distribution:")
    for nd, cnt in ndom_counts.items():
        print(f"    ndom={nd}: {cnt} proteins")

# ---- Build unified domain assignments (CATH preferred, Chainsaw fills gaps) ----
print("\n--- Building unified domain assignments ---")

chainsaw_only = chainsaw_proteins - cath_proteins
both = chainsaw_proteins & cath_proteins
cath_only = cath_proteins - chainsaw_proteins
print(f"  CATH-only: {len(cath_only)}")
print(f"  Chainsaw-only (new): {len(chainsaw_only)}")
print(f"  Both (CATH preferred): {len(both)}")

# Build unified summary
unified_records = []

# Add CATH proteins (using correct column name)
if not cath_summary.empty:
    for _, row in cath_summary.iterrows():
        acc = row.get(cath_summ_acc_col, "")
        unified_records.append({
            "accession": acc,
            "source": "CATH",
            "n_domains": row.get("n_domains", 0),
            "top_superfamily": row.get("top_superfamily", row.get("domain_architecture", "")),
            "cath_class": row.get("cath_class", ""),
        })

# Add Chainsaw-only proteins
n_domain_col = None
for col in ["n_domains", "ndom", "num_domains", "domain_count"]:
    if col in chainsaw.columns:
        n_domain_col = col
        break

for _, row in chainsaw[chainsaw["accession"].isin(chainsaw_only)].iterrows():
    n_dom = int(row[n_domain_col]) if n_domain_col else 1
    unified_records.append({
        "accession": row["accession"],
        "source": "Chainsaw",
        "n_domains": n_dom,
        "top_superfamily": "",
        "cath_class": "",
    })

unified = pd.DataFrame(unified_records)
print(f"\nUnified domain assignments: {len(unified)} proteins")
print(f"  By source: {unified['source'].value_counts().to_dict()}")

# ---- Load substrate lists (using correct column names) ----
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
matrix = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv", sep="\t")
mito = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv", sep="\t")

# GroEL uses 'current_accession'
groel_acc_col = "current_accession" if "current_accession" in groel.columns else "accession"
groel_acc = set(groel[groel_acc_col].dropna().values)

# HSP60, matrix, mito use 'uniprot_id'
hsp60_acc_col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
hsp60_acc = set(hsp60[hsp60_acc_col].dropna().values)

matrix_acc_col = "uniprot_id" if "uniprot_id" in matrix.columns else "accession"
matrix_acc = set(matrix[matrix_acc_col].dropna().values)

mito_acc_col = "uniprot_id" if "uniprot_id" in mito.columns else "accession"
mito_acc = set(mito[mito_acc_col].dropna().values)

print(f"\nSubstrate counts: GroEL={len(groel_acc)} (col: {groel_acc_col}), HSP60={len(hsp60_acc)} (col: {hsp60_acc_col})")
print(f"Background counts: Matrix={len(matrix_acc)} (col: {matrix_acc_col}), Mito={len(mito_acc)} (col: {mito_acc_col})")

# ---- Assign dataset membership ----
def get_datasets(acc):
    ds = []
    if acc in groel_acc: ds.append("groel")
    if acc in hsp60_acc: ds.append("hsp60")
    if acc in matrix_acc: ds.append("matrix_bg")
    if acc in mito_acc: ds.append("mito_bg")
    if not ds: ds.append("proteome_bg")
    return ",".join(ds)

unified["datasets"] = unified["accession"].apply(get_datasets)

# ---- Domain distribution per dataset ----
print("\n--- Domain Distribution by Dataset ---")
datasets_of_interest = ["groel", "hsp60", "matrix_bg", "mito_bg", "proteome_bg"]
dist_records = []

for ds in datasets_of_interest:
    subset = unified[unified["datasets"].str.contains(ds)]
    if len(subset) == 0:
        print(f"\n{ds}: NO PROTEINS FOUND")
        continue
    counts = subset["n_domains"].value_counts().sort_index()
    total = len(subset)
    print(f"\n{ds} (n={total}):")
    for nd, cnt in counts.items():
        pct = 100 * cnt / total
        print(f"  {nd} domains: {cnt} ({pct:.1f}%)")
        dist_records.append({
            "dataset": ds,
            "n_domains": nd,
            "count": cnt,
            "percent": round(pct, 1),
        })

# ---- Load Foldseek clusters ----
clusters_path = f"{RESULTS}/foldseek/analysis/foldseek_clusters_full.tsv"
clusters = pd.read_csv(clusters_path, sep="\t")
print(f"\nFoldseek clusters: {len(clusters)} proteins")
if "cluster_id" in clusters.columns:
    n_clusters = clusters["cluster_id"].nunique()
    print(f"  Unique clusters: {n_clusters}")

# ---- Save outputs ----
os.makedirs(f"{RESULTS}/domains", exist_ok=True)

unified_path = f"{RESULTS}/domains/unified_domain_assignments_full.tsv"
unified.to_csv(unified_path, sep="\t", index=False)
print(f"\nSaved: {unified_path}")

dist_df = pd.DataFrame(dist_records)
dist_path = f"{RESULTS}/domains/domain_distribution_full.tsv"
dist_df.to_csv(dist_path, sep="\t", index=False)
print(f"Saved: {dist_path}")

# Also save the raw Chainsaw data with accession column for reference
chainsaw_ref_path = f"{RESULTS}/domains/chainsaw_full_predictions_annotated.tsv"
chainsaw.to_csv(chainsaw_ref_path, sep="\t", index=False)
print(f"Saved: {chainsaw_ref_path}")

print("\n" + "=" * 70)
print("Module E complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "Module E complete — $(date)"
echo "============================================================"
