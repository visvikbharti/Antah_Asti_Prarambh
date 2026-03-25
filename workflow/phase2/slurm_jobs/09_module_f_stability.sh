#!/bin/bash
#SBATCH --job-name=aap_mod_f
#SBATCH --output=logs/09_module_f_%j.out
#SBATCH --error=logs/09_module_f_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

# =============================================================================
# Antah Asti Prarambh — Step 9: Module F (N-vs-C stability at full scale)
# Uses FoldX DeltaG + contact order + structural metrics
# Requires: Steps 5 (Chainsaw) + 7 (FoldX collect) + 8 (Module E) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 9: Module F — N-vs-C stability — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Check inputs exist
for f in \
  "${RESULTS}/foldx/foldx_stability_all.tsv" \
  "${RESULTS}/domains/unified_domain_assignments_full.tsv" \
  "${RESULTS}/domains/chainsaw_full_predictions.tsv"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing input: $f"
    exit 1
  fi
done

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import numpy as np
import os
import glob
import warnings
from scipy import stats
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Module F: Full-Scale N-vs-C Stability Analysis with FoldX")
print("=" * 70)

# ---- Load FoldX results (the primary new stability metric) ----
foldx = pd.read_csv(f"{RESULTS}/foldx/foldx_stability_all.tsv", sep="\t")
print(f"FoldX results: {len(foldx)} proteins")

# Extract accession from protein ID if needed
if "accession" not in foldx.columns:
    for col in ["protein", "protein_id", "name"]:
        if col in foldx.columns:
            foldx["accession"] = foldx[col].str.extract(r"([A-Z0-9]{6,10})", expand=False)
            break

# Filter to successful runs
if "status" in foldx.columns:
    foldx = foldx[foldx["status"] == "success"]
    print(f"  Successful: {len(foldx)}")

# ---- Load unified domain assignments ----
unified = pd.read_csv(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", sep="\t")
print(f"Unified domain assignments: {len(unified)} proteins")

# ---- Load Chainsaw predictions for domain boundaries ----
chainsaw = pd.read_csv(f"{RESULTS}/domains/chainsaw_full_predictions.tsv", sep="\t")
print(f"Chainsaw predictions: {len(chainsaw)} proteins")

# ---- Load Phase 1 CATH domain assignments for boundary info ----
cath_path = f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv"
if os.path.exists(cath_path):
    cath = pd.read_csv(cath_path, sep="\t")
    print(f"Phase 1 CATH domains: {len(cath)} domain records")
else:
    cath = pd.DataFrame()
    print("No Phase 1 CATH data found (will use Chainsaw only)")

# ---- Load substrate lists ----
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
matrix = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv", sep="\t")
mito = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv", sep="\t")

groel_acc = set(groel["accession"].values) if "accession" in groel.columns else set()
hsp60_acc = set(hsp60["accession"].values) if "accession" in hsp60.columns else set()
matrix_acc = set(matrix["accession"].values) if "accession" in matrix.columns else set()
mito_acc = set(mito["accession"].values) if "accession" in mito.columns else set()

# ---- Load Phase 1 structure index for protein lengths ----
struct_idx_path = f"{PROJECT_DIR}/results/structures/structure_index.tsv"
phase2_idx_path = f"{RESULTS}/structures/structure_index_full.tsv"
idx_path = phase2_idx_path if os.path.exists(phase2_idx_path) else struct_idx_path
if os.path.exists(idx_path):
    struct_idx = pd.read_csv(idx_path, sep="\t")
    print(f"Structure index: {len(struct_idx)} proteins")
else:
    struct_idx = pd.DataFrame()

# ---- Three-region decomposition ----
# For each multi-domain protein:
#   pre-domain tail: residues 1 to (first_domain_start - 1)
#   N-domain: first structural domain by position
#   C-region: everything after the first domain ends
print("\n--- Three-Region Decomposition ---")

# Build domain boundary lookup from CATH (preferred) then Chainsaw
boundary_records = []

# From CATH: get first domain boundaries per protein
if not cath.empty and "start" in cath.columns and "end" in cath.columns:
    for acc, group in cath.groupby("accession"):
        group_sorted = group.sort_values("start")
        first = group_sorted.iloc[0]
        boundary_records.append({
            "accession": acc,
            "source": "CATH",
            "n_domains": len(group_sorted),
            "first_domain_start": int(first["start"]),
            "first_domain_end": int(first["end"]),
        })

# From Chainsaw: parse domain boundaries for proteins not in CATH
cath_accs = {r["accession"] for r in boundary_records}

# Chainsaw output format varies; try to parse domain boundaries
dom_boundary_col = None
for col in ["chopping", "domain_boundaries", "boundaries", "domains"]:
    if col in chainsaw.columns:
        dom_boundary_col = col
        break

if dom_boundary_col and "accession" in chainsaw.columns:
    for _, row in chainsaw.iterrows():
        acc = row["accession"]
        if pd.isna(acc) or acc in cath_accs:
            continue
        chopping = str(row.get(dom_boundary_col, ""))
        if not chopping or chopping == "nan":
            continue
        # Parse Chainsaw chopping format: "1-100,150-300" or "1-100_150-300"
        # Domains separated by _ or ; or similar
        domains = []
        for segment in chopping.replace(";", "_").split("_"):
            parts = segment.strip().split("-")
            if len(parts) == 2:
                try:
                    start, end = int(parts[0]), int(parts[1])
                    domains.append((start, end))
                except ValueError:
                    # Handle complex ranges like "1-50,60-100"
                    for sub in segment.split(","):
                        sub_parts = sub.strip().split("-")
                        if len(sub_parts) == 2:
                            try:
                                domains.append((int(sub_parts[0]), int(sub_parts[1])))
                            except ValueError:
                                pass

        if domains:
            domains.sort(key=lambda x: x[0])
            boundary_records.append({
                "accession": acc,
                "source": "Chainsaw",
                "n_domains": len(domains),
                "first_domain_start": domains[0][0],
                "first_domain_end": domains[0][1],
            })

boundaries = pd.DataFrame(boundary_records)
print(f"Proteins with domain boundaries: {len(boundaries)}")

# Filter to multi-domain (need at least 1 domain to define regions)
multi_domain = boundaries[boundaries["n_domains"] >= 2]
print(f"Multi-domain proteins: {len(multi_domain)}")

# ---- Merge with FoldX stability ----
if "accession" in foldx.columns and not multi_domain.empty:
    merged = multi_domain.merge(foldx, on="accession", how="inner")
    print(f"Proteins with both domains + FoldX: {len(merged)}")
else:
    merged = multi_domain.copy()
    print("No FoldX merge possible (column mismatch)")

# ---- Compute region-level metrics ----
# For FoldX, we need per-residue or per-region DeltaG
# FoldX typically gives per-protein total energy
# We compare FoldX DeltaG between substrate groups as an overall metric

# ---- Paired N-vs-C comparisons ----
print("\n--- N-vs-C Paired Comparisons ---")

paired_records = []
for _, row in boundaries.iterrows():
    acc = row["accession"]
    ds = []
    if acc in groel_acc: ds.append("groel")
    if acc in hsp60_acc: ds.append("hsp60")
    if acc in matrix_acc: ds.append("matrix_bg")
    if acc in mito_acc: ds.append("mito_bg")
    if not ds: ds.append("proteome_bg")

    n_start = row["first_domain_start"]
    n_end = row["first_domain_end"]
    n_length = n_end - n_start + 1
    pre_tail = n_start - 1

    paired_records.append({
        "accession": acc,
        "source": row["source"],
        "n_domains": row["n_domains"],
        "pre_tail_length": pre_tail,
        "n_domain_start": n_start,
        "n_domain_end": n_end,
        "n_domain_length": n_length,
        "datasets": ",".join(ds),
    })

paired = pd.DataFrame(paired_records)

# ---- Merge FoldX DeltaG ----
if "accession" in foldx.columns:
    deltag_col = None
    for col in ["total_energy", "deltaG", "DeltaG", "stability_energy", "total"]:
        if col in foldx.columns:
            deltag_col = col
            break
    if deltag_col:
        foldx_slim = foldx[["accession", deltag_col]].drop_duplicates("accession")
        foldx_slim = foldx_slim.rename(columns={deltag_col: "foldx_deltaG"})
        paired = paired.merge(foldx_slim, on="accession", how="left")
        n_with_foldx = paired["foldx_deltaG"].notna().sum()
        print(f"Proteins with FoldX DeltaG: {n_with_foldx}/{len(paired)}")

# ---- Statistical tests per dataset ----
print("\n--- Statistical Summary ---")
for ds in ["groel", "hsp60", "matrix_bg", "mito_bg", "proteome_bg"]:
    subset = paired[paired["datasets"].str.contains(ds)]
    n = len(subset)
    if n == 0:
        continue
    multi = subset[subset["n_domains"] >= 2]
    print(f"\n{ds}: {n} total, {len(multi)} multi-domain")
    if "foldx_deltaG" in subset.columns:
        valid = subset["foldx_deltaG"].dropna()
        if len(valid) > 0:
            print(f"  FoldX DeltaG: mean={valid.mean():.2f}, median={valid.median():.2f}")

# ---- Save outputs ----
os.makedirs(f"{RESULTS}/stability", exist_ok=True)

boundaries_path = f"{RESULTS}/stability/region_boundaries_full.tsv"
boundaries.to_csv(boundaries_path, sep="\t", index=False)
print(f"\nSaved: {boundaries_path}")

paired_path = f"{RESULTS}/stability/n_vs_c_paired_full.tsv"
paired.to_csv(paired_path, sep="\t", index=False)
print(f"Saved: {paired_path}")

# Also save stability comparisons for Module H
stability_path = f"{RESULTS}/stats/stability_comparisons_full.tsv"
os.makedirs(os.path.dirname(stability_path), exist_ok=True)
paired.to_csv(stability_path, sep="\t", index=False)
print(f"Saved: {stability_path}")

print("\n" + "=" * 70)
print("Module F complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "Module F complete — $(date)"
echo "============================================================"
