#!/bin/bash
#SBATCH --job-name=aap_co_full
#SBATCH --output=logs/13_contact_order_full_%j.out
#SBATCH --error=logs/13_contact_order_full_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00

# =============================================================================
# Antah Asti Prarambh — Step 13: Full-scale Contact Order
# Compute relative contact order for ALL 25,007 proteins
# Uses 3-region decomposition: pre-tail, N-domain, C-region
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 13: Full-scale Contact Order — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

python3 -u << 'PYTHON_SCRIPT'
import os
import sys
import re
import glob
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Full-Scale Contact Order + N-vs-C Three-Region Decomposition")
print("=" * 70)

# ---- Load domain boundaries (CATH preferred, Chainsaw fallback) ----
cath = pd.read_csv(f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv", sep="\t")
chainsaw = pd.read_csv(f"{RESULTS}/domains/chainsaw_full_predictions.tsv", sep="\t")

# Extract accession from Chainsaw chain_id
if "accession" not in chainsaw.columns and "chain_id" in chainsaw.columns:
    chainsaw["accession"] = chainsaw["chain_id"].str.extract(r"AF-([A-Z0-9]+)-F", expand=False)

# Build domain boundary lookup
boundary_records = []

# From CATH
cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
for acc, group in cath.groupby(cath_acc_col):
    group_sorted = group.sort_values("domain_start")
    domains = [(int(r["domain_start"]), int(r["domain_end"])) for _, r in group_sorted.iterrows()]
    boundary_records.append({
        "accession": acc, "source": "CATH", "n_domains": len(domains),
        "domains": domains,
        "first_domain_start": domains[0][0], "first_domain_end": domains[0][1],
    })

cath_accs = {r["accession"] for r in boundary_records}

# From Chainsaw
for _, row in chainsaw.iterrows():
    acc = row.get("accession")
    if pd.isna(acc) or acc in cath_accs:
        continue
    chopping = str(row.get("chopping", ""))
    if not chopping or chopping == "nan":
        continue
    domains = []
    for segment in chopping.replace(";", "_").split("_"):
        parts = segment.strip().split("-")
        if len(parts) == 2:
            try:
                domains.append((int(parts[0]), int(parts[1])))
            except ValueError:
                for sub in segment.split(","):
                    sp = sub.strip().split("-")
                    if len(sp) == 2:
                        try: domains.append((int(sp[0]), int(sp[1])))
                        except ValueError: pass
    if domains:
        domains.sort(key=lambda x: x[0])
        n_dom = int(row["ndom"]) if "ndom" in chainsaw.columns and pd.notna(row.get("ndom")) else len(domains)
        boundary_records.append({
            "accession": acc, "source": "Chainsaw", "n_domains": n_dom,
            "domains": domains,
            "first_domain_start": domains[0][0], "first_domain_end": domains[0][1],
        })

boundaries = pd.DataFrame(boundary_records)
print(f"Domain boundaries: {len(boundaries)} proteins")

# ---- Load substrate lists ----
def get_acc_set(path, cols=["accession", "current_accession", "uniprot_accession", "uniprot_id"]):
    df = pd.read_csv(path, sep="\t")
    for c in cols:
        if c in df.columns:
            return set(df[c].dropna())
    return set()

groel_acc = get_acc_set(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv")
hsp60_acc = get_acc_set(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv")
matrix_acc = get_acc_set(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv")
mito_acc = get_acc_set(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv")
print(f"Substrates: GroEL={len(groel_acc)}, HSP60={len(hsp60_acc)}, Matrix={len(matrix_acc)}, Mito={len(mito_acc)}")

# ---- Find all structure files ----
struct_files = {}
for d in [f"{PROJECT_DIR}/data/raw/alphafold/full/ecoli",
          f"{PROJECT_DIR}/data/raw/alphafold/full/human",
          f"{PROJECT_DIR}/data/raw/alphafold/pilot"]:
    if os.path.exists(d):
        for f in glob.glob(f"{d}/*.cif"):
            m = re.search(r'AF-([A-Z0-9]+)-F', os.path.basename(f))
            if m:
                struct_files[m.group(1)] = f
print(f"Structure files found: {len(struct_files)}")

# ---- Contact order calculation ----
def parse_ca_coords(cif_path):
    """Extract CA coordinates from CIF file."""
    coords = {}
    try:
        with open(cif_path) as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line:
                    parts = line.split()
                    if len(parts) >= 15:
                        try:
                            # CIF format: label_seq_id is residue number
                            # Try to find x, y, z coordinates
                            resnum = int(parts[8]) if parts[8].isdigit() else int(parts[6])
                            x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                            if resnum not in coords:  # take first alt conf
                                coords[resnum] = (x, y, z)
                        except (ValueError, IndexError):
                            pass
    except Exception:
        pass
    return coords

def compute_rco(coords, start=None, end=None, cutoff=8.0, min_sep=6):
    """Compute relative contact order for a region."""
    if start is None:
        residues = sorted(coords.keys())
    else:
        residues = sorted([r for r in coords.keys() if start <= r <= end])

    if len(residues) < 10:
        return np.nan, 0

    L = residues[-1] - residues[0] + 1
    total_sep = 0
    n_contacts = 0

    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:
            if abs(r2 - r1) < min_sep:
                continue
            c1, c2 = coords[r1], coords[r2]
            dist = np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)
            if dist <= cutoff:
                total_sep += abs(r2 - r1)
                n_contacts += 1

    if n_contacts == 0 or L == 0:
        return np.nan, n_contacts

    rco = total_sep / (n_contacts * L)
    return rco, n_contacts

# ---- Process all proteins ----
contact_order_records = []
paired_records = []

total = len(boundaries)
for idx, row in boundaries.iterrows():
    if (idx + 1) % 1000 == 0:
        print(f"  Progress: {idx+1}/{total} ({100*(idx+1)/total:.1f}%)")

    acc = row["accession"]
    if acc not in struct_files:
        continue

    coords = parse_ca_coords(struct_files[acc])
    if len(coords) < 20:
        continue

    # Whole-protein contact order
    rco_whole, nc_whole = compute_rco(coords)

    # Three-region decomposition
    n_start = row["first_domain_start"]
    n_end = row["first_domain_end"]

    # N-domain contact order
    rco_n, nc_n = compute_rco(coords, start=n_start, end=n_end)

    # C-region contact order (everything after first domain)
    c_start = n_end + 1
    c_end = max(coords.keys()) if coords else n_end
    rco_c, nc_c = compute_rco(coords, start=c_start, end=c_end)

    # Dataset assignment
    ds = []
    if acc in groel_acc: ds.append("groel")
    if acc in hsp60_acc: ds.append("hsp60")
    if acc in matrix_acc: ds.append("matrix_bg")
    if acc in mito_acc: ds.append("mito_bg")
    if not ds: ds.append("proteome_bg")

    contact_order_records.append({
        "accession": acc,
        "protein_length": max(coords.keys()) if coords else 0,
        "rco_whole": rco_whole,
        "n_contacts_whole": nc_whole,
        "rco_n_domain": rco_n,
        "n_contacts_n": nc_n,
        "rco_c_region": rco_c,
        "n_contacts_c": nc_c,
        "n_domains": row["n_domains"],
        "source": row["source"],
        "datasets": ",".join(ds),
    })

    # Paired record for multi-domain proteins
    if row["n_domains"] >= 2 and not np.isnan(rco_n) and not np.isnan(rco_c):
        paired_records.append({
            "accession": acc,
            "source": row["source"],
            "n_domains": row["n_domains"],
            "datasets": ",".join(ds),
            "rco_n_domain": rco_n,
            "rco_c_region": rco_c,
            "rco_diff": rco_n - rco_c,
            "n_domain_length": n_end - n_start + 1,
            "c_region_length": c_end - c_start + 1 if c_start <= c_end else 0,
            "pre_tail_length": n_start - 1,
        })

# ---- Save ----
os.makedirs(f"{RESULTS}/stability", exist_ok=True)

co_df = pd.DataFrame(contact_order_records)
co_path = f"{RESULTS}/stability/contact_order_full.tsv"
co_df.to_csv(co_path, sep="\t", index=False)
print(f"\nSaved contact order: {co_path} ({len(co_df)} proteins)")

paired_df = pd.DataFrame(paired_records)
paired_path = f"{RESULTS}/stability/n_vs_c_rco_paired_full.tsv"
paired_df.to_csv(paired_path, sep="\t", index=False)
print(f"Saved N-vs-C paired: {paired_path} ({len(paired_df)} multi-domain proteins)")

# Quick stats
print(f"\n--- Summary ---")
print(f"Total proteins with contact order: {len(co_df)}")
print(f"Multi-domain with N-vs-C pairs: {len(paired_df)}")
for ds in ["groel", "hsp60", "matrix_bg", "mito_bg", "proteome_bg"]:
    sub = paired_df[paired_df["datasets"].str.contains(ds, na=False)]
    if len(sub) > 0:
        print(f"  {ds}: {len(sub)} pairs, median RCO diff = {sub['rco_diff'].median():.4f}")

print("\n" + "=" * 70)
print("Full-scale contact order complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "Contact order complete — $(date)"
echo "============================================================"
