#!/bin/bash
#SBATCH --job-name=aap_cath_full
#SBATCH --output=logs/14_cath_interpro_%j.out
#SBATCH --error=logs/14_cath_interpro_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00

# =============================================================================
# Antah Asti Prarambh — Step 14: Full-scale CATH superfamily classification
# Query InterPro/Gene3D API for all 25,007 proteins to get CATH superfamily
# Rate-limited to 1 request/second; checkpointed for resume
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 14: Full-scale CATH superfamily — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

python3 -u << 'PYTHON_SCRIPT'
import os
import sys
import json
import time
import requests
import pandas as pd
from pathlib import Path

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"

os.makedirs(f"{RESULTS}/domains", exist_ok=True)

print("=" * 70)
print("Full-Scale CATH Superfamily Classification via InterPro API")
print("=" * 70)

# ---- Load all proteins that need classification ----
struct = pd.read_csv(f"{RESULTS}/structures/structure_index_full.tsv", sep="\t")
acc_col = "uniprot_accession" if "uniprot_accession" in struct.columns else "accession"
all_accs = list(struct[acc_col].dropna().unique())
print(f"Total proteins to classify: {len(all_accs)}")

# ---- Load existing CATH (pilot) ----
cath_path = f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv"
existing_cath = {}
if os.path.exists(cath_path):
    cath = pd.read_csv(cath_path, sep="\t")
    cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
    for acc, group in cath.groupby(cath_acc_col):
        domains = []
        for _, row in group.iterrows():
            domains.append({
                "superfamily": row.get("cath_superfamily", ""),
                "superfamily_name": row.get("cath_superfamily_name", ""),
                "cath_class": row.get("cath_class", ""),
                "domain_start": int(row.get("domain_start", 0)),
                "domain_end": int(row.get("domain_end", 0)),
            })
        existing_cath[acc] = domains
    print(f"Existing CATH (pilot): {len(existing_cath)} proteins")

# ---- Checkpoint for resume ----
checkpoint_path = f"{RESULTS}/domains/_cath_full_checkpoint.json"
if os.path.exists(checkpoint_path):
    with open(checkpoint_path) as f:
        checkpoint = json.load(f)
    print(f"Resuming from checkpoint: {len(checkpoint.get('completed', []))} already done")
else:
    checkpoint = {"completed": [], "no_cath": [], "errors": []}

completed = set(checkpoint["completed"])
no_cath = set(checkpoint.get("no_cath", []))

# Skip proteins already processed
work = [acc for acc in all_accs if acc not in completed and acc not in existing_cath and acc not in no_cath]
print(f"Need to query: {len(work)} proteins")

# ---- Query InterPro API ----
BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/cathgene3d/protein/uniprot"

domain_records = []
# Start with existing CATH data
for acc, domains in existing_cath.items():
    for d in domains:
        domain_records.append({
            "accession": acc,
            "cath_superfamily": d["superfamily"],
            "cath_superfamily_name": d["superfamily_name"],
            "cath_class": d["cath_class"],
            "domain_start": d["domain_start"],
            "domain_end": d["domain_end"],
            "source": "CATH_pilot",
        })

success = 0
no_result = 0
errors = 0

for i, acc in enumerate(work):
    if (i + 1) % 100 == 0:
        print(f"  Progress: {i+1}/{len(work)} ({100*(i+1)/len(work):.1f}%) — success={success}, no_cath={no_result}, errors={errors}")
        # Save checkpoint
        checkpoint["completed"] = list(completed)
        checkpoint["no_cath"] = list(no_cath)
        with open(checkpoint_path, "w") as f:
            json.dump(checkpoint, f)

    url = f"{BASE_URL}/{acc}"
    try:
        resp = requests.get(url, timeout=30, headers={"Accept": "application/json"})

        if resp.status_code == 200:
            data = resp.json()
            results = data.get("results", [])
            if results:
                for entry in results:
                    sf_acc = entry.get("metadata", {}).get("accession", "")
                    sf_name = entry.get("metadata", {}).get("name", "")
                    # Extract CATH class from accession (e.g., G3DSA:3.20.20.70 -> class 3)
                    cath_class = ""
                    if ":" in sf_acc:
                        parts = sf_acc.split(":")[1].split(".")
                        if parts:
                            cath_class = parts[0]

                    # Get fragment locations
                    proteins = entry.get("proteins", [])
                    for prot in proteins:
                        locations = prot.get("entry_protein_locations", [])
                        for loc in locations:
                            for frag in loc.get("fragments", []):
                                domain_records.append({
                                    "accession": acc,
                                    "cath_superfamily": sf_acc.replace("G3DSA:", "") if "G3DSA:" in sf_acc else sf_acc,
                                    "cath_superfamily_name": sf_name,
                                    "cath_class": cath_class,
                                    "domain_start": frag.get("start", 0),
                                    "domain_end": frag.get("end", 0),
                                    "source": "InterPro_Gene3D",
                                })
                success += 1
            else:
                no_result += 1
                no_cath.add(acc)
        elif resp.status_code == 204:
            no_result += 1
            no_cath.add(acc)
        elif resp.status_code == 404:
            no_result += 1
            no_cath.add(acc)
        elif resp.status_code == 429:
            # Rate limited — wait and retry
            wait = int(resp.headers.get("Retry-After", 5))
            print(f"  Rate limited, waiting {wait}s...")
            time.sleep(wait)
            continue  # retry this accession
        else:
            errors += 1
            checkpoint.setdefault("errors", []).append(acc)

        completed.add(acc)
        time.sleep(0.35)  # ~3 requests/sec (InterPro allows up to 15/sec for registered users)

    except requests.exceptions.Timeout:
        errors += 1
        time.sleep(2)
    except Exception as e:
        errors += 1
        time.sleep(1)

# Final checkpoint
checkpoint["completed"] = list(completed)
checkpoint["no_cath"] = list(no_cath)
with open(checkpoint_path, "w") as f:
    json.dump(checkpoint, f)

print(f"\nCompleted: {success} with CATH, {no_result} without, {errors} errors")

# ---- Save full-scale CATH ----
if domain_records:
    df = pd.DataFrame(domain_records)
    out_path = f"{RESULTS}/domains/cath_domain_assignments_full.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path} ({len(df)} domain records, {df['accession'].nunique()} proteins)")

    # Summary per protein
    protein_summary = df.groupby("accession").agg(
        n_domains=("cath_superfamily", "count"),
        top_superfamily=("cath_superfamily", "first"),
        source=("source", "first"),
    ).reset_index()
    summary_path = f"{RESULTS}/domains/cath_protein_summary_full.tsv"
    protein_summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Saved: {summary_path} ({len(protein_summary)} proteins)")

    # Coverage stats
    total_queried = len(all_accs)
    with_cath = df["accession"].nunique()
    print(f"\n  CATH coverage: {with_cath}/{total_queried} ({100*with_cath/total_queried:.1f}%)")
    print(f"  No CATH: {len(no_cath)}")

print("\n" + "=" * 70)
print("Full-scale CATH classification complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "CATH classification complete — $(date)"
echo "============================================================"
