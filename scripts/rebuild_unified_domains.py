#!/usr/bin/env python3
"""
B1: Rebuild unified domain assignments with full-scale CATH data.

CATH (InterPro Gene3D) takes priority over Chainsaw when both available.
Produces updated unified_domain_assignments_full.tsv with proper source tracking.
"""

import os
import re
import pandas as pd
import numpy as np

PROJECT_DIR = os.path.expanduser("~/Downloads/Antah_Asti_Prarambh")
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Rebuild Unified Domain Assignments (Full-Scale CATH + Chainsaw)")
print("=" * 70)

# ---- Load full-scale CATH domain assignments ----
cath_path = f"{RESULTS}/domains/cath_domain_assignments_full.tsv"
cath = pd.read_csv(cath_path, sep="\t")
print(f"CATH full-scale: {len(cath)} domain records")

# Get unique CATH proteins (per-protein summary)
cath_summary_path = f"{RESULTS}/domains/cath_protein_summary_full.tsv"
cath_summary = pd.read_csv(cath_summary_path, sep="\t")
cath_proteins = set(cath_summary["accession"].values)
print(f"CATH proteins: {len(cath_proteins)}")

# ---- Load Chainsaw predictions ----
chainsaw_path = f"{RESULTS}/domains/chainsaw_full_predictions.tsv"
chainsaw = pd.read_csv(chainsaw_path, sep="\t")
print(f"Chainsaw: {len(chainsaw)} protein records")

# Extract accession from chain_id (AF-{ACC}-F1-model_v6)
def extract_acc(chain_id):
    m = re.search(r'AF-([A-Z0-9]+)-F', str(chain_id))
    return m.group(1) if m else None

chainsaw["accession"] = chainsaw["chain_id"].apply(extract_acc)
chainsaw = chainsaw.dropna(subset=["accession"])
chainsaw_proteins = set(chainsaw["accession"].values)
print(f"Chainsaw proteins: {len(chainsaw_proteins)}")

# ---- Load substrate / dataset membership lists ----
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
mito = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv", sep="\t")
matrix = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv", sep="\t")

groel_acc_col = "current_accession" if "current_accession" in groel.columns else "accession"
hsp60_acc_col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
mito_acc_col = "accession" if "accession" in mito.columns else mito.columns[0]
matrix_acc_col = "accession" if "accession" in matrix.columns else matrix.columns[0]

groel_set = set(groel[groel_acc_col].dropna())
hsp60_set = set(hsp60[hsp60_acc_col].dropna())
mito_set = set(mito[mito_acc_col].dropna())
matrix_set = set(matrix[matrix_acc_col].dropna())

def get_datasets(acc):
    ds = []
    if acc in groel_set:
        ds.append("groel")
    if acc in hsp60_set:
        ds.append("hsp60")
    if acc in matrix_set:
        ds.append("matrix_bg")
    if acc in mito_set:
        ds.append("mito_bg")
    if not ds:
        ds.append("proteome_bg")
    return ",".join(ds)

# ---- Build unified assignments: CATH priority, Chainsaw fallback ----
records = []

# Process CATH proteins first (they have authoritative domain assignments)
for acc in sorted(cath_proteins):
    acc_cath = cath_summary[cath_summary["accession"] == acc].iloc[0]
    n_domains = int(acc_cath["n_domains"])
    top_sf = acc_cath.get("top_superfamily", "")
    source = acc_cath.get("source", "CATH")

    # Get CATH class from domain-level data
    acc_domains = cath[cath["accession"] == acc]
    cath_class = ""
    if "cath_class" in acc_domains.columns and len(acc_domains) > 0:
        classes = acc_domains["cath_class"].dropna().astype(str).values
        if len(classes) > 0:
            cath_class = classes[0]

    # Superfamily list
    if "cath_superfamily" in acc_domains.columns:
        sfs = acc_domains["cath_superfamily"].dropna().values
        top_sf = "|".join(str(s) for s in sfs) if len(sfs) > 0 else top_sf

    records.append({
        "accession": acc,
        "source": source,
        "n_domains": n_domains,
        "top_superfamily": top_sf,
        "cath_class": cath_class,
        "datasets": get_datasets(acc),
    })

# Chainsaw fallback for proteins not in CATH
chainsaw_only = chainsaw_proteins - cath_proteins
print(f"Chainsaw-only proteins (no CATH): {len(chainsaw_only)}")

for acc in sorted(chainsaw_only):
    row = chainsaw[chainsaw["accession"] == acc].iloc[0]
    n_domains = int(row.get("ndom", 1))

    records.append({
        "accession": acc,
        "source": "Chainsaw",
        "n_domains": n_domains,
        "top_superfamily": "",
        "cath_class": "",
        "datasets": get_datasets(acc),
    })

# ---- Create DataFrame and save ----
unified = pd.DataFrame(records)
unified = unified.sort_values("accession").reset_index(drop=True)

out_path = f"{RESULTS}/domains/unified_domain_assignments_full.tsv"
# Back up old file
if os.path.exists(out_path):
    backup = out_path.replace(".tsv", "_old_pilot_cath.tsv")
    os.rename(out_path, backup)
    print(f"Backed up old unified to: {backup}")

unified.to_csv(out_path, sep="\t", index=False)

# ---- Statistics ----
print(f"\n{'=' * 70}")
print(f"Unified Domain Assignments Rebuilt")
print(f"{'=' * 70}")
print(f"Total proteins: {len(unified)}")
print(f"  CATH (pilot): {(unified['source'] == 'CATH_pilot').sum()}")
print(f"  CATH (InterPro): {(unified['source'] == 'InterPro_Gene3D').sum()}")
print(f"  Chainsaw: {(unified['source'] == 'Chainsaw').sum()}")
print(f"\nCoverage by dataset:")
for ds in ["groel", "hsp60", "matrix_bg", "mito_bg", "proteome_bg"]:
    count = unified["datasets"].str.contains(ds).sum()
    cath_count = unified[unified["datasets"].str.contains(ds)]["source"].isin(
        ["CATH_pilot", "InterPro_Gene3D"]).sum()
    print(f"  {ds}: {count} proteins ({cath_count} with CATH, "
          f"{count - cath_count} Chainsaw-only)")
print(f"\nMulti-domain: {(unified['n_domains'] >= 2).sum()} ({100*(unified['n_domains'] >= 2).mean():.1f}%)")
print(f"\nSaved: {out_path}")
