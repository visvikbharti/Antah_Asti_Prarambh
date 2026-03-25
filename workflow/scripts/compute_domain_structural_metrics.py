#!/usr/bin/env python3
"""
Task E3: Compute per-domain structural metrics from DSSP and AlphaFold CIF files.

For each CATH domain assignment, computes:
  - SS fractions (helix, strand, coil) from DSSP per-residue data
  - pLDDT statistics from AlphaFold CIF B-factors (CA atoms only)

Optimized for low memory usage:
  - DSSP loaded once into dict keyed by accession
  - CIF files parsed with lightweight text parser (no Biopython)
  - Progress printed every 200 domains
"""

import pandas as pd
import os
import sys
from collections import defaultdict

BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

# --- Load domain assignments ---
print("Loading domain assignments...")
domains = pd.read_csv(f"{BASE}/results/domains/cath_domain_assignments.tsv", sep="\t")
print(f"  Total domains: {len(domains)}")

# --- Load DSSP per-residue data into dict ---
print("Loading DSSP per-residue data...")
dssp_dict = defaultdict(dict)  # {accession: {residue_number: ss_code}}

with open(f"{BASE}/results/structures/dssp_per_residue.tsv") as f:
    header = f.readline()  # skip header
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        acc = parts[0]
        try:
            resnum = int(parts[1])
        except ValueError:
            continue
        ss_code = parts[3]
        dssp_dict[acc][resnum] = ss_code

print(f"  Proteins with DSSP data: {len(dssp_dict)}")

# --- Build CIF file path lookup ---
print("Building CIF file lookup...")
struct_index = pd.read_csv(f"{BASE}/results/structures/structure_index.tsv", sep="\t")
cif_lookup = {}
for _, row in struct_index.iterrows():
    acc = row["uniprot_accession"]
    model_path = row["model_path"]
    if pd.isna(model_path):
        continue
    # model_path is relative like "data/raw/alphafold/pilot/AF-XXX-F1-model_v6.cif"
    full_path = os.path.join(BASE, str(model_path))
    if os.path.exists(full_path):
        cif_lookup[acc] = full_path
print(f"  CIF files found: {len(cif_lookup)}")

# --- Lightweight CIF pLDDT parser ---
def parse_cif_all_ca_plddt(filepath):
    """
    Parse CIF file and extract pLDDT (B_iso_or_equiv) for ALL CA atoms.
    Returns dict {residue_number: plddt_value}.
    """
    plddt_map = {}
    in_atom_site = False

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM ") or line.startswith("HETATM "):
                in_atom_site = True
                parts = line.split()
                if len(parts) < 15:
                    continue
                if parts[3] != "CA":
                    continue
                try:
                    resnum = int(parts[8])
                    plddt = float(parts[14])
                    plddt_map[resnum] = plddt
                except (ValueError, IndexError):
                    continue
            elif in_atom_site and line.startswith("#"):
                break

    return plddt_map


# --- Process each domain ---
print("\nProcessing domains...")
results = []
n_domains = len(domains)
missing_dssp = 0
missing_cif = 0

# Cache: store full protein CA pLDDT as {accession: {resnum: plddt}}
cif_cache = {}
helix_codes = {"H", "G", "I"}
strand_codes = {"E", "B"}

for i, (_, dom) in enumerate(domains.iterrows()):
    acc = dom["uniprot_accession"]
    d_start = int(dom["domain_start"])
    d_end = int(dom["domain_end"])
    d_len = d_end - d_start + 1

    # --- DSSP metrics ---
    frac_helix = float("nan")
    frac_strand = float("nan")
    frac_coil = float("nan")

    if acc in dssp_dict:
        residues = dssp_dict[acc]
        n_helix = 0
        n_strand = 0
        n_coil = 0
        n_total = 0

        for resnum in range(d_start, d_end + 1):
            ss = residues.get(resnum, "-")
            n_total += 1
            if ss in helix_codes:
                n_helix += 1
            elif ss in strand_codes:
                n_strand += 1
            else:
                n_coil += 1

        if n_total > 0:
            frac_helix = round(n_helix / n_total, 4)
            frac_strand = round(n_strand / n_total, 4)
            frac_coil = round(n_coil / n_total, 4)
    else:
        missing_dssp += 1

    # --- pLDDT metrics ---
    mean_plddt = float("nan")
    min_plddt = float("nan")
    frac_gt70 = float("nan")
    frac_gt90 = float("nan")

    if acc in cif_lookup:
        # Parse full protein once, cache for reuse across domains
        if acc not in cif_cache:
            cif_cache[acc] = parse_cif_all_ca_plddt(cif_lookup[acc])
            # Limit cache to ~500 proteins to avoid memory issues
            if len(cif_cache) > 500:
                keys_to_remove = list(cif_cache.keys())[:100]
                for k in keys_to_remove:
                    del cif_cache[k]

        plddt_map = cif_cache[acc]
        plddt_vals = [plddt_map[r] for r in range(d_start, d_end + 1) if r in plddt_map]

        if plddt_vals:
            n = len(plddt_vals)
            mean_plddt = round(sum(plddt_vals) / n, 2)
            min_plddt = round(min(plddt_vals), 2)
            frac_gt70 = round(sum(1 for v in plddt_vals if v > 70) / n, 4)
            frac_gt90 = round(sum(1 for v in plddt_vals if v > 90) / n, 4)
    else:
        missing_cif += 1

    results.append({
        "uniprot_accession": acc,
        "domain_index": dom["domain_index"],
        "domain_start": d_start,
        "domain_end": d_end,
        "domain_length": d_len,
        "cath_superfamily": dom["cath_superfamily"],
        "domain_frac_helix": frac_helix,
        "domain_frac_strand": frac_strand,
        "domain_frac_coil": frac_coil,
        "domain_mean_plddt": mean_plddt,
        "domain_min_plddt": min_plddt,
        "domain_frac_plddt_gt70": frac_gt70,
        "domain_frac_plddt_gt90": frac_gt90,
    })

    if (i + 1) % 200 == 0 or (i + 1) == n_domains:
        print(f"  Processed {i + 1}/{n_domains} domains...")

# Clear cache
cif_cache.clear()

# --- Build output DataFrame ---
out_df = pd.DataFrame(results)

# --- Save ---
out_path = f"{BASE}/results/domains/domain_structural_metrics.tsv"
out_df.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path}")

# --- Summary ---
print("\n" + "=" * 60)
print("DOMAIN STRUCTURAL METRICS SUMMARY")
print("=" * 60)
print(f"Total domains processed: {len(out_df)}")
print(f"Missing DSSP data: {missing_dssp}")
print(f"Missing CIF files: {missing_cif}")

valid_plddt = out_df["domain_mean_plddt"].dropna()
valid_helix = out_df["domain_frac_helix"].dropna()

print(f"\nDomains with pLDDT data: {len(valid_plddt)}")
if len(valid_plddt) > 0:
    print(f"  Mean domain pLDDT: {valid_plddt.mean():.1f}")
    print(f"  Median domain pLDDT: {valid_plddt.median():.1f}")
    print(f"  Mean frac_plddt>70: {out_df['domain_frac_plddt_gt70'].dropna().mean():.3f}")
    print(f"  Mean frac_plddt>90: {out_df['domain_frac_plddt_gt90'].dropna().mean():.3f}")

print(f"\nDomains with DSSP data: {len(valid_helix)}")
if len(valid_helix) > 0:
    print(f"  Mean frac_helix: {valid_helix.mean():.3f}")
    print(f"  Mean frac_strand: {out_df['domain_frac_strand'].dropna().mean():.3f}")
    print(f"  Mean frac_coil: {out_df['domain_frac_coil'].dropna().mean():.3f}")
print("=" * 60)
