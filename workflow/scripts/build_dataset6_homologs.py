#!/usr/bin/env python3
"""
Task C6: Build final Dataset 6 — definitive 2-way homolog table.

Merges evidence from RBH (40 pairs) and orthogroup analysis (62 substrate pairs)
into a single unified homolog table.
"""

import pandas as pd
import sys
import os

BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

# --- Load inputs ---
print("Loading input files...")

rbh = pd.read_csv(f"{BASE}/results/homology/rbh_groel_hsp60_annotated.tsv", sep="\t")
print(f"  RBH pairs: {len(rbh)}")

ortho = pd.read_csv(f"{BASE}/results/homology/substrate_orthogroups.tsv", sep="\t")
print(f"  Orthogroups: {len(ortho)}")

groel_sub = pd.read_csv(f"{BASE}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60_sub = pd.read_csv(f"{BASE}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")

# Get sets of substrate accessions for filtering
groel_accessions = set(groel_sub["current_accession"].dropna().unique())
hsp60_accessions = set(hsp60_sub["uniprot_id"].dropna().unique())
print(f"  GroEL substrates: {len(groel_accessions)}")
print(f"  Hsp60 substrates: {len(hsp60_accessions)}")

# --- Process RBH pairs ---
print("\nProcessing RBH pairs...")
rbh_pairs = rbh[["groel_accession", "groel_gene", "groel_class",
                  "hsp60_accession", "hsp60_gene",
                  "evalue", "pident", "bits"]].copy()
rbh_pairs["evidence"] = "rbh"
rbh_pairs["orthogroup_id"] = ""

# Create a set of (groel, hsp60) pairs for quick lookup
rbh_pair_set = set(zip(rbh_pairs["groel_accession"], rbh_pairs["hsp60_accession"]))
print(f"  Unique RBH pairs: {len(rbh_pair_set)}")

# --- Process orthogroup pairs ---
print("\nProcessing orthogroup pairs...")
# Only use 'shared' orthogroups (those that have both groel and hsp60 members)
shared_ortho = ortho[ortho["category"] == "shared"].copy()
print(f"  Shared orthogroups: {len(shared_ortho)}")

ortho_rows = []
for _, row in shared_ortho.iterrows():
    og_id = row["orthogroup_id"]

    # Parse comma-separated accessions
    groel_accs = str(row["groel_accessions"]).split(",") if pd.notna(row["groel_accessions"]) else []
    hsp60_accs = str(row["hsp60_accessions"]).split(",") if pd.notna(row["hsp60_accessions"]) else []
    groel_genes = str(row["groel_genes"]).split(",") if pd.notna(row["groel_genes"]) else []
    hsp60_genes = str(row["hsp60_genes"]).split(",") if pd.notna(row["hsp60_genes"]) else []
    groel_classes = str(row["groel_classes"]).split(",") if pd.notna(row["groel_classes"]) else []

    # Filter to actual substrates
    groel_accs_filtered = [a.strip() for a in groel_accs if a.strip() in groel_accessions]
    hsp60_accs_filtered = [a.strip() for a in hsp60_accs if a.strip() in hsp60_accessions]

    # Build gene/class lookup from the original comma-separated lists
    acc_gene_map = {}
    acc_class_map = {}
    for i, a in enumerate(groel_accs):
        a = a.strip()
        if i < len(groel_genes):
            acc_gene_map[a] = groel_genes[i].strip()
        if i < len(groel_classes):
            acc_class_map[a] = groel_classes[i].strip()

    hsp60_gene_map = {}
    for i, a in enumerate(hsp60_accs):
        a = a.strip()
        if i < len(hsp60_genes):
            hsp60_gene_map[a] = hsp60_genes[i].strip()

    # Create all pairwise combinations of substrate accessions
    for g_acc in groel_accs_filtered:
        for h_acc in hsp60_accs_filtered:
            ortho_rows.append({
                "groel_accession": g_acc,
                "groel_gene": acc_gene_map.get(g_acc, ""),
                "groel_class": acc_class_map.get(g_acc, ""),
                "hsp60_accession": h_acc,
                "hsp60_gene": hsp60_gene_map.get(h_acc, ""),
                "orthogroup_id": og_id,
                "evidence": "orthogroup",
                "evalue": row["best_evalue"],
                "pident": row["mean_pident"],
                "bits": float("nan"),  # orthogroups don't have per-pair bits
            })

ortho_df = pd.DataFrame(ortho_rows)
print(f"  Orthogroup substrate pairs: {len(ortho_df)}")

# --- Merge: identify pairs found by both methods ---
print("\nMerging evidence...")
if len(ortho_df) > 0:
    ortho_pair_set = set(zip(ortho_df["groel_accession"], ortho_df["hsp60_accession"]))
else:
    ortho_pair_set = set()

both_pairs = rbh_pair_set & ortho_pair_set
print(f"  Pairs found by both methods: {len(both_pairs)}")

# Mark RBH pairs that are also in orthogroups as "both"
# and add orthogroup_id from ortho_df
ortho_lookup = {}
if len(ortho_df) > 0:
    for _, row in ortho_df.iterrows():
        key = (row["groel_accession"], row["hsp60_accession"])
        ortho_lookup[key] = row["orthogroup_id"]

for idx, row in rbh_pairs.iterrows():
    key = (row["groel_accession"], row["hsp60_accession"])
    if key in both_pairs:
        rbh_pairs.at[idx, "evidence"] = "both"
        rbh_pairs.at[idx, "orthogroup_id"] = ortho_lookup.get(key, "")

# For orthogroup-only pairs (not in RBH), keep them as "orthogroup"
if len(ortho_df) > 0:
    ortho_only = ortho_df[
        ortho_df.apply(lambda r: (r["groel_accession"], r["hsp60_accession"]) not in rbh_pair_set, axis=1)
    ].copy()
else:
    ortho_only = pd.DataFrame()

print(f"  RBH-only pairs: {len(rbh_pairs[rbh_pairs['evidence'] == 'rbh'])}")
print(f"  Both pairs: {len(rbh_pairs[rbh_pairs['evidence'] == 'both'])}")
print(f"  Orthogroup-only pairs: {len(ortho_only)}")

# --- Combine into final table ---
final = pd.concat([rbh_pairs, ortho_only], ignore_index=True)

# Ensure column order
final = final[["groel_accession", "groel_gene", "groel_class",
               "hsp60_accession", "hsp60_gene", "orthogroup_id",
               "evidence", "evalue", "pident", "bits"]]

# Sort by evidence priority (both > rbh > orthogroup), then by evalue
evidence_order = {"both": 0, "rbh": 1, "orthogroup": 2}
final["_sort"] = final["evidence"].map(evidence_order)
final = final.sort_values(["_sort", "evalue"]).drop(columns=["_sort"]).reset_index(drop=True)

# --- Save ---
out_path = f"{BASE}/data/processed/groel_hsp60_homologs.tsv"
final.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path}")

# --- Summary ---
print("\n" + "=" * 60)
print("DATASET 6 SUMMARY: GroEL-Hsp60 Homolog Pairs")
print("=" * 60)
print(f"\nTotal unique pairs: {len(final)}")
print(f"\nBy evidence type:")
for ev, count in final["evidence"].value_counts().items():
    print(f"  {ev:15s}: {count:4d}")

print(f"\nBy GroEL class:")
for cls, count in final["groel_class"].value_counts().sort_index().items():
    print(f"  Class {cls:5s}: {count:4d}")

print(f"\nUnique GroEL substrates in pairs: {final['groel_accession'].nunique()}")
print(f"Unique Hsp60 substrates in pairs: {final['hsp60_accession'].nunique()}")

# Quick stats
print(f"\nAlignment quality (RBH + both pairs):")
rbh_both = final[final["evidence"].isin(["rbh", "both"])]
print(f"  Mean pident: {rbh_both['pident'].mean():.1f}%")
print(f"  Median evalue: {rbh_both['evalue'].median():.2e}")
print(f"  Mean bits: {rbh_both['bits'].mean():.1f}")
print("=" * 60)
