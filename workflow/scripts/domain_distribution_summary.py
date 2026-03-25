#!/usr/bin/env python3
"""
Task E5: Domain distribution summary per dataset.

Creates a comprehensive summary of domain architecture distribution
across the four datasets (groel, hsp60, mito, matrix).
"""

import pandas as pd
import numpy as np
from collections import Counter
from io import StringIO

BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

# --- Load inputs ---
print("Loading input files...")

domain_metrics = pd.read_csv(f"{BASE}/results/domains/domain_structural_metrics.tsv", sep="\t")
print(f"  Domain metrics: {len(domain_metrics)} domains")

protein_summary = pd.read_csv(f"{BASE}/results/domains/cath_protein_summary.tsv", sep="\t")
print(f"  Protein summary: {len(protein_summary)} proteins")

structure_index = pd.read_csv(f"{BASE}/results/structures/structure_index.tsv", sep="\t")
print(f"  Structure index: {len(structure_index)} proteins")

# --- Map proteins to datasets ---
# source_dataset can be comma-separated (e.g., "matrix,mito")
# Expand so each protein appears once per dataset
protein_datasets = []
for _, row in structure_index.iterrows():
    acc = row["uniprot_accession"]
    datasets = str(row["source_dataset"]).split(",") if pd.notna(row["source_dataset"]) else []
    for ds in datasets:
        ds = ds.strip()
        if ds:
            protein_datasets.append({"uniprot_accession": acc, "dataset": ds})

prot_ds = pd.DataFrame(protein_datasets)
print(f"  Protein-dataset mappings: {len(prot_ds)}")
print(f"  Datasets found: {sorted(prot_ds['dataset'].unique())}")

# Merge protein summary with dataset info
prot_with_ds = prot_ds.merge(protein_summary, on="uniprot_accession", how="left")

# Merge domain metrics with dataset info
domain_with_ds = domain_metrics.merge(prot_ds, on="uniprot_accession", how="left")

# --- CATH class mapping ---
def cath_class_label(sf):
    """Map CATH superfamily to broad class."""
    if pd.isna(sf) or sf == "":
        return "unknown"
    try:
        c = int(str(sf).split(".")[0])
    except (ValueError, IndexError):
        return "unknown"
    labels = {
        1: "mainly-alpha",
        2: "mainly-beta",
        3: "alpha-beta",
        4: "few-SS",
    }
    return labels.get(c, f"class-{c}")


# --- Compute summaries ---
datasets = sorted(prot_ds["dataset"].unique())

summary_rows = []
report_lines = []

report_lines.append("=" * 70)
report_lines.append("DOMAIN DISTRIBUTION REPORT")
report_lines.append("Antah Asti Prarambh — Task E5")
report_lines.append("=" * 70)

for ds in datasets:
    report_lines.append("")
    report_lines.append("-" * 70)
    report_lines.append(f"DATASET: {ds.upper()}")
    report_lines.append("-" * 70)

    # Subset
    ds_prots = prot_with_ds[prot_with_ds["dataset"] == ds].copy()
    ds_domains = domain_with_ds[domain_with_ds["dataset"] == ds].copy()

    n_proteins = ds_prots["uniprot_accession"].nunique()
    n_domains = len(ds_domains)
    n_with_domains = ds_prots[ds_prots["n_domains"] > 0]["uniprot_accession"].nunique()

    summary_rows.append({"dataset": ds, "metric": "n_proteins", "value": str(n_proteins)})
    summary_rows.append({"dataset": ds, "metric": "n_domains_total", "value": str(n_domains)})
    summary_rows.append({"dataset": ds, "metric": "n_proteins_with_domains", "value": str(n_with_domains)})

    report_lines.append(f"\nProteins: {n_proteins}")
    report_lines.append(f"Proteins with CATH domains: {n_with_domains} ({100*n_with_domains/max(n_proteins,1):.1f}%)")
    report_lines.append(f"Total domains: {n_domains}")
    if n_with_domains > 0:
        report_lines.append(f"Mean domains/protein (among assigned): {n_domains/n_with_domains:.2f}")

    # --- 1. Distribution of n_domains ---
    report_lines.append(f"\n  1. Domain count distribution:")
    n_dom_counts = ds_prots["n_domains"].value_counts().sort_index()
    for nd, cnt in n_dom_counts.items():
        label = "single-domain" if nd == 1 else ("multi-domain" if nd > 1 else "no-domain")
        report_lines.append(f"     n_domains={nd}: {cnt} proteins ({label})")
        summary_rows.append({"dataset": ds, "metric": f"n_domains_{nd}_count", "value": str(cnt)})

    # Single vs multi
    single = int((ds_prots["n_domains"] == 1).sum())
    multi = int((ds_prots["n_domains"] > 1).sum())
    no_dom = int((ds_prots["n_domains"] == 0).sum())
    summary_rows.append({"dataset": ds, "metric": "single_domain_proteins", "value": str(single)})
    summary_rows.append({"dataset": ds, "metric": "multi_domain_proteins", "value": str(multi)})
    summary_rows.append({"dataset": ds, "metric": "no_domain_proteins", "value": str(no_dom)})

    # --- 2. Most common CATH superfamilies (top 10) ---
    report_lines.append(f"\n  2. Top 10 CATH superfamilies:")
    if len(ds_domains) > 0:
        sf_counts = ds_domains["cath_superfamily"].value_counts().head(10)
        for rank, (sf, cnt) in enumerate(sf_counts.items(), 1):
            pct = 100 * cnt / len(ds_domains)
            report_lines.append(f"     {rank:2d}. {sf}: {cnt} domains ({pct:.1f}%)")
            summary_rows.append({"dataset": ds, "metric": f"top_sf_rank{rank}", "value": f"{sf}:{cnt}"})
    else:
        report_lines.append("     No domains assigned.")

    # --- 3. CATH class distribution ---
    report_lines.append(f"\n  3. CATH class distribution:")
    if len(ds_domains) > 0:
        ds_domains["cath_class"] = ds_domains["cath_superfamily"].apply(cath_class_label)
        class_counts = ds_domains["cath_class"].value_counts()
        for cls, cnt in class_counts.items():
            pct = 100 * cnt / len(ds_domains)
            report_lines.append(f"     {cls}: {cnt} ({pct:.1f}%)")
            summary_rows.append({"dataset": ds, "metric": f"cath_class_{cls}", "value": str(cnt)})

    # --- 4. Mean domain structural metrics ---
    report_lines.append(f"\n  4. Mean domain structural metrics:")
    metric_cols = [
        ("domain_mean_plddt", "mean_plddt"),
        ("domain_min_plddt", "min_plddt"),
        ("domain_frac_plddt_gt70", "frac_plddt>70"),
        ("domain_frac_plddt_gt90", "frac_plddt>90"),
        ("domain_frac_helix", "frac_helix"),
        ("domain_frac_strand", "frac_strand"),
        ("domain_frac_coil", "frac_coil"),
    ]
    for col, label in metric_cols:
        if col in ds_domains.columns:
            vals = ds_domains[col].dropna()
            if len(vals) > 0:
                m = vals.mean()
                s = vals.std()
                report_lines.append(f"     {label:20s}: {m:.3f} +/- {s:.3f}  (n={len(vals)})")
                summary_rows.append({"dataset": ds, "metric": f"mean_{label}", "value": f"{m:.4f}"})
                summary_rows.append({"dataset": ds, "metric": f"std_{label}", "value": f"{s:.4f}"})

    # --- 5. Domain architecture strings ---
    report_lines.append(f"\n  5. Most common domain architectures (top 10):")
    ds_arch = ds_prots[ds_prots["domain_architecture"].notna() & (ds_prots["domain_architecture"] != "")]
    if len(ds_arch) > 0:
        arch_counts = ds_arch["domain_architecture"].value_counts().head(10)
        for rank, (arch, cnt) in enumerate(arch_counts.items(), 1):
            report_lines.append(f"     {rank:2d}. [{arch}]: {cnt} proteins")
            summary_rows.append({"dataset": ds, "metric": f"top_arch_rank{rank}", "value": f"{arch}:{cnt}"})
    else:
        report_lines.append("     No architectures found.")

report_lines.append("")
report_lines.append("=" * 70)
report_lines.append("END OF REPORT")
report_lines.append("=" * 70)

# --- Save TSV summary ---
summary_df = pd.DataFrame(summary_rows)
tsv_path = f"{BASE}/results/domains/domain_distribution_summary.tsv"
summary_df.to_csv(tsv_path, sep="\t", index=False)
print(f"\nSaved: {tsv_path}")

# --- Save text report ---
report_path = f"{BASE}/results/domains/domain_distribution_report.txt"
with open(report_path, "w") as f:
    f.write("\n".join(report_lines) + "\n")
print(f"Saved: {report_path}")

# --- Print report to stdout ---
print("\n" + "\n".join(report_lines))
