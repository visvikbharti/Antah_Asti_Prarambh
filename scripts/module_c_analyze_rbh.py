#!/usr/bin/env python3
"""
Module C - Step C3 Analysis: Analyze MMseqs2 RBH results, generate summary and annotated TSV.
"""

import pandas as pd
import numpy as np
from datetime import datetime
from Bio import SeqIO

BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

# ─── Input files ──────────────────────────────────────────────────────────────
RBH_TSV     = f"{BASE}/results/homology/rbh_groel_hsp60.tsv"
GROEL_TSV   = f"{BASE}/data/processed/groel_substrates_standardized.tsv"
HSP60_TSV   = f"{BASE}/data/processed/hsp60_tier1_substrates.tsv"
GROEL_FASTA = f"{BASE}/data/processed/groel_substrates.fasta"
HSP60_FASTA = f"{BASE}/data/processed/hsp60_tier1_substrates.fasta"

# ─── Output files ─────────────────────────────────────────────────────────────
ANNOTATED_OUT = f"{BASE}/results/homology/rbh_groel_hsp60_annotated.tsv"
REPORT_OUT    = f"{BASE}/results/homology/rbh_summary_report.txt"


# ═══════════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════════

# RBH results (no header from MMseqs2)
rbh_cols = ["query", "target", "evalue", "bits", "alnlen", "qcov", "tcov", "pident", "qlen", "tlen"]
rbh = pd.read_csv(RBH_TSV, sep="\t", header=None, names=rbh_cols)

# GroEL metadata
groel = pd.read_csv(GROEL_TSV, sep="\t")
groel_meta = groel[["current_accession", "gene_symbol", "groel_class", "protein_name"]].copy()
groel_meta = groel_meta.rename(columns={"current_accession": "groel_accession",
                                         "gene_symbol": "groel_gene",
                                         "protein_name": "groel_protein_name"})

# HSP60 metadata
hsp60 = pd.read_csv(HSP60_TSV, sep="\t")
hsp60_meta = hsp60[["uniprot_id", "gene_name", "description"]].copy()
hsp60_meta = hsp60_meta.rename(columns={"uniprot_id": "hsp60_accession",
                                         "gene_name": "hsp60_gene",
                                         "description": "hsp60_description"})

# Extract accessions from full FASTA header IDs (sp|ACC|ENTRY format)
# The MMseqs2 output uses the full FASTA ID, need to extract accession
def extract_acc(fasta_id):
    """Extract accession from sp|ACC|ENTRY or just return as-is."""
    parts = str(fasta_id).split("|")
    if len(parts) >= 2:
        return parts[1]
    return str(fasta_id)

rbh["groel_accession"] = rbh["query"].apply(extract_acc)
rbh["hsp60_accession"] = rbh["target"].apply(extract_acc)


# ═══════════════════════════════════════════════════════════════════════════════
# Merge annotations
# ═══════════════════════════════════════════════════════════════════════════════

# Merge GroEL info
merged = rbh.merge(groel_meta, on="groel_accession", how="left")

# Merge HSP60 info
merged = merged.merge(hsp60_meta, on="hsp60_accession", how="left")

# Sort by bit score descending
merged = merged.sort_values("bits", ascending=False).reset_index(drop=True)

# Write annotated TSV with required columns
annotated_cols = ["groel_accession", "groel_gene", "groel_class",
                  "hsp60_accession", "hsp60_gene", "evalue", "bits",
                  "pident", "alnlen", "qcov", "tcov"]
annotated = merged[annotated_cols].copy()
annotated.to_csv(ANNOTATED_OUT, sep="\t", index=False)
print(f"Annotated TSV written to: {ANNOTATED_OUT}")
print(f"  Rows: {len(annotated)}")


# ═══════════════════════════════════════════════════════════════════════════════
# Generate summary report
# ═══════════════════════════════════════════════════════════════════════════════

lines = []
def w(s=""):
    lines.append(s)

w("=" * 78)
w("MODULE C - MMseqs2 RECIPROCAL BEST HIT (RBH) ANALYSIS REPORT")
w("=" * 78)
w(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
w()
w("Query:  GroEL substrates (E. coli)   — 252 proteins")
w("Target: HSP60 Tier 1 substrates (H. sapiens) — 266 proteins")
w("Method: MMseqs2 easy-rbh (v18.8, default parameters, e-value < 0.001)")
w()

# 1. Overall RBH counts
w("-" * 78)
w("1. OVERALL RBH RESULTS")
w("-" * 78)
n_rbh = len(rbh)
n_groel_total = 252
n_hsp60_total = 266
n_groel_matched = rbh["groel_accession"].nunique()
n_hsp60_matched = rbh["hsp60_accession"].nunique()

w(f"   Total RBH pairs found:           {n_rbh}")
w(f"   Unique GroEL proteins matched:   {n_groel_matched} / {n_groel_total} ({100*n_groel_matched/n_groel_total:.1f}%)")
w(f"   Unique HSP60 proteins matched:   {n_hsp60_matched} / {n_hsp60_total} ({100*n_hsp60_matched/n_hsp60_total:.1f}%)")
w()

# 2. E-value distribution
w("-" * 78)
w("2. E-VALUE DISTRIBUTION")
w("-" * 78)
evals = merged["evalue"]
w(f"   Min e-value:                      {evals.min():.2e}")
w(f"   Max e-value:                      {evals.max():.2e}")
w(f"   Median e-value:                   {evals.median():.2e}")
w(f"   Mean e-value:                     {evals.mean():.2e}")
w()
# E-value bins
bins = [(0, 1e-100, "< 1e-100 (extremely strong)"),
        (1e-100, 1e-50, "1e-100 to 1e-50 (very strong)"),
        (1e-50, 1e-20, "1e-50 to 1e-20 (strong)"),
        (1e-20, 1e-10, "1e-20 to 1e-10 (moderate)"),
        (1e-10, 1e-5, "1e-10 to 1e-5 (weak)"),
        (1e-5, 1e-3, "1e-5 to 1e-3 (marginal)")]
w("   E-value bin breakdown:")
for lo, hi, label in bins:
    count = ((evals >= lo) & (evals < hi)).sum()
    if count > 0:
        w(f"     {label:42s}: {count:3d} ({100*count/n_rbh:.1f}%)")
w()

# 3. Percent identity distribution
w("-" * 78)
w("3. PERCENT IDENTITY DISTRIBUTION")
w("-" * 78)
pids = merged["pident"]
w(f"   Min:     {pids.min():.1f}%")
w(f"   Max:     {pids.max():.1f}%")
w(f"   Median:  {pids.median():.1f}%")
w(f"   Mean:    {pids.mean():.1f}%")
w(f"   Std:     {pids.std():.1f}%")
w()
pid_bins = [(0, 20, "< 20%"),
            (20, 30, "20-30%"),
            (30, 40, "30-40%"),
            (40, 50, "40-50%"),
            (50, 60, "50-60%"),
            (60, 100.1, ">= 60%")]
w("   Percent identity bin breakdown:")
for lo, hi, label in pid_bins:
    count = ((pids >= lo) & (pids < hi)).sum()
    if count > 0:
        w(f"     {label:20s}: {count:3d} ({100*count/n_rbh:.1f}%)")
w()

# 4. Coverage distribution
w("-" * 78)
w("4. ALIGNMENT COVERAGE DISTRIBUTION")
w("-" * 78)
w(f"   Query coverage (qcov) - GroEL proteins:")
w(f"     Min:    {merged['qcov'].min():.3f}")
w(f"     Max:    {merged['qcov'].max():.3f}")
w(f"     Median: {merged['qcov'].median():.3f}")
w(f"     Mean:   {merged['qcov'].mean():.3f}")
w()
w(f"   Target coverage (tcov) - HSP60 proteins:")
w(f"     Min:    {merged['tcov'].min():.3f}")
w(f"     Max:    {merged['tcov'].max():.3f}")
w(f"     Median: {merged['tcov'].median():.3f}")
w(f"     Mean:   {merged['tcov'].mean():.3f}")
w()

# 5. Breakdown by GroEL class
w("-" * 78)
w("5. BREAKDOWN BY GroEL CLASS")
w("-" * 78)

# Class counts from full dataset
groel_class_counts = groel["groel_class"].value_counts()
w("   GroEL class distribution (full dataset):")
for cls in sorted(groel_class_counts.index):
    w(f"     Class {cls}: {groel_class_counts[cls]} substrates")
w()

# RBH matches per class
rbh_class_counts = merged["groel_class"].value_counts()
w("   RBH matches per GroEL class:")
for cls in sorted(rbh_class_counts.index):
    total_in_class = groel_class_counts.get(cls, 0)
    n_matched = rbh_class_counts.get(cls, 0)
    pct = 100 * n_matched / total_in_class if total_in_class > 0 else 0
    w(f"     Class {cls}: {n_matched:3d} RBH pairs ({pct:.1f}% of class)")
w()

# Class III detail
w("   Class III substrates with human HSP60 orthologs (key finding):")
class3_rbh = merged[merged["groel_class"] == "III"].copy()
n_class3_total = groel_class_counts.get("III", 0)
n_class3_rbh = len(class3_rbh)
w(f"     {n_class3_rbh} out of {n_class3_total} Class III substrates have RBH matches")
w(f"     ({100*n_class3_rbh/n_class3_total:.1f}% of Class III)")
w()
if len(class3_rbh) > 0:
    w("   Class III RBH pairs:")
    for _, row in class3_rbh.sort_values("bits", ascending=False).iterrows():
        w(f"     {row['groel_gene']:10s} ({row['groel_accession']}) <-> "
          f"{row['hsp60_gene']} ({row['hsp60_accession']})  "
          f"bits={row['bits']:.0f}  pident={row['pident']:.1f}%  evalue={row['evalue']:.1e}")
w()

# Stats per class
w("   Average alignment statistics per class:")
w(f"   {'Class':<8s} {'N':>4s} {'Avg bits':>10s} {'Avg pident':>12s} {'Avg qcov':>10s} {'Avg tcov':>10s} {'Avg evalue':>12s}")
for cls in sorted(merged["groel_class"].unique()):
    sub = merged[merged["groel_class"] == cls]
    w(f"   {cls:<8s} {len(sub):>4d} {sub['bits'].mean():>10.1f} {sub['pident'].mean():>11.1f}% {sub['qcov'].mean():>10.3f} {sub['tcov'].mean():>10.3f} {sub['evalue'].mean():>12.2e}")
w()

# 6. Top pairs by bit score
w("-" * 78)
w("6. TOP 20 RBH PAIRS BY BIT SCORE")
w("-" * 78)
w()
w(f"{'Rank':>4s}  {'GroEL Gene':<12s} {'GroEL Acc':<10s} {'Cls':>3s}  {'HSP60 Gene':<12s} {'HSP60 Acc':<12s} {'Bits':>7s} {'pident':>7s} {'E-value':>12s} {'qcov':>6s} {'tcov':>6s}")
w(f"{'----':>4s}  {'----------':<12s} {'---------':<10s} {'---':>3s}  {'----------':<12s} {'---------':<12s} {'-----':>7s} {'------':>7s} {'-------':>12s} {'----':>6s} {'----':>6s}")

top = merged.head(20)
for i, (_, row) in enumerate(top.iterrows(), 1):
    groel_g = str(row['groel_gene'])[:11]
    hsp60_g = str(row['hsp60_gene'])[:11]
    w(f"{i:>4d}  {groel_g:<12s} {row['groel_accession']:<10s} {row['groel_class']:>3s}  "
      f"{hsp60_g:<12s} {row['hsp60_accession']:<12s} {row['bits']:>7.0f} {row['pident']:>6.1f}% "
      f"{row['evalue']:>12.2e} {row['qcov']:>6.3f} {row['tcov']:>6.3f}")
w()

# 7. Biological interpretation
w("-" * 78)
w("7. BIOLOGICAL INTERPRETATION")
w("-" * 78)
w()
w(f"   Of the 252 GroEL substrates, {n_groel_matched} ({100*n_groel_matched/n_groel_total:.1f}%) have")
w(f"   reciprocal best hits among the 266 HSP60 Tier 1 substrates.")
w()
w("   This represents pairs where both:")
w("     (a) The GroEL substrate's best match in the HSP60 set is that human protein, AND")
w("     (b) That human protein's best match in the GroEL set is that E. coli protein.")
w()
w("   These are high-confidence ortholog pairs that are substrates of both the")
w("   bacterial GroEL and human mitochondrial HSP60 chaperonins.")
w()

# Class breakdown interpretation
for cls in ["I", "II", "III"]:
    n_cls = groel_class_counts.get(cls, 0)
    n_cls_rbh = rbh_class_counts.get(cls, 0) if cls in rbh_class_counts.index else 0
    pct = 100 * n_cls_rbh / n_cls if n_cls > 0 else 0
    w(f"   Class {cls}: {n_cls_rbh}/{n_cls} ({pct:.1f}%) have human HSP60 orthologs")

w()
w("   Class III proteins are the most GroEL-dependent (obligate substrates).")
if "III" in rbh_class_counts.index:
    pct3 = 100 * rbh_class_counts["III"] / groel_class_counts["III"]
    w(f"   Finding: {rbh_class_counts['III']} Class III proteins ({pct3:.1f}%) have conserved")
    w(f"   chaperonin dependence from E. coli GroEL to human HSP60.")
w()
w("=" * 78)
w("END OF REPORT")
w("=" * 78)

report_text = "\n".join(lines)
with open(REPORT_OUT, "w") as f:
    f.write(report_text)

print(f"\nSummary report written to: {REPORT_OUT}")
print()
print(report_text)
