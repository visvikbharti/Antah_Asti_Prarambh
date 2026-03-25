#!/usr/bin/env python3
"""
Filter and standardize the HSP60 interactome dataset.
Produces cleaned output files with quality tiers.
"""

import pandas as pd
import numpy as np
import os
import re

# ----- Paths -----
INPUT_FILE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh/data/raw/custom/hsp60_interactome_clean.tsv"
OUTPUT_DIR = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh/data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUT_STANDARDIZED = os.path.join(OUTPUT_DIR, "hsp60_interactome_standardized.tsv")
OUT_TIER1 = os.path.join(OUTPUT_DIR, "hsp60_tier1_substrates.tsv")
OUT_REPORT = os.path.join(OUTPUT_DIR, "hsp60_filtering_report.txt")

# ----- Load data -----
df = pd.read_csv(INPUT_FILE, sep="\t")
print(f"Loaded {len(df)} proteins with {len(df.columns)} columns.")

# ----- Step a) Parse SILAC columns and handle NDIC -----
silac_cols = ["silac_ratio_coIP1", "silac_ratio_coIP2", "silac_ratio_coIP3"]

# Track NDIC counts per column
ndic_counts = {}

# Store raw string values for NDIC detection before conversion
silac_raw = {}
for col in silac_cols:
    silac_raw[col] = df[col].astype(str).str.strip()

# Count NDIC values per column
for col in silac_cols:
    ndic_mask = silac_raw[col].str.upper() == "NDIC"
    ndic_counts[col] = ndic_mask.sum()

# Convert SILAC columns to numeric (NDIC and blanks become NaN temporarily)
silac_numeric = {}
for col in silac_cols:
    silac_numeric[col] = pd.to_numeric(df[col], errors="coerce")

# Compute 95th percentile of observed numeric values per column, then imputation value = 2x that
imputation_values = {}
for col in silac_cols:
    observed = silac_numeric[col].dropna()
    p95 = observed.quantile(0.95)
    imputation_values[col] = 2.0 * p95
    print(f"  {col}: 95th pctl = {p95:.4f}, imputation value (2x) = {imputation_values[col]:.4f}")

# Create imputed columns: numeric stays, NDIC gets imputed, blank stays NaN
silac_imputed = {}
for col in silac_cols:
    imputed = silac_numeric[col].copy()
    ndic_mask = silac_raw[col].str.upper() == "NDIC"
    imputed[ndic_mask] = imputation_values[col]
    silac_imputed[col] = imputed

# ----- Step b) Compute quality metrics -----

# n_silac_replicates: count of non-NaN values in imputed columns
# (NDIC is now a number, truly missing stays NaN)
df["n_silac_replicates"] = pd.DataFrame(silac_imputed).notna().sum(axis=1)

# n_ndic_replicates: count of NDIC per protein
ndic_per_protein = pd.DataFrame({
    col: silac_raw[col].str.upper() == "NDIC" for col in silac_cols
})
df["n_ndic_replicates"] = ndic_per_protein.sum(axis=1)

# median_silac_imputed: median of the 3 imputed SILAC values
imputed_df = pd.DataFrame(silac_imputed)
df["median_silac_imputed"] = imputed_df.median(axis=1)

# has_quantitative_data: True if at least 1 numeric or NDIC value exists
df["has_quantitative_data"] = df["n_silac_replicates"] > 0

# ----- Step c) Flag proteins by category -----

# Define gene-level categories
bait_genes = {"HSPD1", "HSPE1"}
co_chaperone_genes = {"TRAP1", "HSPA9", "GRPEL1", "DNAJA3"}

# Contaminant patterns: immunoglobulins, keratins, tubulins
contaminant_patterns = re.compile(r"^(IGH|IGK|IGL|KRT|TUB)", re.IGNORECASE)

def assign_flag(row):
    gene = str(row["gene_name"]).strip()

    # Bait
    if gene in bait_genes:
        return "bait"

    # Co-chaperone
    if gene in co_chaperone_genes:
        return "co_chaperone"

    # Likely contaminant: explicit IGHG2 + pattern matching
    if gene == "IGHG2" or contaminant_patterns.match(gene):
        return "likely_contaminant"

    # No quantitative data
    if not row["has_quantitative_data"]:
        return "no_quantitative_data"

    # Low enrichment
    if pd.notna(row["median_silac_imputed"]) and row["median_silac_imputed"] < 2:
        return "low_enrichment"

    return "candidate_substrate"

df["flag"] = df.apply(assign_flag, axis=1)

# ----- Step d) Create quality tiers -----

def assign_tier(row):
    if row["flag"] != "candidate_substrate":
        if row["flag"] in ("bait", "co_chaperone", "likely_contaminant"):
            return "Excluded"
        elif row["flag"] == "no_quantitative_data":
            return "No Data"
        elif row["flag"] == "low_enrichment":
            return "Low Enrichment"
        else:
            return "Unassigned"

    is_mitocarta = str(row["mitocarta"]).strip().upper() == "X"
    median_val = row["median_silac_imputed"]

    if is_mitocarta and pd.notna(median_val) and median_val > 5:
        return "Tier 1"
    elif pd.notna(median_val) and median_val > 2:
        return "Tier 2"
    else:
        return "Tier 3"

df["quality_tier"] = df.apply(assign_tier, axis=1)

# ----- Step e) Parse existing annotations -----
# These are already in the dataframe; we just document the encoding.
# present_-HS / present_+HS: "X" = present, "o" = absent
# mitocarta: "X" = in MitoCarta
# matrix: "X" = matrix
# IMS: "X" = IMS, "o" = not

# ----- Build the filtering report -----

report_lines = []
report_lines.append("=" * 72)
report_lines.append("HSP60 INTERACTOME FILTERING REPORT")
report_lines.append("=" * 72)
report_lines.append("")

# Total proteins
report_lines.append(f"Total proteins in dataset: {len(df)}")
report_lines.append("")

# NDIC counts per SILAC column
report_lines.append("--- NDIC (Not Detected In Control) counts per SILAC column ---")
for col in silac_cols:
    report_lines.append(f"  {col}: {ndic_counts[col]} NDIC values")
    report_lines.append(f"    Imputation value (2x 95th pctl): {imputation_values[col]:.4f}")
report_lines.append("")

# Proteins per flag category
report_lines.append("--- Proteins per flag category ---")
flag_counts = df["flag"].value_counts()
for flag_name in ["bait", "co_chaperone", "likely_contaminant", "low_enrichment",
                   "no_quantitative_data", "candidate_substrate"]:
    count = flag_counts.get(flag_name, 0)
    report_lines.append(f"  {flag_name}: {count}")
report_lines.append(f"  TOTAL: {flag_counts.sum()}")
report_lines.append("")

# Proteins per quality tier
report_lines.append("--- Proteins per quality tier ---")
tier_counts = df["quality_tier"].value_counts()
for tier_name in ["Tier 1", "Tier 2", "Tier 3", "Excluded", "Low Enrichment", "No Data"]:
    count = tier_counts.get(tier_name, 0)
    report_lines.append(f"  {tier_name}: {count}")
report_lines.append(f"  TOTAL: {tier_counts.sum()}")
report_lines.append("")

# Median SILAC ratio distribution stats per tier
report_lines.append("--- Median SILAC ratio (imputed) distribution per tier ---")
for tier_name in ["Tier 1", "Tier 2", "Tier 3", "Excluded", "Low Enrichment", "No Data"]:
    subset = df[df["quality_tier"] == tier_name]["median_silac_imputed"].dropna()
    if len(subset) > 0:
        stats = subset.describe(percentiles=[0.25, 0.5, 0.75])
        report_lines.append(f"  {tier_name} (n={len(subset)}):")
        report_lines.append(f"    Min:    {stats['min']:.4f}")
        report_lines.append(f"    25th:   {stats['25%']:.4f}")
        report_lines.append(f"    Median: {stats['50%']:.4f}")
        report_lines.append(f"    75th:   {stats['75%']:.4f}")
        report_lines.append(f"    Max:    {stats['max']:.4f}")
    else:
        report_lines.append(f"  {tier_name}: no quantitative data")
report_lines.append("")

# MitoCarta breakdown per tier
report_lines.append("--- MitoCarta annotation per tier ---")
for tier_name in ["Tier 1", "Tier 2", "Tier 3", "Excluded", "Low Enrichment", "No Data"]:
    subset = df[df["quality_tier"] == tier_name]
    n_total = len(subset)
    if n_total > 0:
        n_mitocarta = (subset["mitocarta"].astype(str).str.strip().str.upper() == "X").sum()
        pct = 100.0 * n_mitocarta / n_total
        report_lines.append(f"  {tier_name}: {n_mitocarta}/{n_total} MitoCarta ({pct:.1f}%)")
    else:
        report_lines.append(f"  {tier_name}: 0 proteins")
report_lines.append("")

# Percent occupancy stats for proteins that have it
report_lines.append("--- Percent occupancy statistics (proteins with data) ---")
pct_occ = pd.to_numeric(df["percent_occupancy"], errors="coerce").dropna()
if len(pct_occ) > 0:
    stats = pct_occ.describe(percentiles=[0.25, 0.5, 0.75])
    report_lines.append(f"  N proteins with percent_occupancy: {len(pct_occ)}")
    report_lines.append(f"    Min:    {stats['min']:.4f}")
    report_lines.append(f"    25th:   {stats['25%']:.4f}")
    report_lines.append(f"    Median: {stats['50%']:.4f}")
    report_lines.append(f"    75th:   {stats['75%']:.4f}")
    report_lines.append(f"    Max:    {stats['max']:.4f}")
report_lines.append("")

# Heat stress response: how many gained/lost in +HS condition
report_lines.append("--- Heat stress response (present_-HS vs present_+HS) ---")
present_minus = df["present_-HS"].astype(str).str.strip().str.upper() == "X"
present_plus = df["present_+HS"].astype(str).str.strip().str.upper() == "X"

both = (present_minus & present_plus).sum()
only_minus = (present_minus & ~present_plus).sum()
only_plus = (~present_minus & present_plus).sum()
neither = (~present_minus & ~present_plus).sum()

report_lines.append(f"  Present in both -HS and +HS: {both}")
report_lines.append(f"  Present only in -HS (lost in HS): {only_minus}")
report_lines.append(f"  Present only in +HS (gained in HS): {only_plus}")
report_lines.append(f"  Present in neither: {neither}")
report_lines.append("")

# List Tier 1 proteins
tier1 = df[df["quality_tier"] == "Tier 1"].sort_values("median_silac_imputed", ascending=False)
report_lines.append(f"--- Tier 1 substrates (n={len(tier1)}) ---")
report_lines.append(f"  {'gene_name':<15} {'uniprot_id':<12} {'median_silac':<14} {'n_reps':<8} {'functional_category'}")
report_lines.append(f"  {'-'*15} {'-'*12} {'-'*14} {'-'*8} {'-'*30}")
for _, row in tier1.iterrows():
    gene = str(row["gene_name"])
    uid = str(row["uniprot_id"])
    med = row["median_silac_imputed"]
    nreps = int(row["n_silac_replicates"])
    func_cat = str(row["functional_category"]) if pd.notna(row["functional_category"]) else ""
    med_str = f"{med:.2f}" if pd.notna(med) else "N/A"
    report_lines.append(f"  {gene:<15} {uid:<12} {med_str:<14} {nreps:<8} {func_cat}")
report_lines.append("")
report_lines.append("=" * 72)
report_lines.append("END OF REPORT")
report_lines.append("=" * 72)

report_text = "\n".join(report_lines)

# ----- Write output files -----

# 1) Full standardized dataset
df.to_csv(OUT_STANDARDIZED, sep="\t", index=False)
print(f"\nWrote standardized dataset: {OUT_STANDARDIZED} ({len(df)} rows)")

# 2) Tier 1 substrates only
tier1_df = df[df["quality_tier"] == "Tier 1"].copy()
tier1_df.to_csv(OUT_TIER1, sep="\t", index=False)
print(f"Wrote Tier 1 substrates: {OUT_TIER1} ({len(tier1_df)} rows)")

# 3) Filtering report
with open(OUT_REPORT, "w") as f:
    f.write(report_text)
print(f"Wrote filtering report: {OUT_REPORT}")

# Print report to stdout
print("\n")
print(report_text)
