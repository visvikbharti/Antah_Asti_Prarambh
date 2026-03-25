#!/usr/bin/env python3
"""
D4: Validate structure quality — flag low-confidence AlphaFold models.

Quality tiers based on mean pLDDT:
  - Very High (>90): Experimental-quality backbone + sidechain confidence
  - High (80-90): Good backbone confidence, most sidechains reliable
  - Moderate (70-80): Backbone mostly reliable, sidechains less certain
  - Low (50-70): Significant disorder or low-confidence regions
  - Very Low (<50): Mostly disordered / unreliable model

Also flags:
  - Fraction of residues with pLDDT > 70 (useful vs total)
  - Per-dataset quality breakdown
  - Proteins with quality concerns for downstream analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")

# Load data
struct_idx = pd.read_csv(BASE / "results/structures/structure_index.tsv", sep="\t")
dssp = pd.read_csv(BASE / "results/structures/dssp_summary.tsv", sep="\t")

# Only structures that exist
has_struct = struct_idx[struct_idx["has_structure"] == True].copy()
no_struct = struct_idx[struct_idx["has_structure"] == False]

print(f"Total proteins in index: {len(struct_idx)}")
print(f"  With structure: {len(has_struct)}")
print(f"  Without structure: {len(no_struct)}")
if len(no_struct) > 0:
    print(f"  Missing: {', '.join(no_struct['uniprot_accession'].tolist())}")

# Assign quality tiers
def assign_tier(plddt):
    if plddt >= 90:
        return "very_high"
    elif plddt >= 80:
        return "high"
    elif plddt >= 70:
        return "moderate"
    elif plddt >= 50:
        return "low"
    else:
        return "very_low"

has_struct["quality_tier"] = has_struct["mean_plddt"].apply(assign_tier)

# Summary
tier_counts = has_struct["quality_tier"].value_counts()
tier_order = ["very_high", "high", "moderate", "low", "very_low"]
print("\n=== Quality Tier Distribution ===")
for tier in tier_order:
    n = tier_counts.get(tier, 0)
    pct = 100.0 * n / len(has_struct)
    print(f"  {tier:12s}: {n:5d} ({pct:5.1f}%)")

# Per-dataset breakdown
print("\n=== Per-Dataset Quality ===")
datasets_all = set()
for ds_str in has_struct["source_dataset"]:
    for ds in str(ds_str).split(","):
        datasets_all.add(ds.strip())

rows = []
for ds in sorted(datasets_all):
    mask = has_struct["source_dataset"].str.contains(ds, na=False)
    subset = has_struct[mask]
    row = {
        "dataset": ds,
        "n_proteins": len(subset),
        "mean_plddt": subset["mean_plddt"].mean(),
        "median_plddt": subset["mean_plddt"].median(),
        "frac_gt70_mean": subset["fraction_plddt_gt70"].mean(),
        "n_very_high": (subset["quality_tier"] == "very_high").sum(),
        "n_high": (subset["quality_tier"] == "high").sum(),
        "n_moderate": (subset["quality_tier"] == "moderate").sum(),
        "n_low": (subset["quality_tier"] == "low").sum(),
        "n_very_low": (subset["quality_tier"] == "very_low").sum(),
    }
    rows.append(row)
    print(f"\n  {ds} (n={row['n_proteins']}):")
    print(f"    Mean pLDDT: {row['mean_plddt']:.1f}, Median: {row['median_plddt']:.1f}")
    print(f"    Mean frac>70: {row['frac_gt70_mean']:.3f}")
    print(f"    Tiers: VH={row['n_very_high']}, H={row['n_high']}, M={row['n_moderate']}, L={row['n_low']}, VL={row['n_very_low']}")

ds_df = pd.DataFrame(rows)

# Flag proteins with quality concerns
# Criteria for flagging:
# 1. mean pLDDT < 50: model is mostly unreliable
# 2. mean pLDDT < 70 AND fraction_plddt_gt70 < 0.5: majority of residues unreliable
# 3. fraction_plddt_gt70 < 0.3: very few usable residues

has_struct["flag_very_low_plddt"] = has_struct["mean_plddt"] < 50
has_struct["flag_majority_unreliable"] = (has_struct["mean_plddt"] < 70) & (has_struct["fraction_plddt_gt70"] < 0.5)
has_struct["flag_few_usable"] = has_struct["fraction_plddt_gt70"] < 0.3
has_struct["any_quality_flag"] = has_struct["flag_very_low_plddt"] | has_struct["flag_majority_unreliable"] | has_struct["flag_few_usable"]

n_flagged = has_struct["any_quality_flag"].sum()
print(f"\n=== Quality Flags ===")
print(f"  Very low pLDDT (<50): {has_struct['flag_very_low_plddt'].sum()}")
print(f"  Majority unreliable (pLDDT<70 & frac>70 <0.5): {has_struct['flag_majority_unreliable'].sum()}")
print(f"  Few usable residues (frac>70 <0.3): {has_struct['flag_few_usable'].sum()}")
print(f"  ANY flag: {n_flagged} ({100.0*n_flagged/len(has_struct):.1f}%)")

# Detail flagged proteins
flagged = has_struct[has_struct["any_quality_flag"]].sort_values("mean_plddt")
if len(flagged) > 0:
    print(f"\n=== Flagged Proteins (n={len(flagged)}) ===")
    for _, row in flagged.iterrows():
        flags = []
        if row["flag_very_low_plddt"]:
            flags.append("VL-pLDDT")
        if row["flag_majority_unreliable"]:
            flags.append("maj-unreliable")
        if row["flag_few_usable"]:
            flags.append("few-usable")
        print(f"  {row['uniprot_accession']:12s}  pLDDT={row['mean_plddt']:5.1f}  frac>70={row['fraction_plddt_gt70']:.3f}  "
              f"residues={int(row['residues_modeled']):5d}  datasets={row['source_dataset']}  flags={','.join(flags)}")

# Impact assessment: how many flagged proteins are in substrate datasets?
print(f"\n=== Impact on Substrate Datasets ===")
for ds in ["groel", "hsp60"]:
    mask = flagged["source_dataset"].str.contains(ds, na=False)
    n = mask.sum()
    total_mask = has_struct["source_dataset"].str.contains(ds, na=False)
    total = total_mask.sum()
    print(f"  {ds}: {n}/{total} flagged ({100.0*n/total:.1f}%)")

# Cross-reference with DSSP: do flagged proteins have unusual SS?
if len(dssp) > 0 and len(flagged) > 0:
    flagged_dssp = dssp[dssp["uniprot_accession"].isin(flagged["uniprot_accession"])]
    if len(flagged_dssp) > 0:
        print(f"\n=== Flagged Proteins: Secondary Structure ===")
        print(f"  Mean helix: {flagged_dssp['frac_helix'].mean()*100:.1f}%")
        print(f"  Mean strand: {flagged_dssp['frac_strand'].mean()*100:.1f}%")
        print(f"  Mean coil: {flagged_dssp['frac_coil'].mean()*100:.1f}%")
        # Compare to unflagged
        unflagged_dssp = dssp[~dssp["uniprot_accession"].isin(flagged["uniprot_accession"])]
        print(f"  (vs unflagged: helix={unflagged_dssp['frac_helix'].mean()*100:.1f}%, "
              f"strand={unflagged_dssp['frac_strand'].mean()*100:.1f}%, "
              f"coil={unflagged_dssp['frac_coil'].mean()*100:.1f}%)")

# Save outputs
out_dir = BASE / "results/structures"

# Full quality index
quality_cols = ["uniprot_accession", "source_dataset", "residues_modeled", "mean_plddt",
                "min_plddt", "fraction_plddt_gt70", "quality_tier",
                "flag_very_low_plddt", "flag_majority_unreliable", "flag_few_usable", "any_quality_flag"]
has_struct[quality_cols].to_csv(out_dir / "structure_quality_validation.tsv", sep="\t", index=False)

# Flagged proteins only
if len(flagged) > 0:
    flagged[quality_cols].to_csv(out_dir / "flagged_low_quality_structures.tsv", sep="\t", index=False)

# Per-dataset summary
ds_df.to_csv(out_dir / "quality_per_dataset.tsv", sep="\t", index=False)

# Write report
report_lines = []
report_lines.append("=" * 70)
report_lines.append("D4: Structure Quality Validation Report")
report_lines.append("=" * 70)
report_lines.append("")
report_lines.append(f"Total structures validated: {len(has_struct)}")
report_lines.append(f"Missing structures: {len(no_struct)} ({', '.join(no_struct['uniprot_accession'].tolist()) if len(no_struct) > 0 else 'none'})")
report_lines.append("")
report_lines.append("Quality Tier Distribution:")
report_lines.append(f"  {'Tier':12s}  {'Count':>5s}  {'Percent':>7s}  {'pLDDT Range':>15s}")
tier_ranges = {"very_high": ">=90", "high": "80-90", "moderate": "70-80", "low": "50-70", "very_low": "<50"}
for tier in tier_order:
    n = tier_counts.get(tier, 0)
    report_lines.append(f"  {tier:12s}  {n:5d}  {100.0*n/len(has_struct):6.1f}%  {tier_ranges[tier]:>15s}")
report_lines.append("")
report_lines.append(f"Proteins with quality flags: {n_flagged} ({100.0*n_flagged/len(has_struct):.1f}%)")
report_lines.append("  Flag criteria:")
report_lines.append("    - Very low pLDDT: mean pLDDT < 50 (model mostly unreliable)")
report_lines.append("    - Majority unreliable: mean pLDDT < 70 AND >50% residues below 70")
report_lines.append("    - Few usable: <30% of residues have pLDDT > 70")
report_lines.append("")
report_lines.append("Recommendation:")
if n_flagged <= 20:
    report_lines.append(f"  {n_flagged} flagged proteins is a small fraction ({100.0*n_flagged/len(has_struct):.1f}%).")
    report_lines.append("  These can be retained in analyses with quality_tier as a covariate,")
    report_lines.append("  or excluded in sensitivity analyses to confirm results are robust.")
else:
    report_lines.append(f"  {n_flagged} flagged proteins — consider sensitivity analysis excluding these.")
report_lines.append("")
report_lines.append("Output files:")
report_lines.append(f"  results/structures/structure_quality_validation.tsv ({len(has_struct)} proteins)")
report_lines.append(f"  results/structures/flagged_low_quality_structures.tsv ({n_flagged} proteins)")
report_lines.append(f"  results/structures/quality_per_dataset.tsv ({len(ds_df)} datasets)")
report_lines.append("=" * 70)

report_text = "\n".join(report_lines)
(out_dir / "quality_validation_report.txt").write_text(report_text)
print(f"\n{report_text}")

print(f"\nDone. Outputs written to {out_dir}/")
