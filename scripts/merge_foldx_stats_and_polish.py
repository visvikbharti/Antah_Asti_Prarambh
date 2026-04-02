#!/usr/bin/env python3
"""
Merge FoldX-specific statistics into the comprehensive Phase 2 test framework.

The previous Module H (pre-FoldX) ran 56 tests across 3 families.
The new Module H (with FoldX) ran 7 FoldX-specific tests.
This script:
  1. Loads both sets
  2. Adds the new FoldX tests to the old comprehensive set
  3. Re-runs hierarchical BH correction on the combined ~59 tests
  4. Generates polished publication figures

Usage:
    python3 scripts/merge_foldx_stats_and_polish.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

warnings.filterwarnings("ignore")

# Paths
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PHASE2 = os.path.join(BASE, "results", "phase2")
STATS = os.path.join(PHASE2, "stats")
FIGS = os.path.join(PHASE2, "figures")

print("=" * 70)
print("Merging FoldX statistics into comprehensive Phase 2 framework")
print("=" * 70)

# ---- Load the previous 56-test comprehensive stats ----
old_path = os.path.join(STATS, "corrected_pvalues_full_prev56.tsv")
old = pd.read_csv(old_path, sep="\t")
print(f"\nPrevious comprehensive stats: {len(old)} tests")
print(f"  Families: {old['family'].value_counts().to_dict()}")

# ---- Load the new FoldX-specific stats ----
new_path = os.path.join(STATS, "corrected_pvalues_full.tsv")
new = pd.read_csv(new_path, sep="\t")
print(f"\nNew FoldX stats: {len(new)} tests")

# ---- Identify truly new tests (not in old set) ----
old_hypotheses = set(old["hypothesis"].values)
new_only = new[~new["hypothesis"].isin(old_hypotheses)]
print(f"New tests to add: {len(new_only)}")
for _, row in new_only.iterrows():
    print(f"  {row['hypothesis']}: p={row['p_value']:.2e}")

# ---- Harmonize column names ----
# Old format: family, hypothesis, test, statistic, p_value, effect_size, ci_low, ci_high, n1, n2, direction, p_bh_within, sig_within_family, family_significant, significant_overall
# New format: family, hypothesis, test, statistic, p_value, effect_size, n_substrate, n_background, p_bh, significant

# Build unified records from new_only
new_records = []
for _, row in new_only.iterrows():
    new_records.append({
        "family": row["family"],
        "hypothesis": row["hypothesis"],
        "test": row["test"],
        "statistic": row["statistic"],
        "p_value": row["p_value"],
        "effect_size": row["effect_size"],
        "ci_low": np.nan,
        "ci_high": np.nan,
        "n1": row.get("n_substrate", np.nan),
        "n2": row.get("n_background", np.nan),
        "direction": "substrate_lower" if row["effect_size"] < 0 else "substrate_higher",
    })

if new_records:
    new_df = pd.DataFrame(new_records)
    # Merge
    combined = pd.concat([old[["family", "hypothesis", "test", "statistic", "p_value",
                               "effect_size", "ci_low", "ci_high", "n1", "n2", "direction"]],
                          new_df], ignore_index=True)
else:
    combined = old[["family", "hypothesis", "test", "statistic", "p_value",
                    "effect_size", "ci_low", "ci_high", "n1", "n2", "direction"]].copy()

print(f"\nCombined: {len(combined)} total tests")

# ---- Hierarchical BH correction ----
# Step 1: BH correction within each family
combined["p_bh_within"] = np.nan
combined["sig_within_family"] = False

for family in combined["family"].unique():
    mask = combined["family"] == family
    pvals = combined.loc[mask, "p_value"].values
    n = len(pvals)
    if n == 0:
        continue

    # BH correction
    sorted_idx = np.argsort(pvals)
    corrected = np.zeros(n)
    for rank_i, idx in enumerate(sorted_idx):
        corrected[idx] = pvals[idx] * n / (rank_i + 1)
    # Enforce monotonicity (step-up)
    for i in range(n - 2, -1, -1):
        corrected[sorted_idx[i]] = min(corrected[sorted_idx[i]], corrected[sorted_idx[i + 1]])
    corrected = np.minimum(corrected, 1.0)

    combined.loc[mask, "p_bh_within"] = corrected
    combined.loc[mask, "sig_within_family"] = corrected < 0.05

# Step 2: Family-level Simes gate
combined["family_significant"] = False
for family in combined["family"].unique():
    mask = combined["family"] == family
    min_p = combined.loc[mask, "p_value"].min()
    n_fam = mask.sum()
    # Simes: family passes if min(p) * n_family < alpha
    simes_p = min_p * n_fam
    combined.loc[mask, "family_significant"] = simes_p < 0.05

# Step 3: Overall significance = within-family AND family gate
combined["significant_overall"] = combined["sig_within_family"] & combined["family_significant"]

# ---- Summary ----
print(f"\n{'=' * 70}")
print("MERGED STATISTICAL RESULTS")
print(f"{'=' * 70}")
print(f"Total tests: {len(combined)}")
print(f"Significant within family: {combined['sig_within_family'].sum()}")
print(f"Family gates passed: {combined['family_significant'].value_counts().get(True, 0)} families")
print(f"Significant overall: {combined['significant_overall'].sum()}")

print("\nBy family:")
for family in combined["family"].unique():
    fam = combined[combined["family"] == family]
    sig = fam["significant_overall"].sum()
    gate = "PASS" if fam["family_significant"].iloc[0] else "FAIL"
    print(f"  {family}: {sig}/{len(fam)} significant (gate: {gate})")

# Print all significant results
print(f"\n{'=' * 70}")
print("SIGNIFICANT RESULTS (BH < 0.05, family gate passed)")
print(f"{'=' * 70}")
sig_results = combined[combined["significant_overall"]].sort_values("p_value")
for _, row in sig_results.iterrows():
    print(f"  {row['hypothesis']}: p={row['p_value']:.2e}, p_BH={row['p_bh_within']:.2e}, "
          f"effect={row['effect_size']:.3f}, dir={row['direction']}")

# ---- Save merged stats ----
out_path = os.path.join(STATS, "corrected_pvalues_full.tsv")
combined.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path} ({len(combined)} tests)")

# ---- Write summary report ----
report_path = os.path.join(STATS, "statistics_summary_full.txt")
with open(report_path, "w") as f:
    f.write("=" * 70 + "\n")
    f.write("ANTAH ASTI PRARAMBH — PHASE 2 FULL-SCALE STATISTICS\n")
    f.write("(Merged: comprehensive 56-test framework + FoldX DeltaG tests)\n")
    f.write(f"Generated: {pd.Timestamp.now()}\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"Total tests: {len(combined)}\n")
    f.write(f"Significant within family: {combined['sig_within_family'].sum()}\n")
    f.write(f"Significant overall: {combined['significant_overall'].sum()}\n\n")

    for family in ["domain_architecture", "stability_asymmetry", "mts_targeting"]:
        fam = combined[combined["family"] == family]
        if fam.empty:
            continue
        gate = "PASS" if fam["family_significant"].iloc[0] else "FAIL"
        f.write(f"\n{'=' * 50}\n")
        f.write(f"{family.upper()} (gate: {gate})\n")
        f.write(f"{'=' * 50}\n\n")
        for _, row in fam.iterrows():
            sig_mark = " ***" if row["significant_overall"] else ""
            f.write(f"{row['hypothesis']}:\n")
            f.write(f"  Test: {row['test']}\n")
            f.write(f"  Statistic: {row['statistic']:.4f}\n")
            f.write(f"  p-value: {row['p_value']:.2e}\n")
            f.write(f"  p-BH (within family): {row['p_bh_within']:.2e}\n")
            f.write(f"  Effect size: {row['effect_size']:.4f}\n")
            f.write(f"  Direction: {row['direction']}\n")
            f.write(f"  n1={row['n1']}, n2={row['n2']}{sig_mark}\n\n")

    f.write("\n" + "=" * 70 + "\n")
    f.write("NOTES:\n")
    f.write("  - FoldX total_energy is NOT ΔG_folding; it's total internal energy\n")
    f.write("  - FoldX was parameterized on experimental structures, not AlphaFold\n")
    f.write("  - pLDDT is model confidence, NOT thermodynamic stability\n")
    f.write("  - Contact order is the primary folding kinetics proxy\n")
    f.write("  - Hierarchical BH: within-family + Simes family gate\n")
    f.write("  - *** = significant at both levels\n")
    f.write("=" * 70 + "\n")

print(f"Saved: {report_path}")

# ---- FoldX-specific summary ----
foldx_path = os.path.join(PHASE2, "foldx", "foldx_stability_all.tsv")
if os.path.exists(foldx_path):
    foldx = pd.read_csv(foldx_path, sep="\t")
    foldx_success = foldx[foldx["status"] == "success"]
    print(f"\n{'=' * 70}")
    print("FoldX SUMMARY")
    print(f"{'=' * 70}")
    print(f"Total proteins: {len(foldx)}")
    print(f"Successful: {len(foldx_success)}")
    print(f"DeltaG: mean={foldx_success['total_energy'].mean():.2f}, "
          f"median={foldx_success['total_energy'].median():.2f}, "
          f"std={foldx_success['total_energy'].std():.2f}")
    print(f"Range: [{foldx_success['total_energy'].min():.2f}, {foldx_success['total_energy'].max():.2f}]")

print(f"\n{'=' * 70}")
print("Merge complete.")
print(f"{'=' * 70}")
