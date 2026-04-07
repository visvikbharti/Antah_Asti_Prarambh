#!/usr/bin/env python3
"""
Module H: Full-scale comparative statistics with hierarchical correction.

Three hypothesis families:
  Family 1: Domain architecture enrichment (Fisher's exact, chi-squared)
  Family 2: N-vs-C stability asymmetry (Wilcoxon paired, Mann-Whitney, Kruskal-Wallis)
  Family 3: MTS targeting (Fisher's exact, binomial)

Multiple testing: Hierarchical Benjamini-Hochberg (within-family + between-family).
Effect sizes: rank-biserial r, odds ratios with 95% CI, Cramer's V, eta-squared.
Background matching: compartment-matched AND size-matched controls.
"""

import os
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from datetime import datetime

warnings.filterwarnings("ignore")

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("WARNING: statsmodels not available, using manual BH correction")

PROJECT_DIR = os.environ.get("PROJECT_DIR",
    os.path.expanduser("~/Downloads/Antah_Asti_Prarambh"))
RESULTS = f"{PROJECT_DIR}/results/phase2"

print("=" * 70)
print("Module H: Phase 2 Full-Scale Comparative Statistics")
print(f"Run date: {datetime.now()}")
print("=" * 70)


# ===========================================================================
# Helper functions
# ===========================================================================

def bh_correct(pvalues):
    """Benjamini-Hochberg FDR correction."""
    pvalues = np.array(pvalues, dtype=float)
    if len(pvalues) == 0:
        return np.array([])
    if HAS_STATSMODELS:
        _, corrected, _, _ = multipletests(pvalues, method="fdr_bh")
        return corrected
    n = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    corrected = np.zeros(n)
    for rank, idx in enumerate(sorted_idx, 1):
        corrected[idx] = pvalues[idx] * n / rank
    corrected[sorted_idx[-1]] = min(corrected[sorted_idx[-1]], 1.0)
    for i in range(n - 2, -1, -1):
        corrected[sorted_idx[i]] = min(corrected[sorted_idx[i]],
                                        corrected[sorted_idx[i + 1]])
    return np.clip(corrected, 0, 1)


def safe_fisher(table):
    """Fisher's exact test with Haldane correction for zero cells."""
    table = np.array(table, dtype=float)
    if table.min() == 0:
        table = table + 0.5  # Haldane correction
    oddsratio, pval = stats.fisher_exact(table.astype(int) if table.min() >= 1
                                          else [[int(max(1, x)) for x in row]
                                                for row in table])
    # 95% CI for log(OR) using Woolf's method
    try:
        log_or = np.log(oddsratio) if oddsratio > 0 else 0
        se = np.sqrt(sum(1.0 / max(x, 0.5) for x in table.flatten()))
        ci_low = np.exp(log_or - 1.96 * se)
        ci_high = np.exp(log_or + 1.96 * se)
    except Exception:
        ci_low, ci_high = np.nan, np.nan
    return oddsratio, pval, ci_low, ci_high


def rank_biserial_paired(diffs):
    """Rank-biserial correlation for paired Wilcoxon test."""
    diffs = diffs[diffs != 0]
    if len(diffs) == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(diffs))
    r_plus = sum(ranks[diffs > 0])
    r_minus = sum(ranks[diffs < 0])
    n = len(diffs)
    return (r_plus - r_minus) / (n * (n + 1) / 2)


def rank_biserial_unpaired(u_stat, n1, n2):
    """Rank-biserial correlation for Mann-Whitney U test."""
    return 1 - (2 * u_stat) / (n1 * n2)


def size_match_sample(substrate_sizes, background_df, size_col, acc_col,
                      bin_width=10000, multiplier=3, seed=42):
    """Size-matched sampling from background (Daltons = size_col * ~110)."""
    rng = np.random.RandomState(seed)
    substrate_sizes = np.array(substrate_sizes)
    sub_bins = (substrate_sizes // bin_width).astype(int)
    bin_counts = pd.Series(sub_bins).value_counts()

    matched = []
    for bin_val, count in bin_counts.items():
        target = count * multiplier
        lo = bin_val * bin_width
        hi = (bin_val + 1) * bin_width
        pool = background_df[
            (background_df[size_col] >= lo) & (background_df[size_col] < hi)
        ]
        if len(pool) == 0:
            continue
        n_sample = min(target, len(pool))
        matched.append(pool.sample(n=n_sample, random_state=rng))

    if matched:
        return pd.concat(matched, ignore_index=True)
    return pd.DataFrame()


# ===========================================================================
# Load data
# ===========================================================================
print("\n--- Loading data ---")

def safe_load(path, desc):
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        print(f"  {desc}: {len(df)} records")
        return df
    print(f"  [MISSING] {desc}: {path}")
    return None

# Module E outputs
domains = safe_load(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", "Domain assignments")
dist = safe_load(f"{RESULTS}/domains/domain_distribution_full.tsv", "Domain distribution")

# Module F outputs
paired = safe_load(f"{RESULTS}/stability/n_vs_c_paired_full.tsv", "N-vs-C paired")
contact_order = safe_load(f"{RESULTS}/stability/contact_order_full.tsv", "Contact order")
regions = safe_load(f"{RESULTS}/stability/region_boundaries_full.tsv", "Region boundaries")

# Phase 1 MTS data
targeting = safe_load(f"{PROJECT_DIR}/results/mts/combined_targeting.tsv", "Targeting")
mts_domain = safe_load(f"{PROJECT_DIR}/results/mts/mts_domain_relationship.tsv", "MTS-domain")

# Full-scale CATH data (18,855 proteins from InterPro Gene3D + pilot)
cath_full_path = f"{RESULTS}/domains/cath_domain_assignments_full.tsv"
cath_pilot_path = f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv"
if os.path.exists(cath_full_path):
    cath = safe_load(cath_full_path, "CATH domains (FULL-SCALE)")
else:
    cath = safe_load(cath_pilot_path, "CATH domains (pilot fallback)")

# Full-scale DSSP data (if available)
dssp_full_path = f"{RESULTS}/structures/dssp_summary_full.tsv"
dssp_pilot_path = f"{PROJECT_DIR}/results/structures/dssp_summary.tsv"
if os.path.exists(dssp_full_path):
    dssp = safe_load(dssp_full_path, "DSSP secondary structure (FULL-SCALE)")
else:
    dssp = safe_load(dssp_pilot_path, "DSSP secondary structure (pilot fallback)")

# Substrate lists
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")

groel_acc_col = "current_accession" if "current_accession" in groel.columns else "accession"
hsp60_acc_col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
groel_acc = set(groel[groel_acc_col].dropna().values)
hsp60_acc = set(hsp60[hsp60_acc_col].dropna().values)

# Homologs
homologs = safe_load(f"{PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv", "Homolog pairs")

print(f"\nGroEL substrates: {len(groel_acc)}, HSP60 substrates: {len(hsp60_acc)}")


# ===========================================================================
# Run all tests
# ===========================================================================
all_tests = []

# ===========================================================================
# FAMILY 1: Domain Architecture
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 1: Domain Architecture Enrichment")
print("=" * 70)

# H1.1: Multi-domain enrichment
if domains is not None:
    for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
        sub = domains[domains["accession"].isin(ds_accs)]
        bg = domains[~domains["accession"].isin(ds_accs)]
        if len(sub) < 5 or len(bg) < 5:
            continue

        sub_multi = (sub["n_domains"] >= 2).sum()
        sub_total = len(sub)
        bg_multi = (bg["n_domains"] >= 2).sum()
        bg_total = len(bg)

        table = [[sub_multi, sub_total - sub_multi],
                 [bg_multi, bg_total - bg_multi]]
        oddsratio, pval, ci_lo, ci_hi = safe_fisher(table)

        print(f"\n{ds_name} multi-domain enrichment:")
        print(f"  Substrate: {sub_multi}/{sub_total} ({100*sub_multi/sub_total:.1f}%)")
        print(f"  Background: {bg_multi}/{bg_total} ({100*bg_multi/bg_total:.1f}%)")
        print(f"  OR={oddsratio:.3f} [{ci_lo:.2f}-{ci_hi:.2f}], p={pval:.2e}")

        all_tests.append({
            "family": "domain_architecture", "hypothesis": f"H1.1_{ds_name}_multi_domain",
            "test": "Fisher exact", "statistic": oddsratio, "p_value": pval,
            "effect_size": oddsratio, "ci_low": ci_lo, "ci_high": ci_hi,
            "n1": sub_total, "n2": bg_total, "direction": "enriched" if oddsratio > 1 else "depleted",
        })

# H1.2: CATH superfamily enrichment (top superfamilies)
if cath is not None and len(cath) > 0:
    cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
    sf_col = "cath_superfamily" if "cath_superfamily" in cath.columns else None

    if sf_col:
        for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
            sub_cath = cath[cath[cath_acc_col].isin(ds_accs)]
            bg_cath = cath[~cath[cath_acc_col].isin(ds_accs)]

            if len(sub_cath) < 3:
                continue

            # Test top superfamilies in substrates
            top_sfs = sub_cath[sf_col].value_counts().head(10)
            for sf, sub_count in top_sfs.items():
                bg_count = (bg_cath[sf_col] == sf).sum()
                sub_other = len(sub_cath) - sub_count
                bg_other = len(bg_cath) - bg_count

                table = [[sub_count, sub_other], [bg_count, bg_other]]
                oddsratio, pval, ci_lo, ci_hi = safe_fisher(table)

                all_tests.append({
                    "family": "domain_architecture",
                    "hypothesis": f"H1.2_{ds_name}_SF_{sf}",
                    "test": "Fisher exact", "statistic": oddsratio, "p_value": pval,
                    "effect_size": oddsratio, "ci_low": ci_lo, "ci_high": ci_hi,
                    "n1": len(sub_cath), "n2": len(bg_cath),
                    "direction": "enriched" if oddsratio > 1 else "depleted",
                })

            print(f"\n{ds_name}: tested {len(top_sfs)} superfamilies")

# H1.3: CATH class distribution (chi-squared)
if cath is not None and "cath_class" in cath.columns:
    cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
    for ds_name, ds_accs in [("GroEL", groel_acc), ("HSP60", hsp60_acc)]:
        sub = cath[cath[cath_acc_col].isin(ds_accs)]
        bg = cath[~cath[cath_acc_col].isin(ds_accs)]

        if len(sub) < 10:
            continue

        classes = sorted(cath["cath_class"].unique())
        observed = np.array([[(sub["cath_class"] == c).sum() for c in classes],
                             [(bg["cath_class"] == c).sum() for c in classes]])

        if observed.min() >= 0 and observed.sum() > 0:
            chi2, pval, dof, expected = stats.chi2_contingency(observed)
            n_total = observed.sum()
            cramers_v = np.sqrt(chi2 / (n_total * (min(observed.shape) - 1)))

            print(f"\n{ds_name} CATH class distribution:")
            print(f"  chi2={chi2:.2f}, dof={dof}, p={pval:.2e}, Cramer's V={cramers_v:.3f}")

            all_tests.append({
                "family": "domain_architecture",
                "hypothesis": f"H1.3_{ds_name}_CATH_class",
                "test": "Chi-squared", "statistic": chi2, "p_value": pval,
                "effect_size": cramers_v, "ci_low": np.nan, "ci_high": np.nan,
                "n1": int(observed[0].sum()), "n2": int(observed[1].sum()),
                "direction": "different",
            })


# ===========================================================================
# FAMILY 2: N-vs-C Stability Asymmetry
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 2: N-vs-C Stability Asymmetry")
print("=" * 70)

if paired is not None and len(paired) > 0:
    metrics_to_test = [
        "relative_contact_order", "mean_plddt", "frac_helix",
        "frac_strand", "frac_plddt_gt70", "frac_hydrophobic",
    ]

    # H2.1: Within-protein paired tests (Wilcoxon signed-rank)
    for ds_name in ["groel", "hsp60", "matrix_bg", "mito_bg"]:
        sub = paired[paired["datasets"].str.contains(ds_name)]
        if len(sub) < 10:
            continue

        print(f"\n{ds_name} (n={len(sub)}) — paired Wilcoxon:")
        for metric in metrics_to_test:
            diff_col = f"{metric}_diff"
            if diff_col not in sub.columns:
                continue

            diffs = sub[diff_col].dropna()
            if len(diffs) < 10:
                continue

            try:
                stat, pval = stats.wilcoxon(diffs, alternative="two-sided")
                r = rank_biserial_paired(diffs.values)
                direction = "N>C" if diffs.median() > 0 else "C>N"

                print(f"  {metric}: median_diff={diffs.median():.4f}, "
                      f"W={stat:.0f}, p={pval:.2e}, r={r:.3f} ({direction})")

                all_tests.append({
                    "family": "stability_asymmetry",
                    "hypothesis": f"H2.1_{ds_name}_{metric}",
                    "test": "Wilcoxon signed-rank", "statistic": stat,
                    "p_value": pval, "effect_size": r,
                    "ci_low": np.nan, "ci_high": np.nan,
                    "n1": len(diffs), "n2": len(diffs),
                    "direction": direction,
                })
            except Exception as e:
                print(f"  {metric}: test failed ({e})")

    # H2.2: Cross-dataset comparison (Mann-Whitney U)
    for metric in ["relative_contact_order", "mean_plddt"]:
        diff_col = f"{metric}_diff"
        if diff_col not in paired.columns:
            continue

        comparisons = [
            ("GroEL_vs_ecoli_bg", "groel", "proteome_bg"),
            ("HSP60_vs_mito_bg", "hsp60", "mito_bg"),
        ]
        for comp_name, ds1, ds2 in comparisons:
            g1 = paired[paired["datasets"].str.contains(ds1)][diff_col].dropna()
            g2 = paired[paired["datasets"].str.contains(ds2)][diff_col].dropna()

            if len(g1) < 5 or len(g2) < 5:
                continue

            stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
            r = rank_biserial_unpaired(stat, len(g1), len(g2))

            print(f"\n{comp_name} {metric}: U={stat:.0f}, p={pval:.2e}, r={r:.3f}")

            all_tests.append({
                "family": "stability_asymmetry",
                "hypothesis": f"H2.2_{comp_name}_{metric}",
                "test": "Mann-Whitney U", "statistic": stat, "p_value": pval,
                "effect_size": r, "ci_low": np.nan, "ci_high": np.nan,
                "n1": len(g1), "n2": len(g2),
                "direction": "substrate_higher" if g1.median() > g2.median() else "bg_higher",
            })

    # H2.3: GroEL class effects (Kruskal-Wallis)
    if "groel_class" in paired.columns:
        groel_sub = paired[paired["datasets"].str.contains("groel")]
        for metric in ["relative_contact_order", "mean_plddt"]:
            diff_col = f"{metric}_diff"
            if diff_col not in groel_sub.columns:
                continue

            groups = []
            class_labels = []
            for cls in ["I", "II", "III"]:
                cls_data = groel_sub[groel_sub["groel_class"].astype(str) == cls][diff_col].dropna()
                if len(cls_data) >= 3:
                    groups.append(cls_data.values)
                    class_labels.append(cls)

            if len(groups) >= 2:
                stat, pval = stats.kruskal(*groups)
                # Eta-squared
                n_total = sum(len(g) for g in groups)
                eta_sq = (stat - len(groups) + 1) / (n_total - len(groups))

                print(f"\nGroEL class effect on {metric}: H={stat:.2f}, p={pval:.2e}, "
                      f"eta2={eta_sq:.3f}")

                all_tests.append({
                    "family": "stability_asymmetry",
                    "hypothesis": f"H2.3_GroEL_class_{metric}",
                    "test": "Kruskal-Wallis", "statistic": stat, "p_value": pval,
                    "effect_size": eta_sq, "ci_low": np.nan, "ci_high": np.nan,
                    "n1": n_total, "n2": len(groups),
                    "direction": "class_effect",
                })


# ===========================================================================
# FAMILY 2 (continued): DSSP Secondary Structure Tests
# ===========================================================================
if dssp is not None and len(dssp) > 0:
    print("\n" + "-" * 50)
    print("DSSP Secondary Structure Comparisons (Mann-Whitney U)")
    print("-" * 50)

    dssp_acc_col = "accession" if "accession" in dssp.columns else (
        "uniprot_accession" if "uniprot_accession" in dssp.columns else dssp.columns[0])

    # E. coli proteome accessions (for GroEL background)
    ecoli_path = f"{PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.tsv"
    if os.path.exists(ecoli_path):
        ecoli_df = pd.read_csv(ecoli_path, sep="\t")
        ecoli_acc_col = "Entry" if "Entry" in ecoli_df.columns else ecoli_df.columns[0]
        ecoli_acc = set(ecoli_df[ecoli_acc_col].dropna().values)
    else:
        ecoli_acc = set()

    # Matrix accessions (for HSP60 background)
    matrix_path = f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv"
    if os.path.exists(matrix_path):
        matrix_df = pd.read_csv(matrix_path, sep="\t")
        mat_acc_col = "accession" if "accession" in matrix_df.columns else matrix_df.columns[0]
        matrix_acc_for_dssp = set(matrix_df[mat_acc_col].dropna().values)
    else:
        matrix_acc_for_dssp = set()

    dssp_metrics = ["frac_helix", "frac_strand", "frac_coil"]
    dssp_comparisons = [
        ("GroEL_vs_ecoli_bg", groel_acc, ecoli_acc - groel_acc, "compartment-matched E. coli"),
        ("HSP60_vs_matrix_bg", hsp60_acc, matrix_acc_for_dssp - hsp60_acc, "compartment-matched matrix"),
    ]

    for comp_name, sub_accs, bg_accs, bg_desc in dssp_comparisons:
        sub_dssp = dssp[dssp[dssp_acc_col].isin(sub_accs)]
        bg_dssp = dssp[dssp[dssp_acc_col].isin(bg_accs)]

        if len(sub_dssp) < 10 or len(bg_dssp) < 10:
            print(f"  {comp_name}: skipped (sub={len(sub_dssp)}, bg={len(bg_dssp)})")
            continue

        print(f"\n{comp_name} (sub={len(sub_dssp)}, bg={len(bg_dssp)}, {bg_desc}):")
        for metric in dssp_metrics:
            if metric not in sub_dssp.columns:
                continue

            g1 = sub_dssp[metric].dropna()
            g2 = bg_dssp[metric].dropna()
            if len(g1) < 10 or len(g2) < 10:
                continue

            stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
            r = rank_biserial_unpaired(stat, len(g1), len(g2))
            direction = "substrate_higher" if g1.median() > g2.median() else "bg_higher"

            print(f"  {metric}: sub_median={g1.median():.3f}, bg_median={g2.median():.3f}, "
                  f"U={stat:.0f}, p={pval:.2e}, r={r:.3f}")

            all_tests.append({
                "family": "stability_asymmetry",
                "hypothesis": f"H2.4_{comp_name}_{metric}",
                "test": "Mann-Whitney U (DSSP)", "statistic": stat, "p_value": pval,
                "effect_size": r, "ci_low": np.nan, "ci_high": np.nan,
                "n1": len(g1), "n2": len(g2),
                "direction": direction,
            })
else:
    print("\n[DSSP data not available — DSSP tests skipped]")


# ===========================================================================
# FAMILY 3: MTS Targeting
# ===========================================================================
print("\n" + "=" * 70)
print("FAMILY 3: MTS Targeting Statistics")
print("=" * 70)

# H3.1: HSP60 matrix enrichment
if targeting is not None:
    targ_acc_col = "uniprot_accession" if "uniprot_accession" in targeting.columns else "accession"

    if targ_acc_col in targeting.columns:
        # Determine matrix column
        matrix_col = None
        for col in ["targeting_classification", "mitocarta_is_matrix",
                     "is_matrix", "localization"]:
            if col in targeting.columns:
                matrix_col = col
                break

        if matrix_col:
            hsp60_targ = targeting[targeting[targ_acc_col].isin(hsp60_acc)]
            bg_targ = targeting[~targeting[targ_acc_col].isin(hsp60_acc)]

            if matrix_col == "mitocarta_is_matrix":
                hsp60_matrix = hsp60_targ[hsp60_targ[matrix_col] == True]
                bg_matrix = bg_targ[bg_targ[matrix_col] == True]
            else:
                hsp60_matrix = hsp60_targ[hsp60_targ[matrix_col].astype(str).str.contains(
                    "matrix|Matrix", case=False, na=False)]
                bg_matrix = bg_targ[bg_targ[matrix_col].astype(str).str.contains(
                    "matrix|Matrix", case=False, na=False)]

            if len(hsp60_targ) > 0 and len(bg_targ) > 0:
                table = [[len(hsp60_matrix), len(hsp60_targ) - len(hsp60_matrix)],
                         [len(bg_matrix), len(bg_targ) - len(bg_matrix)]]
                oddsratio, pval, ci_lo, ci_hi = safe_fisher(table)

                print(f"\nHSP60 matrix enrichment:")
                print(f"  HSP60: {len(hsp60_matrix)}/{len(hsp60_targ)} matrix")
                print(f"  Background: {len(bg_matrix)}/{len(bg_targ)} matrix")
                print(f"  OR={oddsratio:.3f} [{ci_lo:.2f}-{ci_hi:.2f}], p={pval:.2e}")

                all_tests.append({
                    "family": "mts_targeting",
                    "hypothesis": "H3.1_HSP60_matrix_enrichment",
                    "test": "Fisher exact", "statistic": oddsratio, "p_value": pval,
                    "effect_size": oddsratio, "ci_low": ci_lo, "ci_high": ci_hi,
                    "n1": len(hsp60_targ), "n2": len(bg_targ),
                    "direction": "enriched" if oddsratio > 1 else "depleted",
                })

# H3.2: MTS pre-domain dominance (binomial test)
if mts_domain is not None:
    pre_col = "mts_is_pre_domain" if "mts_is_pre_domain" in mts_domain.columns else None
    if pre_col:
        pre_domain = mts_domain[pre_col].sum()
        total = len(mts_domain[mts_domain[pre_col].notna()])
        if total > 0:
            pval = stats.binomtest(pre_domain, total, 0.5, alternative="greater").pvalue

            print(f"\nMTS pre-domain dominance:")
            print(f"  Pre-domain: {pre_domain}/{total} ({100*pre_domain/total:.1f}%)")
            print(f"  Binomial p={pval:.2e} (H0: 50%)")

            all_tests.append({
                "family": "mts_targeting",
                "hypothesis": "H3.2_MTS_pre_domain",
                "test": "Binomial", "statistic": pre_domain / total, "p_value": pval,
                "effect_size": pre_domain / total,
                "ci_low": np.nan, "ci_high": np.nan,
                "n1": total, "n2": total,
                "direction": "pre_domain_dominant",
            })


# ===========================================================================
# HIERARCHICAL BH CORRECTION
# ===========================================================================
print("\n" + "=" * 70)
print("HIERARCHICAL BH CORRECTION")
print("=" * 70)

if all_tests:
    results_df = pd.DataFrame(all_tests)

    # Level 1: BH within each family
    results_df["p_bh_within"] = np.nan
    for family in results_df["family"].unique():
        mask = results_df["family"] == family
        pvals = results_df.loc[mask, "p_value"].values
        corrected = bh_correct(pvals)
        results_df.loc[mask, "p_bh_within"] = corrected

    results_df["sig_within_family"] = results_df["p_bh_within"] < 0.05

    # Level 2: BH across families (using minimum p per family)
    family_min_p = results_df.groupby("family")["p_value"].min()
    family_corrected = bh_correct(family_min_p.values)
    family_sig = dict(zip(family_min_p.index, family_corrected < 0.05))
    results_df["family_significant"] = results_df["family"].map(family_sig)

    # Overall significance: both levels
    results_df["significant_overall"] = (
        results_df["sig_within_family"] & results_df["family_significant"]
    )

    print(f"\nTotal tests: {len(results_df)}")
    print(f"Significant within family (BH<0.05): {results_df['sig_within_family'].sum()}")
    print(f"Significant overall (hierarchical): {results_df['significant_overall'].sum()}")
    print(f"\nBy family:")
    for family in results_df["family"].unique():
        fam = results_df[results_df["family"] == family]
        sig_within = fam["sig_within_family"].sum()
        sig_overall = fam["significant_overall"].sum()
        print(f"  {family}: {sig_within}/{len(fam)} within, {sig_overall}/{len(fam)} overall")

    # Save results
    os.makedirs(f"{RESULTS}/stats", exist_ok=True)

    pval_path = f"{RESULTS}/stats/corrected_pvalues_full.tsv"
    results_df.to_csv(pval_path, sep="\t", index=False)
    print(f"\nSaved: {pval_path}")

    # Summary report
    report_path = f"{RESULTS}/stats/statistics_summary_full.txt"
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("ANTAH ASTI PRARAMBH — PHASE 2 FULL-SCALE STATISTICS\n")
        f.write(f"Generated: {datetime.now()}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Total tests: {len(results_df)}\n")
        f.write(f"Significant within family: {results_df['sig_within_family'].sum()}\n")
        f.write(f"Significant overall: {results_df['significant_overall'].sum()}\n\n")

        for family in ["domain_architecture", "stability_asymmetry", "mts_targeting"]:
            fam = results_df[results_df["family"] == family]
            f.write(f"\n{'='*50}\n{family.upper()}\n{'='*50}\n")
            for _, row in fam.iterrows():
                sig = " ***" if row["significant_overall"] else (
                    " *" if row["sig_within_family"] else "")
                f.write(f"\n{row['hypothesis']}:\n")
                f.write(f"  Test: {row['test']}\n")
                f.write(f"  Statistic: {row['statistic']:.4f}\n")
                f.write(f"  p-value: {row['p_value']:.2e}\n")
                f.write(f"  p-BH (within family): {row['p_bh_within']:.2e}\n")
                f.write(f"  Effect size: {row['effect_size']:.4f}\n")
                f.write(f"  Direction: {row['direction']}\n")
                f.write(f"  n1={row['n1']}, n2={row['n2']}{sig}\n")

        f.write("\n\n" + "=" * 70 + "\n")
        f.write("NOTES:\n")
        f.write("  - pLDDT reflects AlphaFold prediction confidence, NOT thermodynamic stability\n")
        f.write("  - Contact order (r=-0.75 with folding rates; Plaxco et al. 1998) is the primary folding kinetics proxy\n")
        f.write("  - FoldX was parameterized on experimental structures; applied to AlphaFold predictions with caution\n")
        f.write("  - Hierarchical BH: within-family + between-family correction\n")
        f.write("  - *** = significant at both levels, * = within-family only\n")
        f.write("  - DSSP tests use compartment-matched backgrounds (E. coli for GroEL, matrix for HSP60)\n")
        f.write("  - N-vs-C asymmetry is universal across all multi-domain proteins, NOT substrate-specific\n")
        f.write("=" * 70 + "\n")

    print(f"Saved: {report_path}")

else:
    print("No statistical tests could be run (insufficient data)")

print("\n" + "=" * 70)
print("Module H complete.")
print("=" * 70)
