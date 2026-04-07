#!/usr/bin/env python3
"""
Sensitivity Analysis for Key Findings

Tests robustness of the three most critical parameter choices:
  1. SILAC threshold for HSP60 substrate selection
  2. Size-matching bin width for CATH enrichment tests
  3. Background multiplier for enrichment tests

Outputs:
  - results/phase2/stats/sensitivity_analysis.tsv
  - results/phase2/figures/fig8_sensitivity.pdf / .png

Author: Vishal Bharti (with Claude Code assistance)
Date: 2026-04-07
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# Paths
# ============================================================
PROJECT_DIR = os.path.expanduser("~/Downloads/Antah_Asti_Prarambh")

# Input files
HSP60_INTERACTOME = os.path.join(PROJECT_DIR, "data/processed/hsp60_interactome_standardized.tsv")
HSP60_TIER1 = os.path.join(PROJECT_DIR, "data/processed/hsp60_tier1_substrates.tsv")
GROEL_SUBS = os.path.join(PROJECT_DIR, "data/processed/groel_substrates_standardized.tsv")
CATH_DOMAINS = os.path.join(PROJECT_DIR, "results/phase2/domains/cath_domain_assignments_full.tsv")
NVC_PAIRED = os.path.join(PROJECT_DIR, "results/phase2/stability/n_vs_c_paired_full.tsv")
ECOLI_PROTEOME = os.path.join(PROJECT_DIR, "data/raw/uniprot/ecoli_k12_proteome.tsv")
MATRIX_PROTEOME = os.path.join(PROJECT_DIR, "data/processed/human_matrix_proteome.tsv")

# Output paths
STATS_DIR = os.path.join(PROJECT_DIR, "results/phase2/stats")
FIG_DIR = os.path.join(PROJECT_DIR, "results/phase2/figures")
os.makedirs(STATS_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

OUT_TSV = os.path.join(STATS_DIR, "sensitivity_analysis.tsv")
OUT_PDF = os.path.join(FIG_DIR, "fig8_sensitivity.pdf")
OUT_PNG = os.path.join(FIG_DIR, "fig8_sensitivity.png")

# ============================================================
# Colorblind-friendly palette (Okabe-Ito)
# ============================================================
CB_BLUE = "#0072B2"
CB_ORANGE = "#E69F00"
CB_GREEN = "#009E73"
CB_RED = "#D55E00"
CB_PURPLE = "#CC79A7"
CB_CYAN = "#56B4E9"
CB_YELLOW = "#F0E442"
CB_BLACK = "#000000"


# ============================================================
# Helper functions
# ============================================================

def load_tsv(path, description="file"):
    """Load a TSV file with error handling."""
    try:
        df = pd.read_csv(path, sep="\t")
        print(f"  Loaded {description}: {len(df)} rows, {len(df.columns)} cols")
        return df
    except FileNotFoundError:
        print(f"  ERROR: {description} not found at {path}")
        return None
    except Exception as e:
        print(f"  ERROR loading {description}: {e}")
        return None


def fisher_exact_enrichment(n_target_in_group, n_group, n_target_in_bg, n_bg):
    """
    Fisher exact test for enrichment.
    Returns (odds_ratio, p_value).

    Contingency table:
                   target    not-target
    group         a          b
    background    c          d
    """
    a = n_target_in_group
    b = n_group - a
    c = n_target_in_bg
    d = n_bg - c

    # Ensure non-negative
    a, b, c, d = max(0, a), max(0, b), max(0, c), max(0, d)

    try:
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]],
                                                  alternative="greater")
    except Exception:
        odds_ratio, p_value = np.nan, np.nan

    return odds_ratio, p_value


def size_matched_background(substrate_masses, bg_masses, bg_accessions, bg_domains_df,
                            bin_width_kda, multiplier, rng):
    """
    Draw a size-matched background from bg_accessions.

    For each substrate, find bg proteins in the same mass bin and sample
    `multiplier` of them (with replacement if not enough).

    Returns a set of sampled background accessions.
    """
    sampled = []
    for mass in substrate_masses:
        bin_low = (mass // bin_width_kda) * bin_width_kda
        bin_high = bin_low + bin_width_kda
        candidates = bg_accessions[
            (bg_masses >= bin_low) & (bg_masses < bin_high)
        ]
        if len(candidates) == 0:
            # Widen to nearest bin
            continue
        n_to_sample = min(multiplier, len(candidates))
        chosen = rng.choice(candidates.values, size=n_to_sample, replace=False)
        sampled.extend(chosen)

    return set(sampled)


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 72)
    print("SENSITIVITY ANALYSIS")
    print("=" * 72)
    print()

    results = []  # Will become the output TSV

    # ----------------------------------------------------------
    # Load all data
    # ----------------------------------------------------------
    print("--- Loading data ---")

    hsp60_all = load_tsv(HSP60_INTERACTOME, "HSP60 interactome")
    groel = load_tsv(GROEL_SUBS, "GroEL substrates")
    cath = load_tsv(CATH_DOMAINS, "CATH domain assignments")
    nvc = load_tsv(NVC_PAIRED, "N-vs-C paired data")
    ecoli_prot = load_tsv(ECOLI_PROTEOME, "E. coli proteome")
    matrix_prot = load_tsv(MATRIX_PROTEOME, "Matrix proteome")
    print()

    # Validate essential files
    missing = []
    if hsp60_all is None:
        missing.append("HSP60 interactome")
    if groel is None:
        missing.append("GroEL substrates")
    if cath is None:
        missing.append("CATH domains")
    if nvc is None:
        missing.append("N-vs-C paired")
    if ecoli_prot is None:
        missing.append("E. coli proteome")
    if matrix_prot is None:
        missing.append("Matrix proteome")

    if missing:
        print(f"FATAL: Missing essential files: {', '.join(missing)}")
        print("Cannot proceed with sensitivity analysis.")
        sys.exit(1)

    # ----------------------------------------------------------
    # Prepare derived datasets
    # ----------------------------------------------------------

    # GroEL accessions
    groel_accessions = set(groel["current_accession"].dropna().unique())

    # Matrix accessions (background for HSP60)
    matrix_accessions = set(matrix_prot["uniprot_id"].dropna().unique())

    # E. coli proteome: compute approximate mass in kDa from Length
    # (average amino acid mass ~110 Da = 0.110 kDa)
    ecoli_prot = ecoli_prot.copy()
    ecoli_prot["mass_kDa"] = pd.to_numeric(ecoli_prot["Length"], errors="coerce") * 0.110
    ecoli_accessions = set(ecoli_prot["Entry"].dropna().unique())

    # E. coli non-GroEL background
    ecoli_bg = ecoli_prot[~ecoli_prot["Entry"].isin(groel_accessions)].copy()

    # GroEL with mass
    groel_with_mass = groel[["current_accession", "mass_kDa"]].dropna().copy()
    groel_with_mass = groel_with_mass.rename(columns={"current_accession": "accession"})

    # CATH: count domains per protein and identify TIM barrels (3.20.20.70)
    cath_domain_counts = cath.groupby("accession").size().reset_index(name="n_domains")
    cath_superfamilies = cath.groupby("accession")["cath_superfamily"].apply(set).reset_index()
    cath_superfamilies.columns = ["accession", "superfamily_set"]

    # TIM barrel check
    TIM_BARREL = "3.20.20.70"

    def has_tim_barrel(sf_set):
        return TIM_BARREL in sf_set

    cath_superfamilies["has_tim_barrel"] = cath_superfamilies["superfamily_set"].apply(
        has_tim_barrel
    )

    # Multi-domain check
    cath_domain_counts["is_multi_domain"] = cath_domain_counts["n_domains"] >= 2

    # ----------------------------------------------------------
    # ANALYSIS 1: SILAC Threshold Sensitivity
    # ----------------------------------------------------------
    print("=" * 72)
    print("ANALYSIS 1: SILAC Threshold Sensitivity for HSP60 Substrate Selection")
    print("=" * 72)
    print()

    # The original filter (filter_hsp60_interactome.py):
    #   Tier 1 = candidate_substrate + mitocarta=='X' + median_silac_imputed > 5
    #
    # We vary the threshold on median_silac_imputed (as a raw ratio, not log2):
    #   thresholds = [3, 5, 7, 10]
    #
    # The log2 equivalents are approximately [1.58, 2.32, 2.81, 3.32].

    silac_thresholds = [3, 5, 7, 10]

    # Only candidate substrates with mitocarta annotation count as tier1 pool
    hsp60_candidates = hsp60_all[
        (hsp60_all["flag"] == "candidate_substrate") &
        (hsp60_all["mitocarta"].astype(str).str.strip().str.upper() == "X")
    ].copy()

    print(f"  HSP60 candidate substrates with MitoCarta: {len(hsp60_candidates)}")
    print()

    for thresh in silac_thresholds:
        log2_thresh = np.log2(thresh)
        hsp60_at_thresh = hsp60_candidates[
            hsp60_candidates["median_silac_imputed"] > thresh
        ]
        hsp60_ids = set(hsp60_at_thresh["uniprot_id"].dropna().unique())
        n_hsp60 = len(hsp60_ids)

        print(f"  SILAC ratio > {thresh} (log2 > {log2_thresh:.2f}): {n_hsp60} HSP60 substrates")

        # --- Test 1a: Multi-domain enrichment (HSP60 vs matrix background) ---
        # HSP60 substrates: how many are multi-domain?
        hsp60_domains = cath_domain_counts[cath_domain_counts["accession"].isin(hsp60_ids)]
        n_hsp60_with_cath = len(hsp60_domains)
        n_hsp60_multi = int(hsp60_domains["is_multi_domain"].sum())

        # Matrix bg (excluding HSP60 substrates)
        matrix_bg_ids = matrix_accessions - hsp60_ids
        matrix_bg_domains = cath_domain_counts[cath_domain_counts["accession"].isin(matrix_bg_ids)]
        n_matrix_with_cath = len(matrix_bg_domains)
        n_matrix_multi = int(matrix_bg_domains["is_multi_domain"].sum())

        or_multi, p_multi = fisher_exact_enrichment(
            n_hsp60_multi, n_hsp60_with_cath,
            n_matrix_multi, n_matrix_with_cath
        )

        print(f"    Multi-domain enrichment: OR={or_multi:.2f}, p={p_multi:.2e} "
              f"(HSP60: {n_hsp60_multi}/{n_hsp60_with_cath}, "
              f"matrix bg: {n_matrix_multi}/{n_matrix_with_cath})")

        results.append({
            "analysis": "SILAC_threshold",
            "parameter": "median_silac_ratio",
            "value": thresh,
            "test": "HSP60_multi_domain_enrichment",
            "n_substrates": n_hsp60_with_cath,
            "n_background": n_matrix_with_cath,
            "statistic": or_multi,
            "p_value": p_multi,
            "effect_size": or_multi,
        })

        # --- Test 1b: N-vs-C contact order asymmetry for HSP60 ---
        # Filter nvc for HSP60 substrates at this threshold
        nvc_hsp60 = nvc[
            nvc["accession"].isin(hsp60_ids) &
            nvc["relative_contact_order_diff"].notna()
        ]
        n_nvc_hsp60 = len(nvc_hsp60)

        if n_nvc_hsp60 >= 5:
            co_diffs = nvc_hsp60["relative_contact_order_diff"].values
            try:
                stat_wilcox, p_wilcox = stats.wilcoxon(co_diffs, alternative="greater")
            except Exception:
                stat_wilcox, p_wilcox = np.nan, np.nan
            median_diff = np.median(co_diffs)
            print(f"    N>C contact order asymmetry: W={stat_wilcox:.1f}, p={p_wilcox:.2e}, "
                  f"median_diff={median_diff:.4f} (n={n_nvc_hsp60})")
        else:
            stat_wilcox, p_wilcox, median_diff = np.nan, np.nan, np.nan
            print(f"    N>C contact order asymmetry: too few observations (n={n_nvc_hsp60})")

        results.append({
            "analysis": "SILAC_threshold",
            "parameter": "median_silac_ratio",
            "value": thresh,
            "test": "HSP60_NvsC_contact_order_asymmetry",
            "n_substrates": n_nvc_hsp60,
            "n_background": 0,
            "statistic": stat_wilcox,
            "p_value": p_wilcox,
            "effect_size": median_diff,
        })

    print()

    # ----------------------------------------------------------
    # ANALYSIS 2: Size-Matching Bin Width Sensitivity
    # ----------------------------------------------------------
    print("=" * 72)
    print("ANALYSIS 2: Size-Matching Bin Width Sensitivity (GroEL TIM Barrel)")
    print("=" * 72)
    print()

    bin_widths = [5, 10, 15, 20]
    rng = np.random.default_rng(42)  # Reproducible

    # GroEL substrate masses
    groel_masses = groel_with_mass["mass_kDa"].values
    groel_accs = groel_with_mass["accession"].values

    # E. coli bg masses and accessions
    bg_masses_series = ecoli_bg["mass_kDa"].dropna()
    bg_accessions_series = ecoli_bg.loc[bg_masses_series.index, "Entry"]
    # Reset index so we can use .values cleanly
    bg_masses_arr = bg_masses_series.values
    bg_accessions_arr = bg_accessions_series.values

    # GroEL: which have TIM barrel?
    groel_in_cath = cath_superfamilies[cath_superfamilies["accession"].isin(groel_accessions)]
    n_groel_in_cath = len(groel_in_cath)
    n_groel_tim = int(groel_in_cath["has_tim_barrel"].sum())

    print(f"  GroEL substrates with CATH data: {n_groel_in_cath}")
    print(f"  GroEL substrates with TIM barrel: {n_groel_tim}")
    print()

    for bin_width in bin_widths:
        # Draw size-matched background
        sampled_bg = set()
        for mass in groel_masses:
            if np.isnan(mass):
                continue
            bin_low = (mass // bin_width) * bin_width
            bin_high = bin_low + bin_width
            mask = (bg_masses_arr >= bin_low) & (bg_masses_arr < bin_high)
            candidates = bg_accessions_arr[mask]
            if len(candidates) == 0:
                continue
            n_to_sample = min(3, len(candidates))  # Default multiplier=3 for this analysis
            chosen = rng.choice(candidates, size=n_to_sample, replace=False)
            sampled_bg.update(chosen)

        # Background: which have TIM barrel?
        bg_in_cath = cath_superfamilies[cath_superfamilies["accession"].isin(sampled_bg)]
        n_bg_in_cath = len(bg_in_cath)
        n_bg_tim = int(bg_in_cath["has_tim_barrel"].sum())

        or_tim, p_tim = fisher_exact_enrichment(
            n_groel_tim, n_groel_in_cath,
            n_bg_tim, n_bg_in_cath
        )

        print(f"  Bin width {bin_width} kDa: bg_sampled={len(sampled_bg)}, bg_with_CATH={n_bg_in_cath}")
        print(f"    TIM barrel: GroEL {n_groel_tim}/{n_groel_in_cath}, bg {n_bg_tim}/{n_bg_in_cath}")
        print(f"    OR={or_tim:.2f}, p={p_tim:.2e}")

        results.append({
            "analysis": "bin_width",
            "parameter": "bin_width_kDa",
            "value": bin_width,
            "test": "GroEL_TIM_barrel_enrichment",
            "n_substrates": n_groel_in_cath,
            "n_background": n_bg_in_cath,
            "statistic": or_tim,
            "p_value": p_tim,
            "effect_size": or_tim,
        })

    print()

    # ----------------------------------------------------------
    # ANALYSIS 3: Background Multiplier Sensitivity
    # ----------------------------------------------------------
    print("=" * 72)
    print("ANALYSIS 3: Background Multiplier Sensitivity (GroEL TIM Barrel)")
    print("=" * 72)
    print()

    multipliers = [1, 2, 3, 5]
    rng2 = np.random.default_rng(42)  # Fresh RNG for reproducibility

    for mult in multipliers:
        # Draw size-matched background with this multiplier (bin_width=10 kDa, the default)
        sampled_bg = set()
        for mass in groel_masses:
            if np.isnan(mass):
                continue
            bin_low = (mass // 10) * 10
            bin_high = bin_low + 10
            mask = (bg_masses_arr >= bin_low) & (bg_masses_arr < bin_high)
            candidates = bg_accessions_arr[mask]
            if len(candidates) == 0:
                continue
            n_to_sample = min(mult, len(candidates))
            chosen = rng2.choice(candidates, size=n_to_sample, replace=False)
            sampled_bg.update(chosen)

        # Background: which have TIM barrel?
        bg_in_cath = cath_superfamilies[cath_superfamilies["accession"].isin(sampled_bg)]
        n_bg_in_cath = len(bg_in_cath)
        n_bg_tim = int(bg_in_cath["has_tim_barrel"].sum())

        or_tim, p_tim = fisher_exact_enrichment(
            n_groel_tim, n_groel_in_cath,
            n_bg_tim, n_bg_in_cath
        )

        print(f"  Multiplier {mult}x: bg_sampled={len(sampled_bg)}, bg_with_CATH={n_bg_in_cath}")
        print(f"    TIM barrel: GroEL {n_groel_tim}/{n_groel_in_cath}, bg {n_bg_tim}/{n_bg_in_cath}")
        print(f"    OR={or_tim:.2f}, p={p_tim:.2e}")

        results.append({
            "analysis": "bg_multiplier",
            "parameter": "background_multiplier",
            "value": mult,
            "test": "GroEL_TIM_barrel_enrichment",
            "n_substrates": n_groel_in_cath,
            "n_background": n_bg_in_cath,
            "statistic": or_tim,
            "p_value": p_tim,
            "effect_size": or_tim,
        })

    print()

    # ----------------------------------------------------------
    # Save results TSV
    # ----------------------------------------------------------
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved results to: {OUT_TSV}")
    print()

    # ----------------------------------------------------------
    # Generate Figure 8: Sensitivity Analysis
    # ----------------------------------------------------------
    print("--- Generating Figure 8 ---")

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.40, wspace=0.35,
                           left=0.08, right=0.95, top=0.93, bottom=0.08)

    # ---- Panel A: SILAC threshold -- p-values ----
    ax_a = fig.add_subplot(gs[0, 0])

    silac_res = results_df[results_df["analysis"] == "SILAC_threshold"]
    silac_multi = silac_res[silac_res["test"] == "HSP60_multi_domain_enrichment"]
    silac_nvc = silac_res[silac_res["test"] == "HSP60_NvsC_contact_order_asymmetry"]

    ax_a.plot(silac_multi["value"], -np.log10(silac_multi["p_value"].clip(lower=1e-300)),
              "o-", color=CB_BLUE, linewidth=2, markersize=8,
              label="Multi-domain enrichment")
    ax_a.plot(silac_nvc["value"], -np.log10(silac_nvc["p_value"].clip(lower=1e-300)),
              "s-", color=CB_ORANGE, linewidth=2, markersize=8,
              label="N>C contact order")

    ax_a.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.6, linewidth=1)
    ax_a.text(silac_thresholds[-1] + 0.2, -np.log10(0.05) + 0.1,
              "p=0.05", fontsize=8, color="gray", va="bottom")

    ax_a.set_xlabel("SILAC ratio threshold", fontsize=11)
    ax_a.set_ylabel("$-\\log_{10}$(p-value)", fontsize=11)
    ax_a.set_title("A. SILAC threshold sensitivity", fontsize=12, fontweight="bold")
    ax_a.set_xticks(silac_thresholds)
    ax_a.legend(fontsize=9, loc="best", framealpha=0.8)

    # ---- Panel B: SILAC threshold -- sample size and OR ----
    ax_b = fig.add_subplot(gs[0, 1])

    ax_b.bar(np.array(silac_thresholds) - 0.3,
             silac_multi["n_substrates"].values,
             width=0.6, color=CB_CYAN, alpha=0.8, label="n (CATH-assigned HSP60)")

    ax_b2 = ax_b.twinx()
    ax_b2.plot(silac_multi["value"], silac_multi["effect_size"],
               "D-", color=CB_RED, linewidth=2, markersize=7,
               label="Odds ratio")

    ax_b.set_xlabel("SILAC ratio threshold", fontsize=11)
    ax_b.set_ylabel("Sample size (n)", fontsize=11, color=CB_BLUE)
    ax_b2.set_ylabel("Odds ratio (multi-domain)", fontsize=11, color=CB_RED)
    ax_b.set_title("B. Sample size vs effect size", fontsize=12, fontweight="bold")
    ax_b.set_xticks(silac_thresholds)

    # Combined legend
    lines_b1, labels_b1 = ax_b.get_legend_handles_labels()
    lines_b2, labels_b2 = ax_b2.get_legend_handles_labels()
    ax_b.legend(lines_b1 + lines_b2, labels_b1 + labels_b2,
                fontsize=9, loc="best", framealpha=0.8)

    # ---- Panel C: Bin width sensitivity ----
    ax_c = fig.add_subplot(gs[1, 0])

    bw_res = results_df[results_df["analysis"] == "bin_width"]

    color_bw = CB_GREEN
    ax_c.plot(bw_res["value"], -np.log10(bw_res["p_value"].clip(lower=1e-300)),
              "o-", color=color_bw, linewidth=2, markersize=8)
    ax_c.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.6, linewidth=1)

    # Annotate OR values
    for _, row in bw_res.iterrows():
        ax_c.annotate(f"OR={row['effect_size']:.1f}",
                      (row["value"], -np.log10(max(row["p_value"], 1e-300))),
                      textcoords="offset points", xytext=(0, 10),
                      fontsize=8, ha="center", color=color_bw)

    ax_c.set_xlabel("Size-matching bin width (kDa)", fontsize=11)
    ax_c.set_ylabel("$-\\log_{10}$(p-value)", fontsize=11)
    ax_c.set_title("C. Bin width sensitivity (GroEL TIM barrel)", fontsize=12, fontweight="bold")
    ax_c.set_xticks(bin_widths)

    # ---- Panel D: Background multiplier sensitivity ----
    ax_d = fig.add_subplot(gs[1, 1])

    mult_res = results_df[results_df["analysis"] == "bg_multiplier"]

    color_mult = CB_PURPLE
    ax_d.plot(mult_res["value"], -np.log10(mult_res["p_value"].clip(lower=1e-300)),
              "o-", color=color_mult, linewidth=2, markersize=8)
    ax_d.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.6, linewidth=1)

    # Annotate OR values
    for _, row in mult_res.iterrows():
        ax_d.annotate(f"OR={row['effect_size']:.1f}",
                      (row["value"], -np.log10(max(row["p_value"], 1e-300))),
                      textcoords="offset points", xytext=(0, 10),
                      fontsize=8, ha="center", color=color_mult)

    ax_d.set_xlabel("Background multiplier", fontsize=11)
    ax_d.set_ylabel("$-\\log_{10}$(p-value)", fontsize=11)
    ax_d.set_title("D. Multiplier sensitivity (GroEL TIM barrel)", fontsize=12, fontweight="bold")
    ax_d.set_xticks(multipliers)

    # Suptitle
    fig.suptitle("Figure 8: Sensitivity Analysis of Key Parameter Choices",
                 fontsize=14, fontweight="bold", y=0.98)

    # Save
    try:
        fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
        print(f"  Saved: {OUT_PDF}")
    except Exception as e:
        print(f"  WARNING: Could not save PDF: {e}")

    try:
        fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
        print(f"  Saved: {OUT_PNG}")
    except Exception as e:
        print(f"  WARNING: Could not save PNG: {e}")

    plt.close(fig)
    print()

    # ----------------------------------------------------------
    # Print Summary
    # ----------------------------------------------------------
    print("=" * 72)
    print("SENSITIVITY ANALYSIS SUMMARY")
    print("=" * 72)
    print()

    # Analysis 1: SILAC
    print("1. SILAC Threshold Sensitivity")
    print("   ----------------------------")
    for _, row in silac_res.iterrows():
        sig = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01 else
               ("*" if row["p_value"] < 0.05 else "NS"))
        print(f"   Threshold>{row['value']:.0f} | {row['test']:<40s} | "
              f"p={row['p_value']:.2e} {sig} | effect={row['effect_size']:.3f} | "
              f"n_sub={row['n_substrates']:.0f}")
    print()

    # Analysis 2: Bin width
    print("2. Size-Matching Bin Width Sensitivity")
    print("   ------------------------------------")
    for _, row in bw_res.iterrows():
        sig = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01 else
               ("*" if row["p_value"] < 0.05 else "NS"))
        print(f"   Bin={row['value']:.0f} kDa | OR={row['effect_size']:.2f} | "
              f"p={row['p_value']:.2e} {sig} | n_bg={row['n_background']:.0f}")
    print()

    # Analysis 3: Multiplier
    print("3. Background Multiplier Sensitivity")
    print("   ----------------------------------")
    for _, row in mult_res.iterrows():
        sig = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01 else
               ("*" if row["p_value"] < 0.05 else "NS"))
        print(f"   Mult={row['value']:.0f}x | OR={row['effect_size']:.2f} | "
              f"p={row['p_value']:.2e} {sig} | n_bg={row['n_background']:.0f}")
    print()

    # Overall robustness assessment
    print("ROBUSTNESS ASSESSMENT:")
    print("-" * 40)

    # Check if all SILAC thresholds give significant results
    silac_multi_all_sig = all(silac_multi["p_value"] < 0.05)
    silac_nvc_all_sig = all(silac_nvc["p_value"].dropna() < 0.05)
    bw_all_sig = all(bw_res["p_value"] < 0.05)
    mult_all_sig = all(mult_res["p_value"] < 0.05)

    print(f"  HSP60 multi-domain enrichment robust across SILAC thresholds: "
          f"{'YES' if silac_multi_all_sig else 'PARTIALLY (some thresholds NS)'}")
    print(f"  HSP60 N>C contact order robust across SILAC thresholds:       "
          f"{'YES' if silac_nvc_all_sig else 'PARTIALLY (some thresholds NS)'}")
    print(f"  GroEL TIM barrel enrichment robust across bin widths:          "
          f"{'YES' if bw_all_sig else 'PARTIALLY (some bin widths NS)'}")
    print(f"  GroEL TIM barrel enrichment robust across multipliers:         "
          f"{'YES' if mult_all_sig else 'PARTIALLY (some multipliers NS)'}")
    print()

    # p-value range summary
    for label, subset in [("SILAC multi-domain", silac_multi),
                          ("SILAC N>C CO", silac_nvc),
                          ("Bin width TIM", bw_res),
                          ("Multiplier TIM", mult_res)]:
        pvals = subset["p_value"].dropna()
        if len(pvals) > 0:
            print(f"  {label:<25s}: p-range [{pvals.min():.2e}, {pvals.max():.2e}]")

    print()
    print(f"Output files:")
    print(f"  {OUT_TSV}")
    print(f"  {OUT_PDF}")
    print(f"  {OUT_PNG}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
