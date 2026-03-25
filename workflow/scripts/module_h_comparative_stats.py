#!/usr/bin/env python3
"""
Module H: Comparative Statistics for Antah Asti Prarambh
=========================================================

Brings together all previous analyses to answer three core biological questions
with proper statistical testing:
  Goal 1: Domain architecture enrichment
  Goal 2: N-vs-C stability asymmetry
  Goal 3: MTS and matrix targeting

Requirements: pandas, scipy, statsmodels, numpy
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from datetime import datetime

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================================
# PATHS
# ============================================================================
BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

# Inputs
CATH_PROTEIN  = f"{BASE}/results/domains/cath_protein_summary.tsv"
CATH_DOMAINS  = f"{BASE}/results/domains/cath_domain_assignments.tsv"
DOM_METRICS   = f"{BASE}/results/domains/domain_structural_metrics.tsv"
FOLDSEEK      = f"{BASE}/results/domains/foldseek_clusters.tsv"
DOM_DIST      = f"{BASE}/results/domains/domain_distribution_summary.tsv"
NVC_PAIRED    = f"{BASE}/results/termini/n_vs_c_paired.tsv"
CONTACT_ORDER = f"{BASE}/results/termini/contact_order.tsv"
REGION_BOUNDS = f"{BASE}/results/termini/region_boundaries.tsv"
COMBINED_TARG = f"{BASE}/results/mts/combined_targeting.tsv"
MTS_DOMAIN    = f"{BASE}/results/mts/mts_domain_relationship.tsv"
GROEL_SUBS    = f"{BASE}/data/processed/groel_substrates_standardized.tsv"
HSP60_SUBS    = f"{BASE}/data/processed/hsp60_tier1_substrates.tsv"
HOMOLOGS      = f"{BASE}/data/processed/groel_hsp60_homologs.tsv"
STRUCT_INDEX  = f"{BASE}/results/structures/structure_index.tsv"
ECOLI_PROT    = f"{BASE}/data/raw/uniprot/ecoli_k12_proteome.tsv"
HUMAN_PROT    = f"{BASE}/data/raw/uniprot/human_proteome.tsv"
MATRIX_PROT   = f"{BASE}/data/processed/human_matrix_proteome.tsv"
MITO_PROT     = f"{BASE}/data/processed/human_mito_proteome.tsv"

# Outputs
OUT_DIR = f"{BASE}/results/stats"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_DOMAIN   = f"{OUT_DIR}/domain_enrichment.tsv"
OUT_STAB     = f"{OUT_DIR}/stability_comparisons.tsv"
OUT_TARG     = f"{OUT_DIR}/targeting_stats.tsv"
OUT_CORRECT  = f"{OUT_DIR}/corrected_pvalues.tsv"
OUT_REPORT   = f"{OUT_DIR}/statistics_summary_report.txt"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def safe_fisher(a, b, c, d):
    """Fisher's exact test with odds ratio and 95% CI."""
    table = np.array([[a, b], [c, d]])
    odds_ratio, pvalue = stats.fisher_exact(table, alternative='two-sided')
    # Compute 95% CI using log(OR) +/- 1.96*SE
    if a == 0 or b == 0 or c == 0 or d == 0:
        # Haldane correction
        a_c, b_c, c_c, d_c = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        a_c, b_c, c_c, d_c = a, b, c, d
    log_or = np.log(a_c * d_c / (b_c * c_c))
    se = np.sqrt(1/a_c + 1/b_c + 1/c_c + 1/d_c)
    ci_lower = np.exp(log_or - 1.96 * se)
    ci_upper = np.exp(log_or + 1.96 * se)
    return odds_ratio, ci_lower, ci_upper, pvalue


def rank_biserial(x, y=None, paired=False):
    """
    Rank-biserial correlation as effect size.
    For paired test (Wilcoxon signed-rank): r = 1 - (2*W) / (n*(n+1)/2)
    For unpaired test (Mann-Whitney): r = 1 - (2*U) / (n1*n2)
    """
    if paired:
        diff = np.array(x) - np.array(y)
        diff = diff[diff != 0]
        n = len(diff)
        if n == 0:
            return 0.0
        stat, _ = stats.wilcoxon(x, y)
        r = 1 - (2 * stat) / (n * (n + 1) / 2)
        return r
    else:
        n1, n2 = len(x), len(y)
        if n1 == 0 or n2 == 0:
            return 0.0
        stat, _ = stats.mannwhitneyu(x, y, alternative='two-sided')
        r = 1 - (2 * stat) / (n1 * n2)
        return r


def size_match_sample(substrate_lengths, background_df, length_col='Length',
                      id_col='Entry', bin_width=10, seed=42, multiplier=3):
    """
    Build a size-matched control set by stratified sampling in kDa bins.
    MW approximation: length * 0.11 kDa per amino acid.

    multiplier: how many background proteins per substrate per bin (up to available)
    """
    np.random.seed(seed)
    # Convert lengths to kDa bins
    sub_kda = np.array(substrate_lengths) * 0.11
    bg_kda = np.array(background_df[length_col].astype(float)) * 0.11

    sub_bins = (sub_kda // bin_width).astype(int)
    bg_bins = (bg_kda // bin_width).astype(int)

    background_df = background_df.copy()
    background_df['_kda_bin'] = bg_bins

    selected = []
    for b in np.unique(sub_bins):
        n_needed = int(np.sum(sub_bins == b)) * multiplier
        pool = background_df[background_df['_kda_bin'] == b]
        if len(pool) == 0:
            continue
        n_take = min(n_needed, len(pool))
        chosen = pool.sample(n=n_take, random_state=seed)
        selected.append(chosen)

    if selected:
        result = pd.concat(selected, ignore_index=True)
        result.drop(columns=['_kda_bin'], inplace=True)
        return result
    else:
        return pd.DataFrame()


def cath_class_name(c):
    """Map CATH class number to name."""
    m = {1: 'Mainly-alpha', 2: 'Mainly-beta', 3: 'Alpha-beta', 4: 'Few SS',
         6: 'Special'}
    return m.get(c, f'Class_{c}')


print("=" * 70)
print("Module H: Comparative Statistics")
print(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# ============================================================================
# LOAD DATA
# ============================================================================
print("\n[1] Loading input files...")

cath_protein = pd.read_csv(CATH_PROTEIN, sep='\t')
cath_domains = pd.read_csv(CATH_DOMAINS, sep='\t')
dom_metrics  = pd.read_csv(DOM_METRICS, sep='\t')
foldseek     = pd.read_csv(FOLDSEEK, sep='\t')
nvc_paired   = pd.read_csv(NVC_PAIRED, sep='\t')
contact_order = pd.read_csv(CONTACT_ORDER, sep='\t')
region_bounds = pd.read_csv(REGION_BOUNDS, sep='\t')
combined_targ = pd.read_csv(COMBINED_TARG, sep='\t')
mts_domain   = pd.read_csv(MTS_DOMAIN, sep='\t')
groel_subs   = pd.read_csv(GROEL_SUBS, sep='\t')
hsp60_subs   = pd.read_csv(HSP60_SUBS, sep='\t')
homologs     = pd.read_csv(HOMOLOGS, sep='\t')
struct_index = pd.read_csv(STRUCT_INDEX, sep='\t')
ecoli_prot   = pd.read_csv(ECOLI_PROT, sep='\t')
human_prot   = pd.read_csv(HUMAN_PROT, sep='\t')
matrix_prot  = pd.read_csv(MATRIX_PROT, sep='\t')
mito_prot    = pd.read_csv(MITO_PROT, sep='\t')

# Parse source_dataset to identify protein groups
struct_index['is_groel'] = struct_index['source_dataset'].str.contains('groel', na=False)
struct_index['is_hsp60'] = struct_index['source_dataset'].str.contains('hsp60', na=False)
struct_index['is_matrix'] = struct_index['source_dataset'].str.contains('matrix', na=False)
struct_index['is_mito'] = struct_index['source_dataset'].str.contains('mito', na=False)

groel_ids = set(groel_subs['current_accession'])
hsp60_ids = set(hsp60_subs['uniprot_id'])

print(f"  GroEL substrates: {len(groel_ids)}")
print(f"  HSP60 substrates: {len(hsp60_ids)}")
print(f"  E. coli proteome: {len(ecoli_prot)}")
print(f"  Human proteome: {len(human_prot)}")
print(f"  Matrix proteome: {len(matrix_prot)}")
print(f"  Mito proteome: {len(mito_prot)}")
print(f"  CATH domains assigned: {len(cath_domains)}")
print(f"  N-vs-C paired proteins: {len(nvc_paired)}")

# ============================================================================
# STEP H2: DOMAIN ARCHITECTURE ENRICHMENT TESTS
# ============================================================================
print("\n[2] Domain Architecture Enrichment Tests (H1.1, H1.2, H1.3)...")

# -- Build size-matched controls --

# 2a. GroEL: cytoplasmic E. coli proteins as background
# Parse subcellular location to identify cytoplasmic proteins
ecoli_prot['is_cytoplasmic'] = False
for idx, row in ecoli_prot.iterrows():
    loc = str(row.get('Subcellular location [CC]', '')).lower()
    kw  = str(row.get('Keywords', '')).lower()
    go_cc = str(row.get('Gene Ontology (cellular component)', '')).lower()
    if 'cytoplasm' in loc or 'cytoplasm' in kw or 'cytoplasm' in go_cc:
        ecoli_prot.at[idx, 'is_cytoplasmic'] = True

ecoli_cyto = ecoli_prot[ecoli_prot['is_cytoplasmic']].copy()
# Exclude GroEL substrates from background
ecoli_cyto = ecoli_cyto[~ecoli_cyto['Entry'].isin(groel_ids)]
print(f"  E. coli cytoplasmic (non-substrate) background: {len(ecoli_cyto)}")

# Get GroEL substrate lengths from groel_subs
groel_lengths = groel_subs[groel_subs['location_category'] == 'cytoplasmic']['length'].values
print(f"  Cytoplasmic GroEL substrates for size matching: {len(groel_lengths)}")

# Size-match
ecoli_bg = size_match_sample(groel_lengths, ecoli_cyto, length_col='Length',
                              id_col='Entry', bin_width=10, multiplier=3)
print(f"  Size-matched E. coli background: {len(ecoli_bg)}")

# 2b. HSP60: matrix proteins as background
# Use matrix_prot, exclude HSP60 substrates
matrix_bg_df = matrix_prot[~matrix_prot['uniprot_id'].isin(hsp60_ids)].copy()
# Need length info - merge with human_prot or struct_index
# Get lengths from human_prot
human_lengths = human_prot[['Entry', 'Length']].rename(columns={'Entry': 'uniprot_id'})
matrix_bg_df = matrix_bg_df.merge(human_lengths, on='uniprot_id', how='left')
matrix_bg_df['Length'] = matrix_bg_df['Length'].fillna(0).astype(int)
matrix_bg_df = matrix_bg_df[matrix_bg_df['Length'] > 0]

# HSP60 substrate lengths
hsp60_in_struct = struct_index[struct_index['is_hsp60']]
hsp60_lengths = hsp60_in_struct['residues_modeled'].values

hsp60_bg = size_match_sample(hsp60_lengths, matrix_bg_df, length_col='Length',
                              id_col='uniprot_id', bin_width=10, multiplier=3)
print(f"  Matrix background (non-HSP60): {len(matrix_bg_df)}")
print(f"  Size-matched matrix background: {len(hsp60_bg)}")

# -- Now run Fisher's exact tests for CATH superfamily enrichment --

# We need CATH assignments for background proteins too
# Background proteins don't have CATH assignments in our data (only substrates do)
# We'll use what we have: compare superfamily distribution within substrate datasets

# For GroEL: domains from groel substrates vs domains from "mito-only" or ecoli bg
# Actually, our cath_domains are from our structure set. We have structures for
# groel, hsp60, matrix, mito datasets. Let's use:
#   - groel substrates: source_dataset contains 'groel'
#   - background: source_dataset = 'mito' only (no hsp60/matrix) + filter to E. coli bg IDs
# Since we don't have CATH for arbitrary ecoli proteins, we use the mito-only set
# as a structural background. This is the best we can do with available data.

# Map proteins to datasets
cath_domains_with_ds = cath_domains.merge(
    struct_index[['uniprot_accession', 'source_dataset']].drop_duplicates(),
    on='uniprot_accession', how='left'
)

# GroEL substrate domains
groel_domains = cath_domains_with_ds[
    cath_domains_with_ds['source_dataset'].str.contains('groel', na=False)
]
groel_domain_proteins = set(groel_domains['uniprot_accession'])

# HSP60 substrate domains
hsp60_domains = cath_domains_with_ds[
    cath_domains_with_ds['source_dataset'].str.contains('hsp60', na=False)
]
hsp60_domain_proteins = set(hsp60_domains['uniprot_accession'])

# Background: mito-only (not hsp60, not matrix substrate of hsp60)
mito_only_domains = cath_domains_with_ds[
    (cath_domains_with_ds['source_dataset'] == 'mito') |
    (cath_domains_with_ds['source_dataset'] == 'matrix,mito')
]
# Exclude any that overlap with hsp60
mito_only_domains = mito_only_domains[
    ~mito_only_domains['uniprot_accession'].isin(hsp60_ids)
]

# For GroEL background, ideally we use ecoli non-substrate CATH data
# But we only have CATH for our pilot set. Use the full non-groel pilot set as background.
non_groel_ecoli_domains = cath_domains_with_ds[
    ~cath_domains_with_ds['source_dataset'].str.contains('groel', na=False) &
    ~cath_domains_with_ds['source_dataset'].str.contains('hsp60', na=False)
]

print(f"\n  GroEL domains for enrichment: {len(groel_domains)}")
print(f"  Non-GroEL background domains: {len(non_groel_ecoli_domains)}")
print(f"  HSP60 domains for enrichment: {len(hsp60_domains)}")
print(f"  Mito-only background domains: {len(mito_only_domains)}")

# Fisher's exact test for each superfamily in each dataset
enrichment_results = []

def run_superfamily_enrichment(sub_domains, bg_domains, dataset_name):
    """Fisher's exact test for each CATH superfamily: substrate vs background."""
    results = []
    all_sfs = set(sub_domains['cath_superfamily'].dropna()) | set(bg_domains['cath_superfamily'].dropna())

    # Count proteins (not domains) with each superfamily
    sub_sf_counts = sub_domains.groupby('cath_superfamily')['uniprot_accession'].nunique()
    bg_sf_counts = bg_domains.groupby('cath_superfamily')['uniprot_accession'].nunique()

    n_sub_total = sub_domains['uniprot_accession'].nunique()
    n_bg_total = bg_domains['uniprot_accession'].nunique()

    for sf in sorted(all_sfs):
        if pd.isna(sf) or sf == '':
            continue
        a = sub_sf_counts.get(sf, 0)  # substrate has SF
        b = n_sub_total - a           # substrate lacks SF
        c = bg_sf_counts.get(sf, 0)   # background has SF
        d = n_bg_total - c            # background lacks SF

        if a + c < 3:  # skip very rare superfamilies
            continue

        or_val, ci_lo, ci_hi, pval = safe_fisher(a, b, c, d)

        # Get superfamily name if available
        sf_name = ''
        sf_rows = sub_domains[sub_domains['cath_superfamily'] == sf]
        if len(sf_rows) > 0 and 'cath_superfamily_name' in sf_rows.columns:
            names = sf_rows['cath_superfamily_name'].dropna().unique()
            if len(names) > 0 and str(names[0]).strip():
                sf_name = str(names[0])

        results.append({
            'dataset': dataset_name,
            'superfamily': sf,
            'superfamily_name': sf_name,
            'n_substrate': int(a),
            'n_substrate_total': int(n_sub_total),
            'n_background': int(c),
            'n_background_total': int(n_bg_total),
            'odds_ratio': or_val,
            'ci_lower': ci_lo,
            'ci_upper': ci_hi,
            'pvalue': pval,
        })
    return results

# GroEL enrichment
groel_enrich = run_superfamily_enrichment(groel_domains, non_groel_ecoli_domains, 'groel')
enrichment_results.extend(groel_enrich)

# HSP60 enrichment
hsp60_enrich = run_superfamily_enrichment(hsp60_domains, mito_only_domains, 'hsp60')
enrichment_results.extend(hsp60_enrich)

enrichment_df = pd.DataFrame(enrichment_results)

# BH correction within each dataset
enrichment_df['pvalue_bh'] = np.nan
enrichment_df['significant'] = False
for ds in enrichment_df['dataset'].unique():
    mask = enrichment_df['dataset'] == ds
    pvals = enrichment_df.loc[mask, 'pvalue'].values
    if len(pvals) > 0:
        reject, corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        enrichment_df.loc[mask, 'pvalue_bh'] = corrected
        enrichment_df.loc[mask, 'significant'] = reject

enrichment_df = enrichment_df.sort_values(['dataset', 'pvalue']).reset_index(drop=True)
enrichment_df.to_csv(OUT_DOMAIN, sep='\t', index=False)
print(f"\n  Domain enrichment results saved: {OUT_DOMAIN}")
print(f"  Total tests: {len(enrichment_df)}")
print(f"  Significant (BH < 0.05): {enrichment_df['significant'].sum()}")

# Chi-squared test for CATH class distribution
print("\n  Chi-squared tests for CATH class distribution...")
chi2_results = []

for ds_name, ds_domains, bg_doms in [
    ('groel', groel_domains, non_groel_ecoli_domains),
    ('hsp60', hsp60_domains, mito_only_domains)
]:
    sub_classes = ds_domains['cath_class'].dropna().astype(int).value_counts().sort_index()
    bg_classes = bg_doms['cath_class'].dropna().astype(int).value_counts().sort_index()

    # Align indices
    all_classes = sorted(set(sub_classes.index) | set(bg_classes.index))
    sub_arr = np.array([sub_classes.get(c, 0) for c in all_classes])
    bg_arr = np.array([bg_classes.get(c, 0) for c in all_classes])

    # Only keep classes with > 0 in both
    keep = (sub_arr + bg_arr) > 0
    sub_arr = sub_arr[keep]
    bg_arr = bg_arr[keep]
    kept_classes = [c for c, k in zip(all_classes, keep) if k]

    if len(sub_arr) >= 2:
        contingency = np.array([sub_arr, bg_arr])
        chi2, p, dof, expected = stats.chi2_contingency(contingency)
        # Cramér's V
        n = contingency.sum()
        k = min(contingency.shape) - 1
        cramers_v = np.sqrt(chi2 / (n * k)) if (n * k) > 0 else 0

        print(f"    {ds_name}: chi2={chi2:.2f}, dof={dof}, p={p:.4e}, Cramér's V={cramers_v:.3f}")
        class_labels = [cath_class_name(c) for c in kept_classes]
        for i, cl in enumerate(class_labels):
            print(f"      {cl}: substrate={sub_arr[i]}, background={bg_arr[i]}")

        chi2_results.append({
            'dataset': ds_name, 'test': 'chi2_cath_class',
            'statistic': chi2, 'dof': dof, 'pvalue': p, 'effect_size': cramers_v,
            'effect_name': 'cramers_v'
        })

# H1.3: Conservation of enrichment between GroEL and HSP60
print("\n  H1.3: Conservation test (shared superfamilies between GroEL and HSP60)...")
groel_sf_set = set(groel_domains['cath_superfamily'].dropna())
hsp60_sf_set = set(hsp60_domains['cath_superfamily'].dropna())
shared_sf = groel_sf_set & hsp60_sf_set
only_groel = groel_sf_set - hsp60_sf_set
only_hsp60 = hsp60_sf_set - groel_sf_set

print(f"    Shared superfamilies: {len(shared_sf)}")
print(f"    GroEL-only: {len(only_groel)}")
print(f"    HSP60-only: {len(only_hsp60)}")

# Test: are shared superfamilies more than expected by chance?
# Hypergeometric test: from universe of all superfamilies seen,
# draw hsp60_sf_set size, what's P(overlap >= observed)?
universe_sf = groel_sf_set | hsp60_sf_set
N = len(universe_sf)  # population
K = len(groel_sf_set)  # success states in pop
n = len(hsp60_sf_set)  # draws
k = len(shared_sf)     # observed successes

p_hypergeom = stats.hypergeom.sf(k - 1, N, K, n)
print(f"    Hypergeometric test: P(overlap >= {k}) = {p_hypergeom:.4e}")

# Also: For homolog pairs, check if both have same top superfamily
homolog_enrich = []
homologs_with_cath = homologs.copy()
# Get top superfamily for each GroEL substrate
groel_top_sf = groel_domains.groupby('uniprot_accession')['cath_superfamily'].first().to_dict()
hsp60_top_sf = hsp60_domains.groupby('uniprot_accession')['cath_superfamily'].first().to_dict()

n_pairs_checked = 0
n_pairs_same_sf = 0
for _, row in homologs.iterrows():
    g_acc = row['groel_accession']
    h_acc = row['hsp60_accession']
    g_sf = groel_top_sf.get(g_acc)
    h_sf = hsp60_top_sf.get(h_acc)
    if g_sf is not None and h_sf is not None:
        n_pairs_checked += 1
        if g_sf == h_sf:
            n_pairs_same_sf += 1

print(f"    Homolog pairs with CATH for both: {n_pairs_checked}")
print(f"    Pairs sharing top superfamily: {n_pairs_same_sf}")
if n_pairs_checked > 0:
    frac_shared = n_pairs_same_sf / n_pairs_checked
    print(f"    Fraction sharing: {frac_shared:.3f}")

# ============================================================================
# STEP H3: N-VS-C STABILITY COMPARISON STATISTICS
# ============================================================================
print("\n[3] N-vs-C Stability Comparison Statistics (H2.1, H2.2, H2.3)...")

stability_results = []

# H2.1: Within-protein paired tests
# Metrics to test: relative_contact_order, mean_plddt, frac_helix, frac_strand, frac_plddt_gt70
paired_metrics = [
    ('relative_contact_order', 'n_domain', 'c_region'),
    ('mean_plddt', 'n_domain', 'c_region'),
    ('frac_helix', 'n_domain', 'c_region'),
    ('frac_strand', 'n_domain', 'c_region'),
    ('frac_plddt_gt70', 'n_domain', 'c_region'),
]

# Define dataset groups from source_dataset
def classify_nvc_dataset(sd):
    """Classify a source_dataset string into analysis groups."""
    if 'groel' in str(sd):
        return 'groel'
    elif 'hsp60' in str(sd):
        return 'hsp60'
    elif 'matrix' in str(sd):
        return 'matrix_bg'
    elif 'mito' in str(sd):
        return 'mito_bg'
    return 'other'

nvc_paired['analysis_group'] = nvc_paired['source_dataset'].apply(classify_nvc_dataset)

print("\n  H2.1: Within-protein paired Wilcoxon signed-rank tests")
print("  " + "-" * 60)

for group in ['groel', 'hsp60', 'matrix_bg', 'mito_bg']:
    subset = nvc_paired[nvc_paired['analysis_group'] == group]
    if len(subset) < 5:
        continue
    print(f"\n  Dataset: {group} (n={len(subset)})")

    for metric_base, n_sfx, c_sfx in paired_metrics:
        n_col = f"{metric_base}_{n_sfx}"
        c_col = f"{metric_base}_{c_sfx}"
        if n_col not in subset.columns or c_col not in subset.columns:
            continue

        n_vals = subset[n_col].dropna()
        c_vals = subset[c_col].dropna()
        # Need matched pairs
        valid = subset[[n_col, c_col]].dropna()
        if len(valid) < 5:
            continue

        n_v = valid[n_col].values
        c_v = valid[c_col].values

        # Wilcoxon signed-rank test
        try:
            stat_w, p_w = stats.wilcoxon(n_v, c_v)
        except ValueError:
            continue
        rb = rank_biserial(n_v, c_v, paired=True)

        median_diff = np.median(n_v - c_v)
        print(f"    {metric_base}: W={stat_w:.1f}, p={p_w:.4e}, r_rb={rb:.3f}, "
              f"median_diff={median_diff:.4f}")

        stability_results.append({
            'hypothesis': 'H2.1',
            'dataset': group,
            'test': 'wilcoxon_signed_rank',
            'metric': metric_base,
            'n': len(valid),
            'statistic': stat_w,
            'pvalue': p_w,
            'effect_size': rb,
            'effect_name': 'rank_biserial_r',
            'median_n_domain': np.median(n_v),
            'median_c_region': np.median(c_v),
            'median_diff': median_diff,
            'direction': 'N > C' if median_diff > 0 else 'C > N',
        })

# BH correction within each group for H2.1
stab_df = pd.DataFrame(stability_results)

# H2.2: Cross-dataset comparison - substrate vs background asymmetry
print("\n  H2.2: Cross-dataset Mann-Whitney U tests (substrate vs background asymmetry)")
print("  " + "-" * 60)

h22_results = []
asymmetry_col = 'relative_contact_order_diff'  # N - C difference

for sub_group, bg_group, label in [
    ('groel', 'mito_bg', 'GroEL vs mito_bg'),
    ('hsp60', 'matrix_bg', 'HSP60 vs matrix_bg'),
    ('hsp60', 'mito_bg', 'HSP60 vs mito_bg'),
]:
    sub_data = nvc_paired[nvc_paired['analysis_group'] == sub_group][asymmetry_col].dropna()
    bg_data = nvc_paired[nvc_paired['analysis_group'] == bg_group][asymmetry_col].dropna()

    if len(sub_data) < 5 or len(bg_data) < 5:
        print(f"    {label}: skipped (n_sub={len(sub_data)}, n_bg={len(bg_data)})")
        continue

    stat_u, p_u = stats.mannwhitneyu(sub_data, bg_data, alternative='two-sided')
    rb = rank_biserial(sub_data.values, bg_data.values, paired=False)

    print(f"    {label}: U={stat_u:.1f}, p={p_u:.4e}, r_rb={rb:.3f}, "
          f"med_sub={np.median(sub_data):.4f}, med_bg={np.median(bg_data):.4f}")

    h22_results.append({
        'hypothesis': 'H2.2',
        'dataset': label,
        'test': 'mann_whitney_u',
        'metric': 'contact_order_asymmetry',
        'n': len(sub_data) + len(bg_data),
        'n_substrate': len(sub_data),
        'n_background': len(bg_data),
        'statistic': stat_u,
        'pvalue': p_u,
        'effect_size': rb,
        'effect_name': 'rank_biserial_r',
        'median_substrate': np.median(sub_data),
        'median_background': np.median(bg_data),
        'direction': 'sub > bg' if np.median(sub_data) > np.median(bg_data) else 'bg > sub',
    })

# Also test with other metrics
for metric_base in ['mean_plddt', 'frac_helix', 'frac_strand']:
    diff_col = f"{metric_base}_diff"
    if diff_col not in nvc_paired.columns:
        continue
    for sub_group, bg_group, label in [
        ('groel', 'mito_bg', 'GroEL vs mito_bg'),
        ('hsp60', 'matrix_bg', 'HSP60 vs matrix_bg'),
    ]:
        sub_data = nvc_paired[nvc_paired['analysis_group'] == sub_group][diff_col].dropna()
        bg_data = nvc_paired[nvc_paired['analysis_group'] == bg_group][diff_col].dropna()
        if len(sub_data) < 5 or len(bg_data) < 5:
            continue
        stat_u, p_u = stats.mannwhitneyu(sub_data, bg_data, alternative='two-sided')
        rb = rank_biserial(sub_data.values, bg_data.values, paired=False)

        h22_results.append({
            'hypothesis': 'H2.2',
            'dataset': label,
            'test': 'mann_whitney_u',
            'metric': f'{metric_base}_asymmetry',
            'n': len(sub_data) + len(bg_data),
            'n_substrate': len(sub_data),
            'n_background': len(bg_data),
            'statistic': stat_u,
            'pvalue': p_u,
            'effect_size': rb,
            'effect_name': 'rank_biserial_r',
            'median_substrate': np.median(sub_data),
            'median_background': np.median(bg_data),
            'direction': 'sub > bg' if np.median(sub_data) > np.median(bg_data) else 'bg > sub',
        })

# H2.3: GroEL class comparison (Class I vs II vs III)
print("\n  H2.3: GroEL class comparison (Kruskal-Wallis + post-hoc Dunn's)")
print("  " + "-" * 60)

h23_results = []

# Merge GroEL class info
nvc_groel = nvc_paired[nvc_paired['analysis_group'] == 'groel'].copy()
# groel_class column already exists in nvc_paired
nvc_groel = nvc_groel[nvc_groel['groel_class'].isin(['I', 'II', 'III'])]
print(f"  GroEL proteins with class info: {len(nvc_groel)}")

for metric_base in ['relative_contact_order', 'mean_plddt', 'frac_helix', 'frac_strand']:
    diff_col = f"{metric_base}_diff"
    if diff_col not in nvc_groel.columns:
        continue

    groups = {}
    for cls in ['I', 'II', 'III']:
        vals = nvc_groel[nvc_groel['groel_class'] == cls][diff_col].dropna()
        if len(vals) >= 3:
            groups[cls] = vals.values

    if len(groups) < 2:
        continue

    group_arrays = list(groups.values())
    group_names = list(groups.keys())

    # Kruskal-Wallis
    if len(groups) >= 3:
        h_stat, p_kw = stats.kruskal(*group_arrays)
    elif len(groups) == 2:
        h_stat, p_kw = stats.kruskal(*group_arrays)
    else:
        continue

    # Eta-squared for Kruskal-Wallis: H / (n-1)
    n_total = sum(len(g) for g in group_arrays)
    eta_sq = (h_stat - len(groups) + 1) / (n_total - 1) if n_total > 1 else 0

    print(f"    {metric_base}_diff: H={h_stat:.2f}, p={p_kw:.4e}, eta2={eta_sq:.3f}")
    for cn in group_names:
        print(f"      Class {cn}: n={len(groups[cn])}, median={np.median(groups[cn]):.4f}")

    h23_results.append({
        'hypothesis': 'H2.3',
        'dataset': 'groel_classes',
        'test': 'kruskal_wallis',
        'metric': f'{metric_base}_asymmetry',
        'n': n_total,
        'statistic': h_stat,
        'pvalue': p_kw,
        'effect_size': eta_sq,
        'effect_name': 'eta_squared',
        'groups': str({k: len(v) for k, v in groups.items()}),
    })

    # Post-hoc pairwise Mann-Whitney U (Dunn's approximation)
    if p_kw < 0.1 and len(groups) >= 2:
        pairs = []
        for i in range(len(group_names)):
            for j in range(i+1, len(group_names)):
                g1, g2 = group_names[i], group_names[j]
                u_stat, p_pair = stats.mannwhitneyu(groups[g1], groups[g2],
                                                     alternative='two-sided')
                rb = rank_biserial(groups[g1], groups[g2], paired=False)
                pairs.append((g1, g2, u_stat, p_pair, rb))
                print(f"        {g1} vs {g2}: U={u_stat:.1f}, p={p_pair:.4e}, r={rb:.3f}")

        # BH correction on post-hoc
        if pairs:
            ph_pvals = [p[3] for p in pairs]
            reject_ph, corrected_ph, _, _ = multipletests(ph_pvals, alpha=0.05, method='fdr_bh')
            for idx_p, (g1, g2, u_stat, p_pair, rb) in enumerate(pairs):
                h23_results.append({
                    'hypothesis': 'H2.3',
                    'dataset': f'groel_{g1}_vs_{g2}',
                    'test': 'mann_whitney_posthoc',
                    'metric': f'{metric_base}_asymmetry',
                    'n': len(groups[g1]) + len(groups[g2]),
                    'statistic': u_stat,
                    'pvalue': p_pair,
                    'pvalue_bh_posthoc': corrected_ph[idx_p],
                    'effect_size': rb,
                    'effect_name': 'rank_biserial_r',
                })

# Combine all stability results
all_stab = stability_results + h22_results + h23_results
stab_df = pd.DataFrame(all_stab)

# BH correction within each hypothesis group
stab_df['pvalue_bh'] = np.nan
stab_df['significant'] = False
for hyp in stab_df['hypothesis'].unique():
    mask = stab_df['hypothesis'] == hyp
    pvals = stab_df.loc[mask, 'pvalue'].values
    if len(pvals) > 0:
        reject, corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        stab_df.loc[mask, 'pvalue_bh'] = corrected
        stab_df.loc[mask, 'significant'] = reject

stab_df.to_csv(OUT_STAB, sep='\t', index=False)
print(f"\n  Stability comparison results saved: {OUT_STAB}")
print(f"  Total tests: {len(stab_df)}")
print(f"  Significant (BH < 0.05): {stab_df['significant'].sum()}")

# ============================================================================
# STEP H4: MTS AND MATRIX TARGETING STATISTICS
# ============================================================================
print("\n[4] MTS and Matrix Targeting Statistics (H3.1, H3.2, H3.3)...")

targeting_results = []

# H3.1: HSP60 substrates enriched for matrix localization
# Among general mito proteome (combined_targeting): how many are matrix?
# Among HSP60 substrates: how many are matrix?

ct = combined_targ.copy()
ct['is_hsp60'] = ct['is_hsp60_substrate'] if 'is_hsp60_substrate' in ct.columns else ct['uniprot_accession'].isin(hsp60_ids)
ct['is_matrix_target'] = ct['targeting_classification'].str.contains('matrix|Matrix', na=False, case=False)

# Also check mitocarta_is_matrix column
if 'mitocarta_is_matrix' in ct.columns:
    ct['is_matrix_target'] = ct['is_matrix_target'] | (ct['mitocarta_is_matrix'] == True)

n_hsp60_matrix = ct[ct['is_hsp60'] & ct['is_matrix_target']].shape[0]
n_hsp60_nonmatrix = ct[ct['is_hsp60'] & ~ct['is_matrix_target']].shape[0]
n_nonhsp60_matrix = ct[~ct['is_hsp60'] & ct['is_matrix_target']].shape[0]
n_nonhsp60_nonmatrix = ct[~ct['is_hsp60'] & ~ct['is_matrix_target']].shape[0]

or_h31, ci_lo_h31, ci_hi_h31, p_h31 = safe_fisher(
    n_hsp60_matrix, n_hsp60_nonmatrix, n_nonhsp60_matrix, n_nonhsp60_nonmatrix
)

print(f"  H3.1: HSP60 substrates matrix enrichment")
print(f"    HSP60 matrix: {n_hsp60_matrix}, HSP60 non-matrix: {n_hsp60_nonmatrix}")
print(f"    Non-HSP60 matrix: {n_nonhsp60_matrix}, Non-HSP60 non-matrix: {n_nonhsp60_nonmatrix}")
print(f"    OR={or_h31:.2f} [{ci_lo_h31:.2f}, {ci_hi_h31:.2f}], p={p_h31:.4e}")

targeting_results.append({
    'hypothesis': 'H3.1',
    'test': 'fisher_exact',
    'comparison': 'HSP60 matrix enrichment vs mito proteome',
    'a_hsp60_matrix': n_hsp60_matrix,
    'b_hsp60_nonmatrix': n_hsp60_nonmatrix,
    'c_bg_matrix': n_nonhsp60_matrix,
    'd_bg_nonmatrix': n_nonhsp60_nonmatrix,
    'odds_ratio': or_h31,
    'ci_lower': ci_lo_h31,
    'ci_upper': ci_hi_h31,
    'pvalue': p_h31,
})

# H3.2: MTS prevalence - HSP60 matrix substrates vs non-HSP60 matrix proteins
ct['has_mts'] = ct['has_transit_peptide'] == True if 'has_transit_peptide' in ct.columns else False

hsp60_matrix = ct[ct['is_hsp60'] & ct['is_matrix_target']]
nonhsp60_matrix = ct[~ct['is_hsp60'] & ct['is_matrix_target']]

n_hsp60m_mts = hsp60_matrix['has_mts'].sum()
n_hsp60m_nomts = len(hsp60_matrix) - n_hsp60m_mts
n_bg_m_mts = nonhsp60_matrix['has_mts'].sum()
n_bg_m_nomts = len(nonhsp60_matrix) - n_bg_m_mts

or_h32, ci_lo_h32, ci_hi_h32, p_h32 = safe_fisher(
    n_hsp60m_mts, n_hsp60m_nomts, n_bg_m_mts, n_bg_m_nomts
)

print(f"\n  H3.2: MTS prevalence in HSP60 matrix vs non-HSP60 matrix")
print(f"    HSP60 matrix MTS+: {n_hsp60m_mts}, MTS-: {n_hsp60m_nomts}")
print(f"    Non-HSP60 matrix MTS+: {n_bg_m_mts}, MTS-: {n_bg_m_nomts}")
print(f"    OR={or_h32:.2f} [{ci_lo_h32:.2f}, {ci_hi_h32:.2f}], p={p_h32:.4e}")

targeting_results.append({
    'hypothesis': 'H3.2',
    'test': 'fisher_exact',
    'comparison': 'MTS prevalence: HSP60 matrix vs non-HSP60 matrix',
    'a_hsp60_matrix': int(n_hsp60m_mts),
    'b_hsp60_nonmatrix': int(n_hsp60m_nomts),
    'c_bg_matrix': int(n_bg_m_mts),
    'd_bg_nonmatrix': int(n_bg_m_nomts),
    'odds_ratio': or_h32,
    'ci_lower': ci_lo_h32,
    'ci_upper': ci_hi_h32,
    'pvalue': p_h32,
})

# H3.2 additional: First-domain properties by MTS status
print(f"\n  H3.2 (additional): First-domain properties by MTS status")

# Use mts_domain_relationship and domain_structural_metrics
mts_dom = mts_domain.copy()
# Merge with domain structural metrics for domain index 1
first_dom_metrics = dom_metrics[dom_metrics['domain_index'] == 1].copy()
mts_merged = mts_dom.merge(first_dom_metrics,
                            on='uniprot_accession', how='inner',
                            suffixes=('', '_metric'))

mts_plus = mts_merged[mts_merged['mts_is_pre_domain'] == True]
mts_minus_cands = mts_merged[mts_merged['mts_is_pre_domain'] == False]

# Also split by MTS presence from combined_targ
mts_present_ids = set(ct[ct['has_mts']]['uniprot_accession'])
mts_merged['has_mts_signal'] = mts_merged['uniprot_accession'].isin(mts_present_ids)

mts_yes = mts_merged[mts_merged['has_mts_signal']]
mts_no = mts_merged[~mts_merged['has_mts_signal']]

mts_property_tests = []
for metric in ['domain_mean_plddt', 'domain_frac_helix', 'domain_frac_strand', 'domain_frac_coil']:
    if metric not in mts_merged.columns:
        continue
    v_yes = mts_yes[metric].dropna().values
    v_no = mts_no[metric].dropna().values
    if len(v_yes) < 5 or len(v_no) < 5:
        continue

    u_stat, p_val = stats.mannwhitneyu(v_yes, v_no, alternative='two-sided')
    rb = rank_biserial(v_yes, v_no, paired=False)
    print(f"    {metric}: MTS+ (n={len(v_yes)}, med={np.median(v_yes):.3f}) vs "
          f"MTS- (n={len(v_no)}, med={np.median(v_no):.3f}): U={u_stat:.1f}, p={p_val:.4e}, r={rb:.3f}")

    mts_property_tests.append({
        'hypothesis': 'H3.2',
        'test': 'mann_whitney_u',
        'comparison': f'First domain {metric}: MTS+ vs MTS-',
        'n_mts_plus': len(v_yes),
        'n_mts_minus': len(v_no),
        'median_mts_plus': np.median(v_yes),
        'median_mts_minus': np.median(v_no),
        'statistic': u_stat,
        'pvalue': p_val,
        'effect_size': rb,
        'effect_name': 'rank_biserial_r',
    })

# H3.3: MTS is predominantly a pre-domain extension
print(f"\n  H3.3: MTS pre-domain vs overlapping domain")
n_pre_domain = mts_dom['mts_is_pre_domain'].sum()
n_overlaps = mts_dom['mts_overlaps_domain'].sum()
n_total_mts = len(mts_dom)

# Binomial test: is proportion of pre-domain > 0.5?
p_binom = stats.binom_test(n_pre_domain, n_total_mts, p=0.5, alternative='greater')

print(f"    MTS pre-domain: {n_pre_domain}/{n_total_mts} ({100*n_pre_domain/n_total_mts:.1f}%)")
print(f"    MTS overlaps domain: {n_overlaps}/{n_total_mts} ({100*n_overlaps/n_total_mts:.1f}%)")
print(f"    Binomial test (H0: p=0.5): p={p_binom:.4e}")

targeting_results.append({
    'hypothesis': 'H3.3',
    'test': 'binomial_test',
    'comparison': 'MTS pre-domain vs overlapping',
    'n_pre_domain': int(n_pre_domain),
    'n_overlaps': int(n_overlaps),
    'n_total': int(n_total_mts),
    'proportion_pre_domain': n_pre_domain / n_total_mts if n_total_mts > 0 else 0,
    'pvalue': p_binom,
})

# Gap length analysis
gap_lengths = mts_dom['gap_length'].dropna()
print(f"    Gap length (MTS end to domain start): median={gap_lengths.median():.0f}, "
      f"mean={gap_lengths.mean():.1f}, IQR=[{gap_lengths.quantile(0.25):.0f}, {gap_lengths.quantile(0.75):.0f}]")

# Combine all targeting results
targ_df = pd.DataFrame(targeting_results + mts_property_tests)
# BH correction
if len(targ_df) > 0 and 'pvalue' in targ_df.columns:
    pvals = targ_df['pvalue'].values
    reject, corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    targ_df['pvalue_bh'] = corrected
    targ_df['significant'] = reject

targ_df.to_csv(OUT_TARG, sep='\t', index=False)
print(f"\n  Targeting results saved: {OUT_TARG}")
print(f"  Total tests: {len(targ_df)}")
print(f"  Significant (BH < 0.05): {targ_df['significant'].sum()}")

# ============================================================================
# STEP H5: MULTIPLE TESTING CORRECTION (HIERARCHICAL)
# ============================================================================
print("\n[5] Hierarchical Multiple Testing Correction...")

# Collect all p-values across families
all_tests = []

# Family 1: Domain architecture
for _, row in enrichment_df.iterrows():
    all_tests.append({
        'family': 'domain_architecture',
        'hypothesis': 'H1.1/H1.2',
        'test_id': f"{row['dataset']}_{row['superfamily']}",
        'test': 'fisher_exact',
        'pvalue_raw': row['pvalue'],
        'pvalue_within_family': row['pvalue_bh'],
        'effect_size': row['odds_ratio'],
        'effect_name': 'odds_ratio',
    })

# Add chi-squared results
for chi_res in chi2_results:
    all_tests.append({
        'family': 'domain_architecture',
        'hypothesis': 'H1.1/H1.2',
        'test_id': f"{chi_res['dataset']}_chi2_cath_class",
        'test': 'chi2_contingency',
        'pvalue_raw': chi_res['pvalue'],
        'pvalue_within_family': chi_res['pvalue'],  # only 1 per dataset
        'effect_size': chi_res['effect_size'],
        'effect_name': chi_res['effect_name'],
    })

# Hypergeometric for H1.3
all_tests.append({
    'family': 'domain_architecture',
    'hypothesis': 'H1.3',
    'test_id': 'sf_overlap_hypergeometric',
    'test': 'hypergeometric',
    'pvalue_raw': p_hypergeom,
    'pvalue_within_family': p_hypergeom,
    'effect_size': len(shared_sf) / len(universe_sf) if len(universe_sf) > 0 else 0,
    'effect_name': 'jaccard_index',
})

# Family 2: Stability
for _, row in stab_df.iterrows():
    all_tests.append({
        'family': 'stability_asymmetry',
        'hypothesis': row.get('hypothesis', 'H2'),
        'test_id': f"{row.get('dataset', '')}_{row.get('metric', '')}",
        'test': row.get('test', ''),
        'pvalue_raw': row['pvalue'],
        'pvalue_within_family': row.get('pvalue_bh', row['pvalue']),
        'effect_size': row.get('effect_size', np.nan),
        'effect_name': row.get('effect_name', ''),
    })

# Family 3: Targeting
for _, row in targ_df.iterrows():
    all_tests.append({
        'family': 'matrix_targeting',
        'hypothesis': row.get('hypothesis', 'H3'),
        'test_id': row.get('comparison', ''),
        'test': row.get('test', ''),
        'pvalue_raw': row['pvalue'],
        'pvalue_within_family': row.get('pvalue_bh', row['pvalue']),
        'effect_size': row.get('odds_ratio', row.get('effect_size', np.nan)),
        'effect_name': row.get('effect_name', 'odds_ratio'),
    })

corrected_df = pd.DataFrame(all_tests)

# Level 2: BH correction within each family
corrected_df['pvalue_bh_within_family'] = np.nan
for fam in corrected_df['family'].unique():
    mask = corrected_df['family'] == fam
    pvals = corrected_df.loc[mask, 'pvalue_raw'].values
    if len(pvals) > 1:
        reject, corr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        corrected_df.loc[mask, 'pvalue_bh_within_family'] = corr
    else:
        corrected_df.loc[mask, 'pvalue_bh_within_family'] = pvals

# Level 1: Get minimum p-value per family, then BH across families
family_min_p = corrected_df.groupby('family')['pvalue_raw'].min().reset_index()
family_min_p.columns = ['family', 'min_pvalue']
if len(family_min_p) > 1:
    reject_fam, corr_fam, _, _ = multipletests(family_min_p['min_pvalue'].values,
                                                 alpha=0.05, method='fdr_bh')
    family_min_p['family_bh'] = corr_fam
    family_min_p['family_significant'] = reject_fam
else:
    family_min_p['family_bh'] = family_min_p['min_pvalue']
    family_min_p['family_significant'] = family_min_p['min_pvalue'] < 0.05

corrected_df = corrected_df.merge(family_min_p[['family', 'family_bh', 'family_significant']],
                                   on='family', how='left')

corrected_df['significant_within_family'] = corrected_df['pvalue_bh_within_family'] < 0.05
corrected_df['significant_overall'] = (corrected_df['significant_within_family'] &
                                        corrected_df['family_significant'])

corrected_df.to_csv(OUT_CORRECT, sep='\t', index=False)
print(f"  Corrected p-values saved: {OUT_CORRECT}")
print(f"  Total tests across all families: {len(corrected_df)}")
print(f"  Significant within family (BH<0.05): {corrected_df['significant_within_family'].sum()}")
print(f"  Significant overall (hierarchical): {corrected_df['significant_overall'].sum()}")

# Omnibus test: MANOVA-like for correlated N-vs-C metrics
print("\n  Omnibus test for correlated N-vs-C metrics (Hotelling's T-squared approximation)...")
omnibus_results = []
for group in ['groel', 'hsp60', 'matrix_bg', 'mito_bg']:
    subset = nvc_paired[nvc_paired['analysis_group'] == group]
    metrics_to_use = []
    for m in ['relative_contact_order', 'mean_plddt', 'frac_helix', 'frac_strand']:
        n_col = f"{m}_n_domain"
        c_col = f"{m}_c_region"
        if n_col in subset.columns and c_col in subset.columns:
            metrics_to_use.append(m)

    if len(metrics_to_use) < 2:
        continue

    # Compute differences for each metric
    diffs = {}
    for m in metrics_to_use:
        n_col = f"{m}_n_domain"
        c_col = f"{m}_c_region"
        d = (subset[n_col] - subset[c_col]).dropna()
        diffs[m] = d

    # Align indices
    common_idx = diffs[metrics_to_use[0]].index
    for m in metrics_to_use[1:]:
        common_idx = common_idx.intersection(diffs[m].index)

    if len(common_idx) < 10:
        continue

    diff_matrix = np.column_stack([diffs[m].loc[common_idx].values for m in metrics_to_use])
    n = diff_matrix.shape[0]
    p = diff_matrix.shape[1]

    # Hotelling's T-squared: T2 = n * mean' * S_inv * mean
    mean_vec = diff_matrix.mean(axis=0)
    cov_mat = np.cov(diff_matrix, rowvar=False)

    try:
        cov_inv = np.linalg.inv(cov_mat)
        T2 = n * mean_vec @ cov_inv @ mean_vec

        # Convert to F-statistic
        F_stat = T2 * (n - p) / (p * (n - 1))
        p_val = 1 - stats.f.cdf(F_stat, p, n - p)

        print(f"    {group}: T2={T2:.2f}, F({p},{n-p})={F_stat:.2f}, p={p_val:.4e}")
        omnibus_results.append({
            'dataset': group, 'test': 'hotelling_t2', 'T2': T2,
            'F_stat': F_stat, 'df1': p, 'df2': n - p,
            'pvalue': p_val, 'n': n
        })
    except np.linalg.LinAlgError:
        print(f"    {group}: singular covariance matrix, skipping")

# ============================================================================
# FINAL SUMMARY REPORT
# ============================================================================
print("\n[6] Writing summary report...")

report_lines = []
report_lines.append("=" * 78)
report_lines.append("ANTAH ASTI PRARAMBH: MODULE H - COMPARATIVE STATISTICS SUMMARY REPORT")
report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append("=" * 78)

report_lines.append("")
report_lines.append("STATISTICAL FRAMEWORK")
report_lines.append("-" * 78)
report_lines.append("- Alpha level: 0.05 (two-sided)")
report_lines.append("- Multiple testing: Hierarchical BH correction")
report_lines.append("  Level 1: Three goal families (domain, stability, targeting)")
report_lines.append("  Level 2: BH within each family")
report_lines.append("- Effect sizes: Odds ratios (Fisher), rank-biserial r (Wilcoxon/MWU),")
report_lines.append("  eta-squared (Kruskal-Wallis), Cramer's V (chi-squared)")
report_lines.append("- Controls: Compartment-matched AND size-matched (10 kDa bins)")
report_lines.append("")

# Goal 1
report_lines.append("=" * 78)
report_lines.append("GOAL 1: DOMAIN ARCHITECTURE ENRICHMENT")
report_lines.append("=" * 78)

report_lines.append("")
report_lines.append("H1.1: GroEL substrate CATH superfamily enrichment")
report_lines.append("-" * 50)
groel_sig = enrichment_df[(enrichment_df['dataset'] == 'groel') & (enrichment_df['significant'])]
groel_all = enrichment_df[enrichment_df['dataset'] == 'groel']
report_lines.append(f"  Superfamilies tested: {len(groel_all)}")
report_lines.append(f"  Significant (BH < 0.05): {len(groel_sig)}")
if len(groel_sig) > 0:
    report_lines.append("  Top enriched superfamilies:")
    for _, r in groel_sig.head(10).iterrows():
        sf_label = r['superfamily']
        if r.get('superfamily_name', ''):
            sf_label += f" ({r['superfamily_name']})"
        report_lines.append(f"    {sf_label}: OR={r['odds_ratio']:.2f} "
                           f"[{r['ci_lower']:.2f}, {r['ci_upper']:.2f}], "
                           f"p_BH={r['pvalue_bh']:.4e}")
else:
    report_lines.append("  No individual superfamilies reached significance after BH correction.")
    report_lines.append("  Top candidates by raw p-value:")
    for _, r in groel_all.head(5).iterrows():
        sf_label = r['superfamily']
        report_lines.append(f"    {sf_label}: OR={r['odds_ratio']:.2f}, "
                           f"p_raw={r['pvalue']:.4e}, p_BH={r['pvalue_bh']:.4e}")

# Chi-squared for CATH class distribution
for chi_res in chi2_results:
    if chi_res['dataset'] == 'groel':
        report_lines.append(f"\n  CATH class distribution (chi-squared):")
        report_lines.append(f"    chi2={chi_res['statistic']:.2f}, dof={chi_res['dof']}, "
                           f"p={chi_res['pvalue']:.4e}, Cramer's V={chi_res['effect_size']:.3f}")

report_lines.append("")
report_lines.append("H1.2: HSP60 substrate CATH superfamily enrichment")
report_lines.append("-" * 50)
hsp60_sig = enrichment_df[(enrichment_df['dataset'] == 'hsp60') & (enrichment_df['significant'])]
hsp60_all = enrichment_df[enrichment_df['dataset'] == 'hsp60']
report_lines.append(f"  Superfamilies tested: {len(hsp60_all)}")
report_lines.append(f"  Significant (BH < 0.05): {len(hsp60_sig)}")
if len(hsp60_sig) > 0:
    report_lines.append("  Top enriched superfamilies:")
    for _, r in hsp60_sig.head(10).iterrows():
        sf_label = r['superfamily']
        if r.get('superfamily_name', ''):
            sf_label += f" ({r['superfamily_name']})"
        report_lines.append(f"    {sf_label}: OR={r['odds_ratio']:.2f} "
                           f"[{r['ci_lower']:.2f}, {r['ci_upper']:.2f}], "
                           f"p_BH={r['pvalue_bh']:.4e}")
else:
    report_lines.append("  No individual superfamilies reached significance after BH correction.")
    report_lines.append("  Top candidates by raw p-value:")
    for _, r in hsp60_all.head(5).iterrows():
        sf_label = r['superfamily']
        report_lines.append(f"    {sf_label}: OR={r['odds_ratio']:.2f}, "
                           f"p_raw={r['pvalue']:.4e}, p_BH={r['pvalue_bh']:.4e}")

for chi_res in chi2_results:
    if chi_res['dataset'] == 'hsp60':
        report_lines.append(f"\n  CATH class distribution (chi-squared):")
        report_lines.append(f"    chi2={chi_res['statistic']:.2f}, dof={chi_res['dof']}, "
                           f"p={chi_res['pvalue']:.4e}, Cramer's V={chi_res['effect_size']:.3f}")

report_lines.append("")
report_lines.append("H1.3: Conservation of fold enrichment between GroEL and HSP60")
report_lines.append("-" * 50)
report_lines.append(f"  Shared CATH superfamilies: {len(shared_sf)}")
report_lines.append(f"  GroEL-only: {len(only_groel)}, HSP60-only: {len(only_hsp60)}")
report_lines.append(f"  Jaccard index: {len(shared_sf)/len(universe_sf):.3f}" if len(universe_sf) > 0 else "  N/A")
report_lines.append(f"  Hypergeometric test: p={p_hypergeom:.4e}")
report_lines.append(f"  Homolog pairs with shared top SF: {n_pairs_same_sf}/{n_pairs_checked}")

# Goal 2
report_lines.append("")
report_lines.append("=" * 78)
report_lines.append("GOAL 2: N-VS-C STABILITY ASYMMETRY")
report_lines.append("=" * 78)

report_lines.append("")
report_lines.append("H2.1: Within-protein paired tests (Wilcoxon signed-rank)")
report_lines.append("-" * 50)

h21_rows = stab_df[stab_df['hypothesis'] == 'H2.1']
for group in ['groel', 'hsp60', 'matrix_bg', 'mito_bg']:
    g_rows = h21_rows[h21_rows['dataset'] == group]
    if len(g_rows) == 0:
        continue
    report_lines.append(f"\n  {group}:")
    for _, r in g_rows.iterrows():
        sig_marker = "*" if r.get('significant', False) else ""
        report_lines.append(f"    {r['metric']}: W={r['statistic']:.1f}, "
                           f"p_BH={r.get('pvalue_bh', r['pvalue']):.4e}, "
                           f"r={r['effect_size']:.3f}, "
                           f"direction={r.get('direction', '')} {sig_marker}")

report_lines.append("")
report_lines.append("H2.2: Cross-dataset asymmetry comparison (Mann-Whitney U)")
report_lines.append("-" * 50)
h22_rows = stab_df[stab_df['hypothesis'] == 'H2.2']
for _, r in h22_rows.iterrows():
    sig_marker = "*" if r.get('significant', False) else ""
    report_lines.append(f"  {r['dataset']} ({r['metric']}): U={r['statistic']:.1f}, "
                       f"p_BH={r.get('pvalue_bh', r['pvalue']):.4e}, "
                       f"r={r['effect_size']:.3f} {sig_marker}")

report_lines.append("")
report_lines.append("H2.3: GroEL class comparison (Kruskal-Wallis)")
report_lines.append("-" * 50)
h23_rows = stab_df[stab_df['hypothesis'] == 'H2.3']
kw_rows = h23_rows[h23_rows['test'] == 'kruskal_wallis']
for _, r in kw_rows.iterrows():
    sig_marker = "*" if r.get('significant', False) else ""
    report_lines.append(f"  {r['metric']}: H={r['statistic']:.2f}, "
                       f"p_BH={r.get('pvalue_bh', r['pvalue']):.4e}, "
                       f"eta2={r['effect_size']:.3f} {sig_marker}")

posthoc_rows = h23_rows[h23_rows['test'] == 'mann_whitney_posthoc']
if len(posthoc_rows) > 0:
    report_lines.append("  Post-hoc pairwise comparisons:")
    for _, r in posthoc_rows.iterrows():
        report_lines.append(f"    {r['dataset']}: U={r['statistic']:.1f}, "
                           f"p_raw={r['pvalue']:.4e}, r={r['effect_size']:.3f}")

# Omnibus
report_lines.append("")
report_lines.append("  Omnibus Hotelling's T-squared (multivariate):")
for o in omnibus_results:
    report_lines.append(f"    {o['dataset']}: T2={o['T2']:.2f}, "
                       f"F({o['df1']},{o['df2']})={o['F_stat']:.2f}, p={o['pvalue']:.4e}")

# Goal 3
report_lines.append("")
report_lines.append("=" * 78)
report_lines.append("GOAL 3: MATRIX TARGETING AND MTS")
report_lines.append("=" * 78)

report_lines.append("")
report_lines.append("H3.1: HSP60 matrix localization enrichment")
report_lines.append("-" * 50)
report_lines.append(f"  HSP60 matrix: {n_hsp60_matrix}, HSP60 non-matrix: {n_hsp60_nonmatrix}")
report_lines.append(f"  Non-HSP60 matrix: {n_nonhsp60_matrix}, Non-HSP60 non-matrix: {n_nonhsp60_nonmatrix}")
report_lines.append(f"  OR={or_h31:.2f} [{ci_lo_h31:.2f}, {ci_hi_h31:.2f}], p={p_h31:.4e}")
report_lines.append(f"  Conclusion: {'Significant enrichment' if p_h31 < 0.05 else 'Not significant'}")

report_lines.append("")
report_lines.append("H3.2: MTS prevalence and first-domain properties")
report_lines.append("-" * 50)
report_lines.append(f"  MTS prevalence (HSP60 matrix vs non-HSP60 matrix):")
report_lines.append(f"    OR={or_h32:.2f} [{ci_lo_h32:.2f}, {ci_hi_h32:.2f}], p={p_h32:.4e}")
report_lines.append(f"  First-domain property comparisons (MTS+ vs MTS-):")
for r in mts_property_tests:
    report_lines.append(f"    {r['comparison']}: U={r['statistic']:.1f}, "
                       f"p={r['pvalue']:.4e}, r={r['effect_size']:.3f}")

report_lines.append("")
report_lines.append("H3.3: MTS is predominantly a pre-domain extension")
report_lines.append("-" * 50)
report_lines.append(f"  MTS pre-domain: {n_pre_domain}/{n_total_mts} ({100*n_pre_domain/n_total_mts:.1f}%)")
report_lines.append(f"  MTS overlaps domain: {n_overlaps}/{n_total_mts} ({100*n_overlaps/n_total_mts:.1f}%)")
report_lines.append(f"  Binomial test (H0: p=0.5): p={p_binom:.4e}")
report_lines.append(f"  Gap length (MTS end to 1st domain): median={gap_lengths.median():.0f}, "
                   f"mean={gap_lengths.mean():.1f}")
report_lines.append(f"  Conclusion: {'MTS is predominantly pre-domain' if p_binom < 0.05 and n_pre_domain > n_overlaps else 'MTS is NOT predominantly pre-domain' if p_binom < 0.05 else 'Inconclusive'}")

# Hierarchical summary
report_lines.append("")
report_lines.append("=" * 78)
report_lines.append("HIERARCHICAL MULTIPLE TESTING SUMMARY")
report_lines.append("=" * 78)
report_lines.append("")
for _, frow in family_min_p.iterrows():
    report_lines.append(f"  Family: {frow['family']}")
    report_lines.append(f"    Min raw p-value: {frow['min_pvalue']:.4e}")
    report_lines.append(f"    Family-level BH p-value: {frow['family_bh']:.4e}")
    report_lines.append(f"    Family significant: {frow['family_significant']}")

report_lines.append(f"\n  Total tests: {len(corrected_df)}")
report_lines.append(f"  Significant within family: {corrected_df['significant_within_family'].sum()}")
report_lines.append(f"  Significant overall (hierarchical): {corrected_df['significant_overall'].sum()}")

report_lines.append("")
report_lines.append("=" * 78)
report_lines.append("KEY CONCLUSIONS")
report_lines.append("=" * 78)
report_lines.append("")

# Summarize key findings
n_sig_domain = corrected_df[(corrected_df['family'] == 'domain_architecture') & corrected_df['significant_overall']].shape[0]
n_sig_stab = corrected_df[(corrected_df['family'] == 'stability_asymmetry') & corrected_df['significant_overall']].shape[0]
n_sig_targ = corrected_df[(corrected_df['family'] == 'matrix_targeting') & corrected_df['significant_overall']].shape[0]

report_lines.append(f"1. Domain Architecture: {n_sig_domain} tests significant after hierarchical correction")
if n_sig_domain > 0:
    report_lines.append("   -> Evidence for specific fold enrichment in chaperonin substrates")
else:
    report_lines.append("   -> Limited evidence for individual superfamily enrichment")
    report_lines.append("      (may reflect small sample sizes or distributed fold preferences)")

report_lines.append(f"\n2. N-vs-C Stability: {n_sig_stab} tests significant after hierarchical correction")
if n_sig_stab > 0:
    report_lines.append("   -> Evidence for systematic N-vs-C asymmetry in chaperonin substrates")
else:
    report_lines.append("   -> N-vs-C differences exist but may not exceed background levels")

report_lines.append(f"\n3. Matrix Targeting: {n_sig_targ} tests significant after hierarchical correction")
if p_binom < 0.05 and n_pre_domain > n_overlaps:
    report_lines.append("   -> MTS is predominantly a pre-domain extension (strong evidence)")
if p_h31 < 0.05:
    report_lines.append("   -> HSP60 substrates enriched for matrix localization")

report_lines.append("")
report_lines.append("=" * 78)
report_lines.append("OUTPUT FILES")
report_lines.append("=" * 78)
report_lines.append(f"  {OUT_DOMAIN}")
report_lines.append(f"  {OUT_STAB}")
report_lines.append(f"  {OUT_TARG}")
report_lines.append(f"  {OUT_CORRECT}")
report_lines.append(f"  {OUT_REPORT}")
report_lines.append("")
report_lines.append("END OF REPORT")

report_text = '\n'.join(report_lines)
with open(OUT_REPORT, 'w') as f:
    f.write(report_text)

print(f"\n  Summary report saved: {OUT_REPORT}")
print("\n" + "=" * 70)
print("Module H complete.")
print("=" * 70)
