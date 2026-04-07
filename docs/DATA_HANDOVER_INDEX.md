# Antah Asti Prarambh — Data Handover Index

**Prepared for:** Collaborator
**Prepared by:** Vishal Bharti, CSIR-IGIB
**Date:** 2026-04-01 (updated)
**Version:** Phase 2 COMPLETE (including FoldX thermodynamic stability — 25,007 proteins)

---

## Quick Start

All results are in `results/phase2/`. The most important files to look at first:

1. **`stats/corrected_pvalues_full.tsv`** — All 62 statistical tests with p-values, effect sizes, and significance flags
2. **`figures/fig1–fig6.pdf`** — 6 publication-quality figures (use PDFs for vector quality)
3. **`stability/n_vs_c_paired_full.tsv`** — The core paired comparison dataset (2,648 multi-domain proteins)
4. **`stats/statistics_summary_full.txt`** — Human-readable summary of all results

---

## File-by-File Guide

### Figures (`results/phase2/figures/`)

| File | What it shows | Key takeaway |
|------|--------------|--------------|
| `fig1_domain_architecture.pdf` | CATH class distribution, top enriched superfamilies per dataset, domain count histogram | GroEL enriched for TIM barrels (OR=22.6) |
| `fig2_n_vs_c_stability.pdf` | Split violin plots of contact order and pLDDT for N-domain vs C-region across all datasets | N > C contact order is universal (p = 7.1×10⁻¹⁸ in mito background) |
| `fig3_groel_class_comparison.pdf` | N–C contact order difference boxplots for GroEL Class I, II, III | No class effect (KW p = 0.77) |
| `fig4_mts_targeting.pdf` | Sub-mitochondrial localization bars, MTS-domain gap histogram, cleavage vs domain start scatter | 84.4% MTS are pre-domain, median 12-residue gap |
| `fig5_orthology.pdf` | Orthogroup category breakdown, N-domain RCO scatter for homolog pairs | r = 0.84, p = 5.3×10⁻¹³ for RCO conservation |
| `fig6_summary.pdf` | Key findings overview: N–C asymmetry significance bars, test result heatmap, MTS pie chart | Visual summary of all three goals |

PNG versions (300 DPI) are provided alongside PDFs for presentations.

### Statistics (`results/phase2/stats/`)

| File | Format | Rows | How to use |
|------|--------|:----:|-----------|
| `corrected_pvalues_full.tsv` | TSV | 62 | **Master results table.** Each row is one statistical test. Key columns: `family` (H1/H2/H3), `hypothesis`, `test`, `statistic`, `p_value`, `effect_size`, `ci_low`/`ci_high`, `n1`/`n2`, `direction`, `p_bh_within`, `significant_overall`. Filter on `significant_overall == True` for the 45 significant results. |
| `statistics_summary_full.txt` | Text | 534 lines | Human-readable version of all tests grouped by family. `***` marks significant results. Good for reading through results without opening a spreadsheet. |
| `stability_comparisons_full.tsv` | TSV | 2,648 | Per-protein paired comparison data. Same as `n_vs_c_paired_full.tsv` but formatted for statistical analysis. |

### Domain Analysis (`results/phase2/domains/`)

| File | Rows | Key columns | How to use |
|------|:----:|-------------|-----------|
| `unified_domain_assignments_full.tsv` | 25,019 | `accession`, `source` (CATH/Chainsaw), `n_domains`, `top_superfamily`, `cath_class`, `datasets` | **Master domain table.** Every protein with its domain assignment source (75.3% CATH, 24.7% Chainsaw), domain count, and which datasets it belongs to. The `datasets` column contains comma-separated labels (groel, hsp60, matrix_bg, mito_bg, proteome_bg). |
| `chainsaw_full_predictions.tsv` | 25,007 | `accession`, `ndom`, `chopping`, `confidence` | Raw Chainsaw ML predictions. The `chopping` column gives domain boundaries in PDB-style notation. |
| `domain_distribution_full.tsv` | 56 | `dataset`, `n_domains`, `count`, `percent` | Domain count distribution per dataset. Each row = one (dataset, domain_count) combination. Ready for bar chart plotting. |
| `chainsaw_full_predictions_annotated.tsv` | 25,007 | Same + superfamily annotations | Chainsaw predictions with CATH superfamily labels added. |

### N-vs-C Structural Analysis (`results/phase2/stability/`)

| File | Rows | Key columns | How to use |
|------|:----:|-------------|-----------|
| `n_vs_c_paired_full.tsv` | 2,648 | `accession`, `source`, `n_domains`, `datasets`, `groel_class`, `*_n_domain`, `*_c_region`, `*_diff` (for each metric) | **Core analysis file.** Each row is a multi-domain protein with paired N-domain vs C-region values for: length, net charge, hydrophobicity, helix/strand fraction, pLDDT, contact order. The `*_diff` columns = N − C. |
| `region_boundaries_full.tsv` | 5,322 | `accession`, `n_domain_start`, `n_domain_end`, `c_region_start`, `c_region_end`, `pre_domain_length` | Three-region model boundaries. Includes single-domain proteins (only N-domain, no C-region). |
| `contact_order_full.tsv` | 11,824 | `accession`, `domain_id`, `region_type`, `relative_contact_order`, `n_contacts`, `n_residues` | Per-domain/region contact order. `region_type` = n_domain or c_region. For deep-dive analysis. |
| `structure_metrics_full.tsv` | 11,824 | `accession`, `domain_id`, `mean_plddt`, `frac_helix`, `frac_strand`, `frac_coil` | Per-domain structural metrics from STRIDE secondary structure assignment. |
| `sequence_metrics_full.tsv` | 11,824 | `accession`, `domain_id`, `frac_hydrophobic`, `net_charge`, `length` | Per-domain sequence composition metrics. |

### Structural Clustering (`results/phase2/foldseek/analysis/`)

| File | Rows | How to use |
|------|:----:|-----------|
| `foldseek_clusters_full.tsv` | 16,242 | Cluster-level summary. Each row is one cluster with size, representative protein, and dataset labels. |
| `combined_cluster_membership.tsv` | 27,063 | Per-protein cluster assignment. Join on `accession` with other tables to analyze cluster membership. |
| `foldseek_full_summary.txt` | — | Human-readable clustering statistics: total counts, size distribution, top 30 largest clusters. |

### FoldX Stability (`results/phase2/foldx/`) — COMPLETE

| File | Rows | Description |
|------|-----:|-------------|
| `foldx_stability_all.tsv` | 25,007 | Aggregated FoldX total energy for all proteins. Columns: accession, status, total_energy (kcal/mol). All 25,007 proteins succeeded. Mean=309.6, Median=111.2, Range=[-242.7, 11442.8]. GroEL substrates significantly lower (median -38.6 vs -15.2 bg, p=2.9e-3 compartment-matched). |
| `foldx_array.slurm` | — | SLURM array job script (501 tasks, 50 proteins each). |
| `collect_results.slurm` | — | Collection job script. |

**Note:** Per-protein JSON intermediates (`per_protein/`) are retained on HPC only. The aggregated TSV contains all needed data.

---

## Dataset Labels

Throughout the result files, proteins are tagged with dataset membership. Here is the mapping:

| Label | Dataset | Size | Description |
|-------|---------|:----:|-------------|
| `groel` | Dataset 4 | 252 | *E. coli* GroEL substrates (Kerner et al. 2005) |
| `hsp60` | Dataset 5 | 266 | Human HSP60 Tier 1 substrates |
| `matrix_bg` | Dataset 7 | 525 | Mitochondrial matrix proteins (MitoCarta 3.0) |
| `mito_bg` | Dataset 3 | 1,136 | Full mitochondrial proteome |
| `proteome_bg` | Datasets 1+2 | ~24,000 | Combined *E. coli* + human proteome (background) |

A protein can belong to multiple datasets (e.g., an HSP60 substrate may also be in `matrix_bg` and `mito_bg`).

---

## Column Name Conventions

- `*_n_domain` / `*_c_region` — Metric value for the N-terminal domain / C-terminal region
- `*_diff` — Difference (N − C), positive = N higher
- `frac_*` — Fraction (0–1 scale)
- `mean_plddt` — AlphaFold confidence score (0–100), NOT thermodynamic stability
- `relative_contact_order` — Plaxco et al. RCO (higher = more complex topology = slower folding)
- `p_bh_within` — Benjamini-Hochberg corrected p-value within hypothesis family
- `significant_overall` — True/False after hierarchical correction at both levels
- `groel_class` — GroEL dependence class: I (spontaneous), II (partial), III (obligate)

---

## Reproducing the Analysis

All scripts are in `workflow/phase2/`:
- `run_foldx.py` — FoldX orchestrator
- `slurm_jobs/09_module_f_stability.sh` — N-vs-C analysis (embedded Python)
- `slurm_jobs/10_module_h_stats.sh` — Statistical tests (embedded Python)
- `slurm_jobs/11_module_i.sh` — Figure generation (embedded Python)
- `scripts/module_i_polished.py` — Local figure polishing script
- `config.yaml` — All pipeline parameters

---

## Known Limitations

1. **pLDDT ≠ stability.** AlphaFold pLDDT is model confidence, not thermodynamic stability. FoldX total energy (25,007 proteins) addresses this: GroEL substrates have median -38.6 vs -15.2 for E. coli background (p=2.9e-3, d=-0.07, compartment-matched). Note: FoldX was parameterized on experimental structures, not AlphaFold models — caveat for publication.
2. **Chainsaw domain predictions** are ML-based; where CATH assignments exist, they are preferred.
3. **MTS analysis** uses CATH domain boundaries (preferred over Chainsaw in the unified assignments).
4. **8 proteins** lack AlphaFold structures: P07203, P30042, P36969, Q16881, Q5THJ4, Q86UA3, Q9BVL4, Q9NNW7.

---

*All data derived from actual computational results. No simulated or fabricated data.*
