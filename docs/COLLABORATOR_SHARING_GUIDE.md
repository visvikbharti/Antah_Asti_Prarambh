# Antah Asti Prarambh — Complete Sharing Package for Collaborator

**Prepared by:** Vishal Bharti, CSIR-IGIB
**Date:** 2026-03-22

This document lists everything to share with your collaborator for **full reproducibility** — raw data, processed datasets, all analysis scripts + HPC pipeline, results, figures, and documentation.

---

## What to Share (Summary)

| Category | What | Size |
|----------|------|------|
| Raw data | Literature source files, UniProt proteomes, MitoCarta, AlphaFold structures | ~510 MB |
| Processed data | 7 cleaned datasets + FASTA files | ~1.7 MB |
| Analysis scripts | 20 Python analysis scripts | ~780 KB |
| HPC pipeline scripts | Snakefile, config, 7 Python scripts, 19 SLURM jobs | ~930 KB |
| Core analysis results | Domain, homology, stability, MTS, statistics, figures | ~90 MB |
| Full-scale results | Full-scale domains, stability, stats, figures | ~19 MB |
| Documentation | 14 markdown files | ~440 KB |
| **Total** | | **~620 MB** |

**Note:** The AlphaFold structures (465 MB) are the bulk of the size. If bandwidth is a concern, share everything except `data/raw/alphafold/pilot/` (the locally-cached structure set) and tell the collaborator to re-download using the provided script.

---

## Complete File List — Organized by Category

### A. RAW DATA (`data/raw/`)

These are the unmodified source files. The collaborator needs these to verify data provenance.

#### A1. Literature source files (root directory)
```
kerner_2005_groel_interactors_table_s3.csv    — Original GroEL substrate table (Kerner et al. 2005)
kerner_2005_groel_interactors_clean.csv       — Cleaned GroEL table
12192_2020_1080_MOESM4_ESM.xlsx               — HSP60 interactome supplementary (Bruderer et al. 2020)
hsp60_interactome_clean.tsv                   — Cleaned HSP60 interactome
```

#### A2. UniProt proteomes (`data/raw/uniprot/`)
```
ecoli_k12_proteome.fasta       — E. coli K-12 reference proteome (4,403 proteins)
ecoli_k12_proteome.tsv         — E. coli K-12 metadata (accession, gene, length, etc.)
human_proteome.fasta           — Human reference proteome (20,416 proteins)
human_proteome.tsv             — Human proteome metadata
```

#### A3. MitoCarta (`data/raw/mitocarta/`)
```
Human.MitoCarta3.0.xls         — MitoCarta 3.0 (9.7 MB) — mitochondrial proteome annotations
```

#### A4. AlphaFold structures (`data/raw/alphafold/pilot/`)
```
1,382 CIF files                — AF-{ACC}-F1-model_v4.cif format (465 MB total)
                               — AlphaFold DB v6 predictions for all 7 datasets
```
*Note: This directory holds the locally-cached subset. The full structure set (25,007) is on HPC only (~15 GB) and not transferred locally.*

#### A5. Duplicate copies in `data/raw/custom/`
```
Same 4 literature files as root — kept for pipeline path references
```

---

### B. PROCESSED DATA (`data/processed/`)

These are the 7 cleaned, standardized datasets ready for analysis. Critical for reproducibility.

```
groel_substrates_standardized.tsv    — Dataset 4: 252 GroEL substrates (demerged accessions)
groel_substrates.fasta               — Dataset 4: protein sequences
hsp60_tier1_substrates.tsv           — Dataset 5: 266 HSP60 Tier 1 substrates (SILAC-filtered)
hsp60_tier1_substrates.fasta         — Dataset 5: protein sequences
hsp60_interactome_standardized.tsv   — Full HSP60 interactome (315 proteins, all tiers)
human_mito_proteome.tsv              — Dataset 3: 1,136 mitochondrial proteins (MitoCarta 3.0)
human_matrix_proteome.tsv            — Dataset 7: 525 matrix-only proteins
groel_hsp60_homologs.tsv             — Dataset 6: 69 cross-species homolog pairs
hsp60_filtering_report.txt           — QC report: exclusions, tier assignment rationale
mitocarta_summary_report.txt         — MitoCarta v2→v3 reclassification report
```

---

### C. ANALYSIS SCRIPTS

#### C1. Data preparation scripts (`scripts/`)
```
validate_uniprot_accessions.py       — UniProt REST API validation (demerge resolution)
filter_hsp60_interactome.py          — SILAC-based HSP60 filtering pipeline
module_c_extract_fasta.py            — FASTA extraction for MMseqs2 input
module_c_analyze_rbh.py              — Reciprocal best hit analysis
```

#### C2. Analysis pipeline scripts (`workflow/scripts/`)
```
download_alphafold_pilot.py          — Module A: AlphaFold structure download
run_dssp.py                          — Module B: DSSP secondary structure
get_cath_domains.py                  — Module D1: CATH domain assignment (Gene3D API)
run_chainsaw_e2.py                   — Module D2/E2: Chainsaw ML domain prediction
validate_structure_quality.py        — Module D4: Structure quality validation
compute_domain_structural_metrics.py — Module E: Domain-level metrics
analyze_foldseek.py                  — Module E3: Foldseek cluster analysis
run_orthology.py                     — Module C: OrthoFinder + MMseqs2 orthology
build_dataset6_homologs.py           — Module C: Dataset 6 construction
domain_distribution_summary.py      — Module E: Domain class distribution
parse_mitocarta.py                   — Module G: MitoCarta parsing
module_f_n_vs_c_analysis.py          — Module F: N-vs-C terminus stability
module_f_extension_chainsaw.py       — Module F: Chainsaw-based extension
module_g_mts_analysis.py             — Module G: MTS targeting analysis
module_h_comparative_stats.py        — Module H: Statistical tests
generate_figures.py                  — Module I: Publication figures
```

---

### D. HPC PIPELINE SCRIPTS

#### D1. Pipeline orchestration (`workflow/phase2/`)
```
Snakefile                            — Snakemake workflow (717 lines, all rules)
config.yaml                          — Pipeline configuration (paths, params, resources)
README.md                            — HPC pipeline documentation
download_alphafold_full.py           — Bulk AlphaFold download (E. coli + Human)
run_foldseek_full.py                 — Foldseek search + clustering orchestrator
run_foldx.py                         — FoldX orchestrator (generate/chunk/collect modes)
```

#### D2. Full-scale analysis scripts (`workflow/phase2/scripts/`)
```
module_f_full.py                     — Module F: N-vs-C on 25,007 proteins
module_h_full.py                     — Module H: 62 statistical tests, hierarchical BH
module_i_full.py                     — Module I: 6 publication figures (HPC version)
module_i_polished.py                 — Module I: Polished figures (local version)
```

#### D3. SLURM job scripts (`workflow/phase2/slurm_jobs/`)

These are the actual HPC submission scripts. Each is self-contained with embedded Python where needed.

```
00_setup.sh                          — Environment setup (conda, paths, directories)
01_download_alphafold.sh             — AlphaFold bulk download (4,371 ecoli + 20,636 human)
02_foldseek_createdb.sh              — Foldseek database creation
03_foldseek_search.sh                — Foldseek all-vs-all structural search
04_foldseek_cluster.sh               — Foldseek clustering (16,242 clusters)
05_chainsaw.sh                       — Chainsaw domain prediction (25,007 proteins)
05a_decompress_human.sh              — Decompress .cif.gz → .cif for human structures
06_foldx_generate.sh                 — Generate FoldX array job script
07_foldx_collect.sh                  — Collect FoldX per-protein results into TSV
08_module_e_domains.sh               — Module E: Unified domain integration
09_module_f.sh                       — Module F: N-vs-C analysis (variant 1)
09_module_f_stability.sh             — Module F: N-vs-C analysis (variant 2, embedded Python)
10_module_h.sh                       — Module H: Statistics (variant 1)
10_module_h_stats.sh                 — Module H: Statistics (variant 2, embedded Python)
11_module_i.sh                       — Module I: Figures (variant 1)
11_module_i_figures.sh               — Module I: Figures (variant 2, embedded Python)
submit_pipeline.sh                   — Sequential pipeline submission wrapper
submit_analysis.sh                   — Analysis chain submission wrapper
test_stride_fix.sh                   — STRIDE binary compatibility test
```

**Important for collaborator:** The `a` variants (09_module_f_stability.sh, 10_module_h_stats.sh, 11_module_i_figures.sh) contain embedded Python scripts and are the ones that were actually run for final results.

#### D4. Environment (`workflow/envs/`)
```
(check for any conda environment YAML files)
```

---

### E. CORE ANALYSIS RESULTS

#### E1. Structures (`results/structures/`)
```
structure_index.tsv                  — 1,382 structure metadata (accession, pLDDT, length)
dssp_summary.tsv                     — Per-protein DSSP (helix%, strand%, coil%)
dssp_per_residue.tsv                 — Per-residue secondary structure (7.3 MB)
dssp/                                — 1,382 individual DSSP files
structure_quality_validation.tsv     — Quality scores and flags
flagged_low_quality_structures.tsv   — 63 flagged structures
quality_per_dataset.tsv              — Quality breakdown by dataset
quality_validation_report.txt        — QC narrative
```

#### E2. Domains (`results/domains/`)
```
cath_domain_assignments.tsv          — CATH domain assignments (Gene3D)
cath_protein_summary.tsv             — Per-protein CATH summary
chainsaw_domain_predictions.tsv      — Chainsaw ML predictions
chainsaw_raw_output.tsv              — Chainsaw raw output
ml_domain_assignments.tsv            — Unified CATH + Chainsaw (1,387 proteins)
domain_structural_metrics.tsv        — Per-domain contact order, pLDDT
domain_distribution_summary.tsv      — Domain count distributions
foldseek_clusters.tsv                — 1,155 structural clusters
foldseek_cluster_summary.txt         — Cluster analysis report
```

#### E3. Homology (`results/homology/`)
```
rbh_groel_hsp60.tsv                  — 40 RBH pairs (MMseqs2)
rbh_groel_hsp60_annotated.tsv        — Annotated RBH pairs
orthogroups_ecoli_human.tsv          — 422 orthogroups
substrate_orthogroups.tsv            — 34 shared substrate orthogroups
orthology_comparison.tsv             — RBH vs orthogroup comparison
rbh_summary_report.txt               — RBH analysis report
orthology_summary_report.txt         — Orthology report
```

#### E4. N-vs-C Terminus (`results/termini/`)
```
region_boundaries.tsv                — Three-region model boundaries
region_boundaries_extended.tsv       — Extended with Chainsaw proteins
n_vs_c_paired.tsv                    — Paired N-vs-C metrics
n_vs_c_paired_extended.tsv           — Extended paired comparisons
contact_order.tsv                    — Per-domain contact order
contact_order_extended.tsv           — Extended contact order
sequence_metrics.tsv                 — Sequence composition per region
structure_metrics.tsv                — Structural metrics per region
```

#### E5. MTS Targeting (`results/mts/`)
```
combined_targeting.tsv               — Integrated targeting predictions
mts_domain_relationship.tsv          — MTS vs domain boundary analysis
uniprot_transit_signal_cache.tsv     — UniProt transit peptide annotations
targeting_summary_report.txt         — MTS summary
```

#### E6. Statistics (`results/stats/`)
```
corrected_pvalues.tsv                — Statistical tests with hierarchical BH correction
statistics_summary_report.txt        — Human-readable summary
stability_comparisons.tsv            — Stability comparison results
domain_enrichment.tsv                — Superfamily enrichment (OR, CI)
targeting_stats.tsv                  — MTS statistics
```

#### E7. Figures (`results/figures/`)
```
fig1_domain_architecture.{pdf,png}   — CATH class distribution, superfamilies
fig2_n_vs_c_stability.{pdf,png}      — Contact order split violins
fig3_groel_class_comparison.{pdf,png} — GroEL class I/II/III comparison
fig4_mts_targeting.{pdf,png}         — MTS-domain spatial relationship
fig5_orthology.{pdf,png}             — Cross-species RCO conservation
fig6_summary.{pdf,png}               — Key findings summary
```

---

### F. FULL-SCALE RESULTS

#### F1. Domains (`results/phase2/domains/`)
```
chainsaw_full_predictions.tsv            — 25,007 proteins, Chainsaw ML predictions
chainsaw_full_predictions_annotated.tsv  — Annotated with superfamilies
unified_domain_assignments_full.tsv      — 25,258 unified (CATH + Chainsaw)
domain_distribution_full.tsv             — Domain count by dataset
```

#### F2. Structural Clustering (`results/phase2/foldseek/analysis/`)
```
foldseek_clusters_full.tsv               — 16,242 clusters
combined_cluster_membership.tsv          — 27,063 per-protein assignments
foldseek_full_summary.txt                — Clustering summary
```

#### F3. N-vs-C Stability (`results/phase2/stability/`)
```
n_vs_c_paired_full.tsv                   — 2,648 paired comparisons (30 columns)
region_boundaries_full.tsv               — 5,322 region boundaries
contact_order_full.tsv                   — 11,824 contact order values
structure_metrics_full.tsv               — 11,824 pLDDT + SS per domain
sequence_metrics_full.tsv                — 11,824 composition per domain
```

#### F4. Statistics (`results/phase2/stats/`)
```
corrected_pvalues_full.tsv               — 62 tests, hierarchical BH correction
statistics_summary_full.txt              — Human-readable summary (534 lines)
stability_comparisons_full.tsv           — Detailed paired comparison results
```

#### F5. Figures (`results/phase2/figures/`)
```
fig1_domain_architecture.{pdf,png}       — Updated for full-scale data
fig2_n_vs_c_stability.{pdf,png}          — Updated for full-scale data
fig3_groel_class_comparison.{pdf,png}    — Updated for full-scale data
fig4_mts_targeting.{pdf,png}             — Updated for full-scale data
fig5_orthology.{pdf,png}                 — Updated for full-scale data
fig6_summary.{pdf,png}                   — Updated for full-scale data
```

#### F6. Structures metadata (`results/phase2/structures/`)
```
structure_index_full.tsv                 — 25,007 structure metadata
missing_accessions_ecoli.txt             — 8 missing E. coli structures
missing_accessions_human.txt             — 10 missing human structures
```

#### F7. FoldX — COMPLETE (April 1, 2026)
```
results/phase2/foldx/foldx_stability_all.tsv — Aggregated stability (25,007 proteins, 1.1 MB)
  Columns: accession, status, total_energy (kcal/mol)
  Key: GroEL median=-38.6, HSP60 median=74.6, proteome median=119.2
```

---

### G. DOCUMENTATION (`docs/`)

```
PROJECT_PLAN.md                      — Project strategy, goals, timeline
PRIMARY_HYPOTHESES.md                — Pre-registered hypotheses (H1, H2, H3)
METHODS_AND_PROTOCOLS.md             — Detailed computational protocols
DOCUMENTATION.md                     — General project documentation (114 KB)
RESULTS_NARRATIVE.md                 — Scientific interpretation of findings (32 KB)
PHASE1_VERIFICATION.md              — Analysis completion checklist
PHASE2_RESULTS_REPORT.md            — Full-scale results for collaborator (16 KB)
HPC_DEPLOYMENT_GUIDE.md             — HPC cluster setup + job submission
SESSION_CONTINUITY.md               — Cross-session state management
SESSION4_DOCUMENTATION.md           — Session 4: bug fixes
SESSION5_DOCUMENTATION.md           — Session 5: analysis chain + FoldX
COLLABORATOR_PRESENTATION.md        — Meeting presentation document
DATA_HANDOVER_INDEX.md              — File-by-file guide to results
COLLABORATOR_SHARING_GUIDE.md       — This document
TODO.md                              — Outstanding tasks
```

---

## How to Share

### Option 1: Full project (recommended)

Tar the entire project directory excluding AlphaFold structures and DSSP individual files:

```bash
cd ~/Downloads/
tar czf Antah_Asti_Prarambh_share.tar.gz \
    --exclude='data/raw/alphafold' \
    --exclude='results/structures/dssp' \
    --exclude='results/phase2/foldseek/clusters' \
    --exclude='results/phase2/foldseek/databases' \
    --exclude='results/phase2/foldseek/search' \
    --exclude='results/phase2/foldseek/combined_structures' \
    --exclude='results/homology/_mmseqs_ortho_work' \
    --exclude='.claude' \
    --exclude='.DS_Store' \
    --exclude='foldx5_Linux.zip' \
    Antah_Asti_Prarambh/
```

**Estimated size: ~30-40 MB** (without structures)

The collaborator can re-download AlphaFold structures using the provided scripts. Include a note about this.

### Option 2: Essential files only (~15 MB)

If you want a minimal package:

```bash
mkdir -p Antah_Asti_Prarambh_essential/
# Copy essential directories
cp -r data/processed/ Antah_Asti_Prarambh_essential/data_processed/
cp -r results/phase2/ Antah_Asti_Prarambh_essential/results_phase2/
cp -r results/stats/ Antah_Asti_Prarambh_essential/results_stats/
cp -r results/figures/ Antah_Asti_Prarambh_essential/figures/
cp -r docs/ Antah_Asti_Prarambh_essential/docs/
cp -r workflow/ Antah_Asti_Prarambh_essential/workflow/
cp -r scripts/ Antah_Asti_Prarambh_essential/scripts/
# Copy root data files
cp *.csv *.xlsx *.tsv Antah_Asti_Prarambh_essential/ 2>/dev/null
```

### Option 3: HPC scripts separately

If the collaborator wants to run the HPC pipeline on their own cluster:

```bash
tar czf Antah_Asti_Prarambh_HPC_pipeline.tar.gz \
    workflow/phase2/ \
    docs/HPC_DEPLOYMENT_GUIDE.md \
    data/processed/
```

They will need: Python 3.9+, Foldseek, Chainsaw, STRIDE (heiniglab), FoldX 5.1, gemmi, pandas, scipy, matplotlib, seaborn.

---

## Collaborator Quick-Start

Tell your collaborator to read these files in this order:

1. **`docs/COLLABORATOR_PRESENTATION.md`** — Full overview, findings, figures
2. **`docs/DATA_HANDOVER_INDEX.md`** — What each result file contains
3. **`docs/PRIMARY_HYPOTHESES.md`** — The pre-registered hypotheses being tested
4. **`docs/METHODS_AND_PROTOCOLS.md`** — Detailed computational methods
5. **`results/phase2/stats/statistics_summary_full.txt`** — All 62 tests at a glance
6. **`results/phase2/figures/fig1-6.pdf`** — The publication figures

For pipeline reproduction:
1. **`docs/HPC_DEPLOYMENT_GUIDE.md`** — HPC setup
2. **`workflow/phase2/README.md`** — Pipeline documentation
3. **`workflow/phase2/config.yaml`** — All parameters

---

## What is NOT Included (HPC-only files)

These large files exist only on the HPC and are not shared:

| File | Size | Why not shared |
|------|------|---------------|
| AlphaFold structures (full) | ~15 GB | Re-downloadable via script |
| Foldseek databases | ~5 GB | Intermediate, regenerable |
| FoldX work directories | ~10 GB | Intermediate, regenerable |
| SLURM log files | ~500 MB | Session-specific, not needed |
| Chainsaw batch intermediates | ~2 GB | Intermediate, already merged |

All of these can be regenerated by running the pipeline from the provided scripts.
