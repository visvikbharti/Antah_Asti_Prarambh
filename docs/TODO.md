# Antah Asti Prarambh - Task Tracker

## Legend
- [ ] Not started
- [~] In progress
- [x] Complete

---

## Phase 0: Setup & Data Cleaning

### Module A: Data Acquisition & Cleaning
- [x] A1: Set up project directory structure
- [x] A2: Copy raw data files into project structure
- [x] A3: GroEL substrates — UniProt ID remapping (252 proteins)
  - 252/252 resolved; 149 demerged to K-12 entries; 0 obsolete
  - 181 cytoplasmic, 71 non-cytoplasmic (42 unknown)
  - 211 have SCOP folds parsed from Table S3
  - 4 accessions not in K-12 reference proteome (plasmid/strain-specific)
  - Output: `data/processed/groel_substrates_standardized.tsv`
- [x] A4: HSP60 interactome — SILAC-based filtering
  - 266 Tier 1 (MitoCarta+ & SILAC>5), 49 Tier 2, 10 excluded
  - NDIC values imputed at 2x 95th percentile
  - Output: `data/processed/hsp60_tier1_substrates.tsv`
- [x] A5: MitoCarta 3.0 — download and process
  - 1,136 mito proteins; 525 matrix-localized
  - 70 localization changes vs MC2 documented
  - Output: `data/processed/human_mito_proteome.tsv`, `human_matrix_proteome.tsv`

### Module B: Reference Proteome Acquisition
- [x] B1: Download E. coli K-12 proteome from UniProt (FASTA + TSV)
  - 4,403 reviewed Swiss-Prot entries
  - Output: `data/raw/uniprot/ecoli_k12_proteome.fasta`, `.tsv`
- [x] B2: Download Human proteome from UniProt (FASTA + TSV)
  - 20,416 reviewed Swiss-Prot entries
  - Output: `data/raw/uniprot/human_proteome.fasta`, `.tsv`
- [x] B3: Install MMseqs2 (v18.8cc5c) and Foldseek (v10.941cd33)
  - Installed in conda env `proteomics`
  - Path: `/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/`
- [x] B4: Cross-check GroEL substrates against E. coli proteome
  - 248/252 found; 4 are plasmid/strain-specific (expected)
- [x] B5: Create project plan document
  - Output: `docs/PROJECT_PLAN.md` (1,097 lines)

---

## Phase 1: Pilot Analysis (~900 proteins, all LOCAL)

### Module C: Orthology / Homology Layer
- [x] C1: Extract FASTA sequences for GroEL substrates (252 proteins)
  - 248 from K-12 proteome + 4 fetched from UniProt API
  - Output: `data/processed/groel_substrates.fasta`
- [x] C2: Extract FASTA sequences for HSP60 Tier 1 substrates (266 proteins)
  - 264 from human proteome + 2 fetched from UniProt API (A0A087WU62, G3V325)
  - Output: `data/processed/hsp60_tier1_substrates.fasta`
- [x] C3: Run MMseqs2 reciprocal best hit (GroEL vs HSP60)
  - **40 RBH pairs found** (15.9% of GroEL substrates)
  - By class: I=7/38 (18.4%), II=24/126 (19.0%), III=8/84 (9.5%)
  - Median identity: 35.8%, median qcov: 0.924
  - Top pairs: sucA/OGDH (652 bits), sdhA/SDHA (591), fusA/GFM1 (572)
  - Output: `results/homology/rbh_groel_hsp60.tsv`, `rbh_groel_hsp60_annotated.tsv`
  - Report: `results/homology/rbh_summary_report.txt`
- [x] C4: Run orthology analysis (MMseqs2 all-vs-all + connected-component clustering)
  - OrthoFinder not available; used MMseqs2 bidirectional search as fallback
  - **422 orthogroups** (E. coli ↔ Human), 2,895 reciprocal pairs
  - **34 shared orthogroups** with both GroEL + HSP60 substrates
  - 51 GroEL-only orthogroups, 58 HSP60-only orthogroups
  - **62 substrate pairs** via orthogroups (vs 40 from simple RBH — 29 additional)
  - Output: `results/homology/orthogroups_ecoli_human.tsv`, `substrate_orthogroups.tsv`
  - Report: `results/homology/orthology_summary_report.txt`
  - Script: `workflow/scripts/run_orthology.py`
- [x] C5: Compare RBH vs orthogroup results
  - 33/40 RBH pairs confirmed by orthogroups; 7 lost to stricter coverage filtering
  - 29 additional pairs found via many-to-many orthology (paralog families)
  - Output: `results/homology/orthology_comparison.tsv`
- [x] C6: Build Dataset 6 (2-way homologs)
  - Merged RBH + orthogroup evidence: **69 total pairs**
  - 33 found by both methods, 7 RBH-only, 29 orthogroup-only
  - By class: II=38, III=15, I=13, ambiguous=3
  - 48 unique GroEL, 56 unique HSP60 proteins paired
  - Output: `data/processed/groel_hsp60_homologs.tsv`
  - Script: `workflow/scripts/build_dataset6_homologs.py`

### Module D: Structure Acquisition & Indexing
- [x] D1: Download AlphaFold structures for pilot proteins
  - 1,390 unique proteins across all datasets
  - **1,382 downloaded** (v6), 8 failed (no AlphaFold model)
  - Failed: P07203, P30042, P36969, Q16881, Q5THJ4, Q86UA3, Q9BVL4, Q9NNW7
  - Disk: 466 MB in `data/raw/alphafold/pilot/`
  - Script: `workflow/scripts/download_alphafold_pilot.py`
- [x] D2: Build structure index
  - 1,382 structures with pLDDT stats
  - Mean pLDDT: 85.8, Median: 87.6, fraction >70: 0.831
  - 100% coverage for GroEL and HSP60; 99.3% for mito
  - Output: `results/structures/structure_index.tsv`
- [x] D3: Run DSSP on all pilot structures (secondary structure assignment)
  - 1,382/1,382 processed, zero failures. Runtime: 203 seconds.
  - Mean SS: 43.5% helix, 14.2% strand, 42.2% coil (helix-rich, consistent with mito proteins)
  - Output: `results/structures/dssp_summary.tsv` (per-protein), `dssp_per_residue.tsv` (518K residues)
  - Individual DSSP files: `results/structures/dssp/` (75 MB)
  - Script: `workflow/scripts/run_dssp.py`
- [x] D4: Validate structure quality (flag low-confidence models)
  - 1,382 structures validated; 77.4% high/very-high quality (pLDDT ≥ 80)
  - 63 flagged (4.6%): 4 very-low (<50), 21 few-usable (<30% residues >70)
  - Impact: 0/252 GroEL flagged, 3/266 HSP60 flagged — excellent core quality
  - Flagged proteins mostly from broader mito dataset; enriched for coil (57%)
  - Output: `results/structures/structure_quality_validation.tsv`, `flagged_low_quality_structures.tsv`
  - Report: `results/structures/quality_validation_report.txt`
  - Script: `workflow/scripts/validate_structure_quality.py`

### Module E: Structural Domain Architecture
- [x] E1: Obtain CATH domain assignments for pilot proteins
  - Used InterPro/Gene3D API (CATH mapped via InterPro)
  - **1,151/1,390 proteins (82.8%) have CATH assignments**
  - 2,141 total domains; mean 1.86 domains/protein
  - 50.7% single-domain, 29.4% two-domain, 10.5% three-domain
  - Top superfamilies: P-loop NTPases (105), Mito carrier (61), Rossmann (61)
  - 239 proteins lack CATH → need Chainsaw/Merizo (Step E2)
  - Coverage by dataset: GroEL 98%, HSP60 90.6%, matrix 89.9%, mito 79.3%
  - Output: `results/domains/cath_domain_assignments.tsv`, `cath_protein_summary.tsv`
  - Script: `workflow/scripts/get_cath_domains.py`
- [x] E2: For proteins without CATH assignments, run Chainsaw domain segmentation
  - 236/239 processed (3 had no AlphaFold structure)
  - Chainsaw v3 (Wells et al., 2024): deep learning domain boundary prediction
  - 18 proteins: 0 domains (small/disordered), 143: 1 domain, 54: 2 domains, 17: 3 domains, 4: 4 domains
  - Mean confidence: 0.821. Runtime: 8.5 min (CPU)
  - **Unified coverage: 1,387/1,390 (99.8%)** — CATH 1,151 + Chainsaw 236
  - Only 3 proteins unassigned (no AlphaFold model)
  - Output: `results/domains/chainsaw_domain_predictions.tsv`, `ml_domain_assignments.tsv`
  - Report: `results/domains/chainsaw_report.txt`
  - Script: `workflow/scripts/run_chainsaw_e2.py`
- [x] E3: Compute per-domain structural metrics
  - 2,141 domains processed; mean domain pLDDT: 92.1 (high confidence)
  - SS: 44.1% helix, 19.1% strand, 36.8% coil
  - 95.4% of domain residues have pLDDT > 70
  - Output: `results/domains/domain_structural_metrics.tsv`
  - Script: `workflow/scripts/compute_domain_structural_metrics.py`
- [x] E4: Run Foldseek clustering on pilot structures
  - 1,155 structural clusters from 1,382 proteins
  - 999 singletons (86.5%), 149 small (2-5), 7 medium (6-20)
  - **24 shared clusters** between GroEL + HSP60 substrates
  - Top shared: aldehyde dehydrogenases, thiolases, peroxiredoxins
  - Output: `results/domains/foldseek_clusters.tsv`, `foldseek_cluster_summary.txt`
  - Script: `workflow/scripts/analyze_foldseek.py`
- [x] E5: Compute domain architecture distribution per dataset
  - Alpha-beta dominates all datasets (60-71%)
  - GroEL top: TIM barrels (3.20.20.70) tied with P-loop NTPases
  - HSP60/mito top: P-loop NTPases (3.40.50.300)
  - Mean domains/protein: ~1.8-2.0 across datasets
  - Output: `results/domains/domain_distribution_summary.tsv`, `domain_distribution_report.txt`
  - Script: `workflow/scripts/domain_distribution_summary.py`

### Module F: N-domain vs C-region Stability Analysis
- [x] F1: Define three regions per protein
  - 1,151 proteins: 583 single-domain, 568 multi-domain
  - Output: `results/termini/region_boundaries.tsv`
- [x] F2: Compute sequence-derived metrics per region
  - 4,150 region rows (charge, hydrophobicity, aromaticity, polarity)
  - Output: `results/termini/sequence_metrics.tsv`
- [x] F3: Compute structure-derived metrics per region
  - 4,150 region rows (SS fractions, pLDDT confidence)
  - Output: `results/termini/structure_metrics.tsv`
- [x] F4: Compute contact order (Plaxco-Simons relative CO)
  - 5,296 region records. N-domains have HIGHER CO than C-regions (p=4.5e-20)
  - Output: `results/termini/contact_order.tsv`
- [x] F5: Within-protein paired comparison (multi-domain, n=567)
  - **KEY FINDINGS (Wilcoxon signed-rank):**
  - N-domains have higher relative contact order (+0.050, p=4.5e-20) — strongest signal
  - N-domains have higher pLDDT confidence (+2.1, p=3.6e-9)
  - N-domains have more beta-strand (+1.9%, p=6e-4)
  - No significant difference in hydrophobicity or charge
  - C-regions are ~70 residues longer on average
  - Output: `results/termini/n_vs_c_paired.tsv`
  - Script: `workflow/scripts/module_f_n_vs_c_analysis.py`

### Module G: Mitochondrial Targeting Analysis
- [x] G1: Transit peptide analysis (UniProt annotations)
  - 494/1,139 proteins (43.4%) have annotated transit peptides
  - Used UniProt REST API for transit peptide + signal peptide features
  - Output: `results/mts/uniprot_transit_signal_cache.tsv`
- [x] G2: Signal peptide exclusion
  - Only 15/1,139 (1.3%) have signal peptides — confirms mito targeting
- [x] G3: MitoCarta sub-mitochondrial localization (from existing data)
- [x] G4: Integrated targeting classification
  - HSP60 Tier 1: 124 high-confidence matrix (46.6%), 71 IM (26.7%), 56 non-canonical matrix (21.1%)
  - 168/266 HSP60 substrates (63.2%) have transit peptides
  - Output: `results/mts/combined_targeting.tsv`
- [x] G5: MTS vs domain boundary relationship
  - **84.4% (368/436): MTS is separate pre-domain extension** (not overlapping first domain)
  - 15.6% (68/436): MTS partially overlaps first domain (mean overlap 10.3 residues)
  - Median gap MTS-to-first-domain: 18 residues
  - Output: `results/mts/mts_domain_relationship.tsv`
  - Report: `results/mts/targeting_summary_report.txt`
  - Script: `workflow/scripts/module_g_mts_analysis.py`

### Module H: Comparative Statistics
- [x] H1: Pre-register primary hypotheses
  - 9 hypotheses across 3 goals documented
  - Output: `docs/PRIMARY_HYPOTHESES.md`
- [x] H2: Domain architecture enrichment tests
  - 242 Fisher's exact tests. 8 significant after hierarchical correction.
  - GroEL enriched: TIM barrel (OR=9.2, p=2.3e-6), Winged helix (OR=42.8, p=2.3e-6)
  - HSP60: mito carrier domain depleted (OR=0.16)
  - 79.7% of homolog pairs share the same top superfamily
  - Output: `results/stats/domain_enrichment.tsv`
- [x] H3: N-vs-C stability comparison statistics
  - 33 tests. 11 significant after correction.
  - N>C contact order confirmed (GroEL r=0.39, HSP60 r=0.52)
  - **Asymmetry is NOT substrate-specific** — background shows same pattern (p>0.14)
  - No difference between GroEL classes I/II/III (Kruskal-Wallis p>0.62)
  - Output: `results/stats/stability_comparisons.tsv`
- [x] H4: MTS and matrix targeting statistics
  - HSP60 substrates enriched for matrix (OR=3.29, p=1.6e-16)
  - MTS more prevalent in HSP60 substrates (OR=1.54, p=0.029)
  - MTS is pre-domain in 84.4% (binomial p=3.4e-51)
  - Output: `results/stats/targeting_stats.tsv`
- [x] H5: Multiple testing correction
  - 281 total tests. 22 significant at hierarchical level.
  - All 3 goal families passed family-level BH gate.
  - Output: `results/stats/corrected_pvalues.tsv`
  - Report: `results/stats/statistics_summary_report.txt`
  - Script: `workflow/scripts/module_h_comparative_stats.py`

### Module I: Visualization
- [x] I1: Domain distribution plots (Fig 1: CATH class, top superfamilies, domain counts)
- [x] I2: N-vs-C comparison plots (Fig 2: violin plots + heatmap; Fig 3: GroEL class comparison)
- [x] I3: MTS and targeting plots (Fig 4: targeting classification, gap histogram, scatter)
- [x] I4: Orthology visualization (Fig 5: Venn diagram + RCO correlation r=0.82)
- [x] I5: Summary figure (Fig 6: key findings from all 3 goals)
- All 6 figures in PDF+PNG at 300 DPI: `results/figures/`
- Script: `workflow/scripts/generate_figures.py`

### Module J: Documentation
- [x] J1: Master project documentation
  - 1,255 lines covering abstract, background, datasets, methods, results, decisions, limitations
  - Output: `docs/DOCUMENTATION.md`
- [x] J2: Methods and protocols (reproducibility guide)
  - ~700 lines with exact commands, data flow, QC checklist (70+ items)
  - Output: `docs/METHODS_AND_PROTOCOLS.md`
- [x] J3: Results narrative (manuscript-style)
  - All statistics, effect sizes, biological interpretations, figure references
  - Output: `docs/RESULTS_NARRATIVE.md`
- [x] J4: Pre-registered hypotheses document
  - 9 hypotheses across 3 goals
  - Output: `docs/PRIMARY_HYPOTHESES.md`

---

## Phase 2: Full-Scale (after pilot validation)

### Prerequisites
- [x] All Phase 1 modules complete and validated
- [x] 5 success criteria met — verified in `docs/PHASE1_VERIFICATION.md`
  - Domain coverage 99.8%, 3+ stability metrics, 82.5% orthology concordance
  - Cohen's d 0.26–0.61 for significant hypotheses, 100% MTS coverage

### Phase 1 Extended Analysis
- [x] F-ext: Re-run Module F with Chainsaw domains (642 multi-domain proteins, was 567)
  - New finding: N-domains significantly more hydrophobic (p=0.001, was borderline p=0.09)
  - All previous signals strengthened (contact order p=1.05e-20, pLDDT p=6.33e-10)
  - Output: `results/termini/region_boundaries_extended.tsv`, `n_vs_c_paired_extended.tsv`

### HPC Pipeline (prepared, ready to deploy)
- [x] Write Phase 2 HPC scripts and Snakemake workflow
  - `workflow/phase2/config.yaml` — Central configuration
  - `workflow/phase2/download_alphafold_full.py` — AlphaFold bulk download (~22 GB)
  - `workflow/phase2/run_foldseek_full.py` — Full-scale clustering (64 GB RAM)
  - `workflow/phase2/run_foldx.py` — FoldX stability with SLURM array jobs
  - `workflow/phase2/Snakefile` — 10-rule dependency graph
  - `workflow/phase2/README.md` — HPC deployment instructions

### Scale-Up Tasks (run on HPC)
- [ ] Deploy pipeline to HPC, configure paths in config.yaml
- [ ] Download full AlphaFold bulk structures (Human + E. coli, ~22 GB)
- [ ] Run Foldseek clustering at full scale (64 GB RAM, 16 CPUs)
- [ ] Run Chainsaw domain assignment on full proteomes
- [ ] Run FoldX stability calculations at scale (SLURM array, ~500 jobs)
- [ ] Repeat all Module E-I analyses on full proteomes
- [ ] Generate publication-ready figures

---

## Notes
- 4 GroEL accessions (P69408, Q99390, P62593, P29368) are plasmid/strain-specific, not in K-12 reference proteome
- Tools in conda env `proteomics`: `conda activate proteomics` before using mmseqs/foldseek
- MitoFates web server may be offline; TargetP 2.0 substitutes
- Disk constraint: 18 GB free. Full AlphaFold needs external storage or HPC.
