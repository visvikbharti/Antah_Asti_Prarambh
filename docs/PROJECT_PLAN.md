# Antah Asti Prarambh: Comparative Structural Proteomics Project Plan

**Project:** Comparative structural proteomics of chaperonin substrates
**Last updated:** 2026-03-11
**Status:** Phase 1 (Pilot) — Data acquisition partly complete

---

## Table of Contents

1. [Scientific Objectives](#1-scientific-objectives)
2. [Datasets](#2-datasets)
3. [Critical Scientific Decisions](#3-critical-scientific-decisions)
4. [Compute Environment](#4-compute-environment)
5. [Directory Structure](#5-directory-structure)
6. [Module A: Data Acquisition & Cleaning](#module-a-data-acquisition--cleaning)
7. [Module B: Dataset Construction](#module-b-dataset-construction)
8. [Module C: Orthology / Homology Layer](#module-c-orthology--homology-layer)
9. [Module D: Structure Acquisition & Indexing](#module-d-structure-acquisition--indexing)
10. [Module E: Structural Domain Architecture (CATH-based)](#module-e-structural-domain-architecture-cath-based)
11. [Module F: N-domain vs C-region Stability Analysis](#module-f-n-domain-vs-c-region-stability-analysis)
12. [Module G: Mitochondrial Targeting Analysis](#module-g-mitochondrial-targeting-analysis)
13. [Module H: Comparative Statistics](#module-h-comparative-statistics)
14. [Module I: Visualization & Figures](#module-i-visualization--figures)
15. [Dependency Graph & Implementation Order](#dependency-graph--implementation-order)
16. [Pilot vs Full-Scale Comparison](#pilot-vs-full-scale-comparison)
17. [Risk Register](#risk-register)

---

## 1. Scientific Objectives

This project performs a comparative structural proteomics analysis of chaperonin substrates across two organisms:

- **GroEL** (E. coli chaperonin, Kerner et al. 2005)
- **HSP60/HSPD1** (human mitochondrial chaperonin, Bruderer et al. 2020)

### Three primary aims:

| # | Aim | Question |
|---|-----|----------|
| 1 | **Structural domain distribution** | Do chaperonin substrates have distinct domain architectures (fold complexity, multi-domain arrangements) compared to non-substrate backgrounds? Are these patterns conserved between GroEL and HSP60? |
| 2 | **N-terminus vs C-terminus stability** | Are N-terminal structural domains of chaperonin substrates less stable or more disordered than C-terminal regions? Does co-translational interaction with chaperonins correlate with N-terminal instability? |
| 3 | **Mitochondrial matrix targeting** | Do HSP60 substrates have distinctive mitochondrial targeting sequence (MTS) properties? Do orthologous GroEL substrates that became mitochondrial acquire new targeting features? |

---

## 2. Datasets

Seven datasets are required. Current status:

| # | Dataset | Source | N | Status | File |
|---|---------|--------|---|--------|------|
| 1 | E. coli K-12 MG1655 proteome | UniProt UP000000625 | ~4,400 | **TO DOWNLOAD** | `data/raw/uniprot/ecoli_k12_proteome.fasta` |
| 2 | Human reference proteome | UniProt UP000005640 | ~20,600 | **TO DOWNLOAD** | `data/raw/uniprot/human_proteome.fasta` |
| 3 | Human mito proteome | MitoCarta 3.0 | 1,136 | **DONE** | `data/processed/human_mito_proteome.tsv` |
| 4 | GroEL substrates | Kerner 2005 Table S3 | 252 | **DONE** | `data/processed/groel_substrates_standardized.tsv` |
| 5 | HSP60 interactome (Tier 1) | Bruderer 2020 + filtering | 266 | **DONE** | `data/processed/hsp60_tier1_substrates.tsv` |
| 6 | 2-way orthology map (4 vs 5) | OrthoFinder / MMseqs2 | TBD | **TO COMPUTE** | `results/homology/ortholog_pairs.tsv` |
| 7 | Human mito matrix-only | MitoCarta 3.0 (is_matrix=1) | 525 | **DONE** | `data/processed/human_matrix_proteome.tsv` |

### Completed dataset details

**Dataset 3 — Human mito proteome** (1,136 proteins):
- Parsed from `data/raw/mitocarta/Human.MitoCarta3.0.xls` by `workflow/scripts/parse_mitocarta.py`
- Columns: `uniprot_id`, `gene_symbol`, `protein_name`, `mitocarta_score`, `sub_mito_localization`, `pathways`, `is_matrix`, `is_im`, `is_ims`, `is_om`
- Sub-compartment breakdown: Matrix 523, MIM 359, MOM 110, unknown 56, IMS 51, Membrane 34

**Dataset 4 — GroEL substrates** (252 proteins):
- Standardized from Kerner 2005 by `scripts/validate_uniprot_accessions.py`
- All accessions resolved against current UniProt via REST API; demerged accessions mapped to K-12 entries
- Columns: `original_accession`, `current_accession`, `accession_status`, `entry_name`, `gene_symbol`, `protein_name`, `organism`, `length`, `reviewed`, `groel_class` (I/II/III), `mass_kDa`, `subcellular_location`, `location_category`, `scop_folds`, `description_clean`

**Dataset 5 — HSP60 Tier 1** (266 proteins):
- Filtered from 325 total interactors by `scripts/filter_hsp60_interactome.py`
- Tier 1 criteria: MitoCarta+ AND median imputed SILAC ratio > 5
- Excluded (10 total): bait (HSPD1, HSPE1), co-chaperones (TRAP1, HSPA9, GRPEL1, DNAJA3), contaminants (IGHG2, TUBA1B, TUBB4B, TUBB)
- NDIC imputation: 2x 95th percentile of observed values per SILAC replicate
- Full filtering report: `data/processed/hsp60_filtering_report.txt`

**Dataset 7 — Matrix proteome** (525 proteins):
- Subset of Dataset 3 where `is_matrix == 1`
- 70 localization changes identified between MitoCarta 2.0 and 3.0 (52 Matrix->MIM reclassifications)
- Discrepancy report: `data/processed/mitocarta_summary_report.txt`

---

## 3. Critical Scientific Decisions

These decisions were made during the analysis phase and are binding for all downstream modules.

### Decision 1: Domain boundaries — CATH over InterPro

| Aspect | Choice | Rationale |
|--------|--------|-----------|
| Primary domain definition | **CATH proteome-level assignments** | CATH defines structural domains from 3D fold topology. Consistent, non-overlapping boundaries. Directly maps to AlphaFold structures. |
| Fallback for proteins without CATH | **Chainsaw / Merizo** | ML-based domain parsers that operate on AlphaFold PDB files. Use when CATH assignment unavailable. |
| InterPro role | **Complementary layer only** | InterPro domains mix sequence/structural/functional definitions. Overlapping boundaries. Used for functional annotation enrichment, NOT for domain counting or boundary analysis. |

### Decision 2: Stability measurement — NOT pLDDT alone

| Metric | Role | Notes |
|--------|------|-------|
| pLDDT | **Confidence / disorder proxy only** | AlphaFold's per-residue confidence. Low pLDDT indicates disorder or flexibility, NOT thermodynamic instability. |
| FoldX ddG | **Thermodynamic stability estimate** | Computed on AlphaFold structures. Caveat: FoldX was parameterized on experimental structures. |
| Contact order (absolute & relative) | **Folding kinetics proxy** | Higher contact order = slower folding. Relevant to chaperonin dependency. |
| Hydrophobic packing density | **Core quality metric** | Residue-level solvation energy or buried hydrophobic surface area. |
| B-factor equivalent (pLDDT-derived) | **Local flexibility** | Converted from pLDDT for compatibility with legacy structural tools. |

### Decision 3: HSP60 substrate filtering

- **Tier 1 = 266 proteins**: Must be in MitoCarta 3.0 AND have median imputed SILAC > 5
- **Excluded proteins**:
  - Bait: HSPD1 (HSP60 itself), HSPE1 (HSP10 co-chaperonin)
  - Co-chaperones: TRAP1, HSPA9 (mortalin/mtHSP70), GRPEL1, DNAJA3
  - Contaminants: IGHG2, tubulins (TUBA1B, TUBB4B, TUBB)
- **Tier 2** (49 proteins): SILAC > 2 but not MitoCarta+ or SILAC <= 5. Not used in primary analysis.
- **NDIC handling**: Imputed at 2x 95th percentile of observed values (high confidence enrichment)

### Decision 4: Orthology method

| Method | Role |
|--------|------|
| **OrthoFinder** | Primary. Run on full E. coli + human proteomes. Extract orthogroups, then intersect with substrate lists (Datasets 4 & 5). |
| **MMseqs2 RBH** | Supplementary validation. Reciprocal best hits between Dataset 4 and Dataset 5 directly. |
| **Intersection logic** | An ortholog pair is "conserved substrate" only if BOTH the E. coli gene is in Dataset 4 AND the human gene is in Dataset 5. |

### Decision 5: Controls — compartment-matched AND size-matched

| Organism | Control Tier | Population | Purpose |
|----------|-------------|------------|---------|
| E. coli | Primary | Cytoplasmic-only (excl. periplasm, membrane, secreted) | GroEL operates in cytoplasm |
| Human | Tier 1 (stringent) | Matrix-only (Dataset 7, 525 proteins) | HSP60 operates in matrix |
| Human | Tier 2 (intermediate) | All mito (Dataset 3, 1,136 proteins) | Broader mitochondrial context |
| Human | Tier 3 (permissive) | Full human proteome (Dataset 2) | Genome-wide baseline |

Size-matched controls: For each substrate set, generate 1,000 random draws from the background, matching the protein length distribution (within +/- 20% bins). This controls for the known correlation between protein size and domain count, fold complexity, and stability metrics.

### Decision 6: N-terminal analysis — three-region decomposition

For multi-domain proteins:
1. **Pre-domain disordered tail**: Residues 1 to (first domain start - 1). May include MTS.
2. **First structural domain**: As defined by CATH (or Chainsaw/Merizo fallback).
3. **Remainder**: All domains and linkers after the first domain.

For single-domain proteins: Handled separately. Compare N-terminal quarter vs C-terminal quarter of the single domain.

### Decision 7: MTS prediction strategy

| Source | Role |
|--------|------|
| **MitoCarta 3.0** | Ground truth for known mito proteins. No predictor needed. |
| **TargetP 2.0** | Primary predictor for novel/non-MitoCarta proteins. |
| **SignalP 6.0** | Distinguish MTS from signal peptides (similar N-terminal signals). |
| **DeepMito** | Sub-mitochondrial localization prediction. |
| **MitoFates** | Optional, server may be offline. |

### Decision 8: MitoCarta version

- Re-derived all annotations from MitoCarta **3.0** (Rath et al. 2021), not 2.0 (used in the original Bruderer 2020 HSP60 paper).
- 70 localization changes affect HSP60 interactors: 52 proteins reclassified from Matrix to MIM (primarily OXPHOS subunits and respiratory chain components).
- 6 proteins lost from MitoCarta entirely; 2 gained.
- All downstream analysis uses 3.0 annotations exclusively.

### Decision 9: Statistical framework

- **Hierarchical testing**: Pre-registered primary hypotheses tested first, exploratory tests second.
- **Primary analyses**: Multivariate (PCA on domain/stability feature matrix, regularized logistic regression for substrate vs background classification).
- **Secondary analyses**: Univariate comparisons (Mann-Whitney U, permutation tests) for individual features.
- **Multiple testing correction**: Permutation-based FDR (Westfall-Young) preferred over Bonferroni/BH for correlated features.
- **Effect sizes**: Report Cohen's d or rank-biserial correlation alongside p-values.

---

## 4. Compute Environment

### Local machine (Phase 1 — Pilot)

| Resource | Value |
|----------|-------|
| Hardware | Apple M1 (arm64), 8 GB RAM |
| Free disk | ~18 GB |
| OS | macOS (Darwin 25.2.0) |
| Python | 3.9 (Anaconda) |
| Shell | zsh |

### Installed software

| Package | Version | Purpose |
|---------|---------|---------|
| biopython | 1.78 | Sequence parsing, PDB I/O |
| pandas | 2.2.2 | Data manipulation |
| scipy | 1.9.2 | Statistical tests |
| matplotlib | (installed) | Plotting |
| seaborn | (installed) | Statistical visualization |
| snakemake | 7.32.4 | Workflow management |
| mkdssp | (installed) | Secondary structure assignment |
| openpyxl | (installed) | Excel file I/O |

### To install (conda bioconda)

```bash
conda install -c bioconda mmseqs2    # Sequence clustering, RBH
conda install -c bioconda foldseek   # Structural search/clustering
```

### HPC (Phase 2 — Full-scale only)

- Available at institute. Use only for:
  - Foldseek all-vs-all clustering of full proteomes
  - AlphaFold bulk structure retrieval/processing
  - InterProScan on full proteomes
  - OrthoFinder on full proteomes (if local RAM insufficient)

---

## 5. Directory Structure

```
Antah_Asti_Prarambh/
├── config/                          # Snakemake config, parameter files
│   ├── config.yaml                  # (TO CREATE) Master config
│   └── excluded_proteins.yaml       # (TO CREATE) Exclusion lists
├── data/
│   ├── raw/
│   │   ├── alphafold/               # AlphaFold PDB/mmCIF files
│   │   │   ├── ecoli/               # E. coli structures
│   │   │   └── human/               # Human structures
│   │   ├── custom/
│   │   │   ├── 12192_2020_1080_MOESM4_ESM.xlsx   # Bruderer 2020 supplement
│   │   │   ├── hsp60_interactome_clean.tsv        # Pre-cleaned HSP60 data
│   │   │   ├── kerner_2005_groel_interactors_clean.csv
│   │   │   └── kerner_2005_groel_interactors_table_s3.csv
│   │   ├── mitocarta/
│   │   │   └── Human.MitoCarta3.0.xls
│   │   ├── uniprot/                 # (TO POPULATE)
│   │   │   ├── ecoli_k12_proteome.fasta
│   │   │   ├── ecoli_k12_proteome.tsv
│   │   │   ├── human_proteome.fasta
│   │   │   └── human_proteome.tsv
│   │   └── cath/                    # (TO CREATE)
│   │       ├── cath-domain-list.txt
│   │       └── cath-b-newest-all.gz
│   ├── interim/                     # Intermediate computations
│   │   ├── orthofinder/             # OrthoFinder working directory
│   │   ├── mmseqs2/                 # MMseqs2 databases and results
│   │   ├── foldx/                   # FoldX energy calculations
│   │   └── dssp/                    # DSSP secondary structure files
│   └── processed/
│       ├── groel_substrates_standardized.tsv         # DONE
│       ├── hsp60_interactome_standardized.tsv        # DONE
│       ├── hsp60_tier1_substrates.tsv                # DONE
│       ├── hsp60_filtering_report.txt                # DONE
│       ├── human_mito_proteome.tsv                   # DONE
│       ├── human_matrix_proteome.tsv                 # DONE
│       ├── mitocarta_summary_report.txt              # DONE
│       ├── ecoli_cytoplasmic_background.tsv          # (TO CREATE)
│       ├── human_matrix_background.tsv               # (TO CREATE, size-matched)
│       ├── domain_architecture_groel.tsv             # (TO CREATE)
│       ├── domain_architecture_hsp60.tsv             # (TO CREATE)
│       ├── stability_features_groel.tsv              # (TO CREATE)
│       ├── stability_features_hsp60.tsv              # (TO CREATE)
│       └── mts_features_hsp60.tsv                    # (TO CREATE)
├── docs/
│   └── PROJECT_PLAN.md              # This file
├── logs/                            # Snakemake / script logs
├── results/
│   ├── dataset_build/               # Dataset assembly QC reports
│   ├── domains/                     # Domain architecture analysis outputs
│   ├── homology/                    # Orthology mapping outputs
│   ├── mts/                         # MTS prediction/analysis outputs
│   ├── stats/                       # Statistical test results
│   ├── termini/                     # N-term vs C-term stability outputs
│   └── figures/                     # Publication-quality figures
├── scripts/
│   ├── filter_hsp60_interactome.py                   # DONE
│   ├── validate_uniprot_accessions.py                # DONE
│   ├── download_uniprot_proteomes.py                 # (TO CREATE)
│   ├── download_alphafold_structures.py              # (TO CREATE)
│   ├── fetch_cath_assignments.py                     # (TO CREATE)
│   ├── run_chainsaw_merizo.py                        # (TO CREATE)
│   ├── compute_domain_architecture.py                # (TO CREATE)
│   ├── compute_stability_features.py                 # (TO CREATE)
│   ├── compute_contact_order.py                      # (TO CREATE)
│   ├── compute_mts_features.py                       # (TO CREATE)
│   ├── run_orthofinder.py                            # (TO CREATE)
│   ├── run_mmseqs2_rbh.py                            # (TO CREATE)
│   ├── build_controls.py                             # (TO CREATE)
│   ├── run_statistics.py                             # (TO CREATE)
│   └── generate_figures.py                           # (TO CREATE)
└── workflow/
    ├── envs/                        # Conda environment YAMLs
    ├── rules/                       # Snakemake rule files
    │   ├── data_acquisition.smk     # (TO CREATE)
    │   ├── structures.smk           # (TO CREATE)
    │   ├── domains.smk              # (TO CREATE)
    │   ├── stability.smk            # (TO CREATE)
    │   ├── homology.smk             # (TO CREATE)
    │   ├── mts.smk                  # (TO CREATE)
    │   └── statistics.smk           # (TO CREATE)
    └── scripts/
        └── parse_mitocarta.py       # DONE
```

---

## Module A: Data Acquisition & Cleaning

### Purpose
Obtain and standardize all raw input datasets. Validate identifiers, resolve obsolete accessions, and produce clean TSV files with consistent column naming.

### Status: MOSTLY DONE

| Task | Status | Script | Notes |
|------|--------|--------|-------|
| GroEL substrates (Kerner 2005) | DONE | `scripts/validate_uniprot_accessions.py` | 252 proteins. All accessions resolved via UniProt REST API. Demerged accessions mapped to K-12 entries. |
| HSP60 interactome (Bruderer 2020) | DONE | `scripts/filter_hsp60_interactome.py` | 325 total -> 266 Tier 1. NDIC imputation, tiered quality assignment. |
| MitoCarta 3.0 | DONE | `workflow/scripts/parse_mitocarta.py` | 1,136 mito proteins, 525 matrix. Cross-referenced with HSP60 file. |
| UniProt E. coli proteome | TO DO | `scripts/download_uniprot_proteomes.py` | See Module B |
| UniProt human proteome | TO DO | `scripts/download_uniprot_proteomes.py` | See Module B |

### Input files (raw)

| File | Path |
|------|------|
| Kerner 2005 clean CSV | `data/raw/custom/kerner_2005_groel_interactors_clean.csv` |
| Kerner 2005 Table S3 CSV | `data/raw/custom/kerner_2005_groel_interactors_table_s3.csv` |
| Bruderer 2020 supplement | `data/raw/custom/12192_2020_1080_MOESM4_ESM.xlsx` |
| HSP60 interactome TSV | `data/raw/custom/hsp60_interactome_clean.tsv` |
| MitoCarta 3.0 XLS | `data/raw/mitocarta/Human.MitoCarta3.0.xls` |

### Output files (processed)

| File | Path | Rows |
|------|------|------|
| GroEL substrates | `data/processed/groel_substrates_standardized.tsv` | 252 |
| HSP60 full standardized | `data/processed/hsp60_interactome_standardized.tsv` | 325 |
| HSP60 Tier 1 | `data/processed/hsp60_tier1_substrates.tsv` | 266 |
| HSP60 filtering report | `data/processed/hsp60_filtering_report.txt` | — |
| Human mito proteome | `data/processed/human_mito_proteome.tsv` | 1,136 |
| Human matrix proteome | `data/processed/human_matrix_proteome.tsv` | 525 |
| MitoCarta report | `data/processed/mitocarta_summary_report.txt` | — |

### Dependencies
None (first module).

### Compute
Local only.

---

## Module B: Dataset Construction

### Purpose
Download full reference proteomes from UniProt. Build compartment-matched and size-matched background sets for statistical controls.

### Status: PARTLY DONE (mito/matrix subsets exist; proteomes not yet downloaded)

### Task B1: Download E. coli K-12 MG1655 proteome

**Script:** `scripts/download_uniprot_proteomes.py`
**Method:** UniProt REST API or direct FTP download
**Proteome ID:** UP000000625
**Taxonomy:** 83333 (E. coli K-12)
**Format:** Download both FASTA (for OrthoFinder/MMseqs2) and TSV (with annotations)

```
# TSV query fields:
# accession, id, gene_names, protein_name, length, mass,
# cc_subcellular_location, go_c, go_f, reviewed
```

**Output files:**
- `data/raw/uniprot/ecoli_k12_proteome.fasta`
- `data/raw/uniprot/ecoli_k12_proteome.tsv`

**Expected count:** ~4,400 reviewed proteins (Swiss-Prot only; exclude TrEMBL)

### Task B2: Download human reference proteome

**Proteome ID:** UP000005640
**Taxonomy:** 9606
**Note:** Download canonical sequences only (no isoforms) for Phase 1. Isoforms in Phase 2 if needed.

**Output files:**
- `data/raw/uniprot/human_proteome.fasta`
- `data/raw/uniprot/human_proteome.tsv`

**Expected count:** ~20,600 reviewed proteins

### Task B3: Build E. coli cytoplasmic background

**Script:** `scripts/build_controls.py`
**Input:** `data/raw/uniprot/ecoli_k12_proteome.tsv` + `data/processed/groel_substrates_standardized.tsv`
**Logic:**
1. Parse `cc_subcellular_location` from UniProt TSV.
2. Keep proteins annotated as "Cytoplasm" or "Cytosol".
3. Exclude proteins annotated as "Periplasm", "Cell outer membrane", "Cell inner membrane", "Secreted".
4. Exclude proteins already in Dataset 4 (GroEL substrates).
5. For "unknown" location: include in a separate "all-cytoplasmic-plus-unknown" set for sensitivity analysis.

**Output:**
- `data/processed/ecoli_cytoplasmic_background.tsv`
- `data/processed/ecoli_cyto_plus_unknown_background.tsv`

### Task B4: Build size-matched control draws

**Script:** `scripts/build_controls.py` (same script, separate function)
**Method:**
1. For each substrate set (GroEL, HSP60 Tier 1), compute the protein length distribution in 50-residue bins.
2. From the corresponding background, draw 1,000 bootstrap samples matching this length distribution (within +/- 20% per bin).
3. Store draws as a compressed numpy array or HDF5 for efficient downstream permutation testing.

**Output:**
- `data/interim/size_matched_draws_groel.npz`
- `data/interim/size_matched_draws_hsp60.npz`

### Tools
- `requests` or `urllib` (UniProt API)
- `pandas`
- `numpy`

### Compute
Local.

### Dependencies
Module A must be complete.

---

## Module C: Orthology / Homology Layer

### Purpose
Establish ortholog relationships between E. coli and human proteins. Identify "conserved chaperonin substrates" — proteins that are GroEL substrates in E. coli AND HSP60 substrates in human.

### Task C1: OrthoFinder (primary orthology)

**Script:** `scripts/run_orthofinder.py`
**Input:**
- `data/raw/uniprot/ecoli_k12_proteome.fasta`
- `data/raw/uniprot/human_proteome.fasta`

**Method:**
1. Prepare input FASTA files in OrthoFinder format (one file per species).
2. Run OrthoFinder with DIAMOND for sequence search (default).
3. Extract orthogroups and 1:1 orthologs.
4. Intersect with substrate lists:
   - For each orthogroup containing at least one GroEL substrate AND one HSP60 Tier 1 substrate, flag as "conserved substrate orthogroup."
   - For 1:1 ortholog pairs: flag pairs where BOTH members are substrates.

**Parameters:**
- OrthoFinder default DIAMOND e-value (1e-3)
- MCL inflation parameter: default (1.5)

**Output:**
- `data/interim/orthofinder/` (OrthoFinder working directory)
- `results/homology/orthogroups_all.tsv`
- `results/homology/ortholog_pairs_1to1.tsv`
- `results/homology/conserved_substrate_orthogroups.tsv`
- `results/homology/conserved_substrate_pairs.tsv`

**Pilot vs Full-scale:**
- Pilot: Run OrthoFinder on full proteomes locally (should fit in 8 GB RAM with DIAMOND). If RAM is insufficient, use MMseqs2 RBH only.
- Full-scale: Same, but on HPC if needed.

### Task C2: MMseqs2 reciprocal best hits (supplementary)

**Script:** `scripts/run_mmseqs2_rbh.py`
**Input:**
- `data/processed/groel_substrates_standardized.tsv` (extract FASTA from UniProt IDs)
- `data/processed/hsp60_tier1_substrates.tsv` (extract FASTA from UniProt IDs)

**Method:**
1. Create MMseqs2 databases for each substrate list.
2. Run `mmseqs rbh` between the two databases.
3. Filter by e-value < 1e-5, coverage > 50%.
4. Compare RBH pairs against OrthoFinder 1:1 orthologs for concordance.

**Parameters:**
- `--min-seq-id 0.2` (permissive, since E. coli vs human is ancient divergence)
- `-e 1e-5`
- `-c 0.5 --cov-mode 0` (bidirectional coverage)

**Output:**
- `data/interim/mmseqs2/` (MMseqs2 databases)
- `results/homology/mmseqs2_rbh_pairs.tsv`
- `results/homology/orthology_concordance.tsv`

### Tools
- OrthoFinder (to install: `conda install -c bioconda orthofinder`)
- MMseqs2 (to install: `conda install -c bioconda mmseqs2`)
- DIAMOND (installed automatically with OrthoFinder)

### Compute
- Pilot: Local (OrthoFinder on 2 proteomes should take ~20 min on M1; MMseqs2 RBH on ~500 proteins is trivial)
- Full-scale: HPC if OrthoFinder exceeds local RAM

### Dependencies
Module B (requires full proteome FASTA files).

---

## Module D: Structure Acquisition & Indexing

### Purpose
Download AlphaFold predicted structures for all proteins in Datasets 4, 5, and control sets. Index them for downstream domain parsing, stability calculations, and visualization.

### Task D1: Download AlphaFold structures (Pilot)

**Script:** `scripts/download_alphafold_structures.py`
**Input:** UniProt IDs from Datasets 4, 5, 7, and E. coli cytoplasmic background
**Method:**
1. For each UniProt accession, download the AlphaFold predicted structure from `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb`.
2. Also download the per-residue pLDDT JSON (or extract from PDB B-factor column).
3. Verify file integrity (check PDB can be parsed by Biopython).
4. Log missing structures (no AlphaFold model available).

**Output:**
- `data/raw/alphafold/ecoli/AF-{ID}-F1-model_v4.pdb` (one per protein)
- `data/raw/alphafold/human/AF-{ID}-F1-model_v4.pdb` (one per protein)
- `data/interim/alphafold_download_log.tsv` (ID, status, file_size, pLDDT_mean)

**Pilot scope:**
- GroEL substrates: 252 structures
- HSP60 Tier 1: 266 structures
- Matrix background: 525 structures
- E. coli cytoplasmic background: ~3,500 structures (if available)
- **Total pilot: ~4,500 structures, estimated ~2-4 GB disk**

**Full-scale:** Full E. coli + full human proteomes (~25,000 structures, ~10-15 GB). Use AlphaFold bulk download from EBI FTP.

### Task D2: Run DSSP for secondary structure assignment

**Script:** Part of `scripts/compute_stability_features.py`
**Input:** AlphaFold PDB files from D1
**Method:** Run `mkdssp` on each PDB file.
**Output:** `data/interim/dssp/{UNIPROT_ID}.dssp`

### Task D3: Build structure index

**Script:** Part of `scripts/download_alphafold_structures.py`
**Output:** `data/processed/structure_index.tsv`

Columns: `uniprot_id`, `pdb_path`, `pdb_available`, `n_residues`, `mean_plddt`, `median_plddt`, `frac_disordered` (pLDDT < 50), `frac_confident` (pLDDT > 70)

### Tools
- `requests` (HTTP download)
- `biopython` (PDB parsing)
- `mkdssp` (secondary structure)

### Compute
- Pilot: Local. ~4,500 downloads at ~100 KB each = ~450 MB. Rate-limited to ~2 requests/sec.
- Full-scale: Use AlphaFold EBI FTP bulk download on HPC.

### Dependencies
Module A (needs UniProt IDs from processed datasets).
Module B (needs background protein lists for structure download).

---

## Module E: Structural Domain Architecture (CATH-based)

### Purpose
Assign structural domain boundaries to every protein. Compute domain architecture features (domain count, fold types, multi-domain topology). Compare substrates vs backgrounds.

### Task E1: Fetch CATH proteome-level assignments

**Script:** `scripts/fetch_cath_assignments.py`
**Input:** UniProt IDs from all datasets
**Method:**
1. Download CATH-Gene3D proteome-level domain assignments from the CATH FTP:
   - `ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest/sequence-data/`
   - Or use CATH API: `https://www.cathdb.info/version/latest/api/rest/uniprot_to_funfam/{UNIPROT_ID}`
2. For each protein, extract:
   - Number of CATH domains
   - CATH superfamily codes (e.g., 3.40.50.300)
   - Domain boundaries (residue ranges)
   - CATH class (mainly alpha, mainly beta, alpha-beta, few secondary structures)
3. Map SCOP folds from Kerner 2005 (already in `groel_substrates_standardized.tsv`) to CATH equivalents where possible.

**Output:**
- `data/raw/cath/cath_assignments_pilot.tsv` (raw CATH data for pilot proteins)
- `data/processed/cath_domain_boundaries.tsv`

Columns: `uniprot_id`, `n_domains`, `domain_id`, `cath_superfamily`, `cath_class`, `start_residue`, `end_residue`, `domain_length`

### Task E2: Chainsaw/Merizo fallback for proteins without CATH

**Script:** `scripts/run_chainsaw_merizo.py`
**Input:** AlphaFold PDB files for proteins not covered by CATH
**Method:**
1. Identify proteins with AlphaFold structures but no CATH assignment.
2. Run Chainsaw (https://github.com/JudeWells/chainsaw) on the PDB files.
3. As additional validation, run Merizo (https://github.com/psipred/merizo) on the same PDBs.
4. Accept domain boundaries where Chainsaw and Merizo agree (within 10 residues).
5. For disagreements, use Chainsaw as primary (better benchmarked on CATH).

**Output:**
- `data/interim/chainsaw_domains.tsv`
- `data/interim/merizo_domains.tsv`
- `data/processed/fallback_domain_boundaries.tsv`

### Task E3: Compute domain architecture features

**Script:** `scripts/compute_domain_architecture.py`
**Input:**
- `data/processed/cath_domain_boundaries.tsv`
- `data/processed/fallback_domain_boundaries.tsv`
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- Background protein lists

**Features computed per protein:**

| Feature | Description |
|---------|-------------|
| `n_domains` | Total number of structural domains |
| `domain_class_distribution` | Fraction of each CATH class (alpha, beta, alpha-beta) |
| `has_tandem_repeats` | Whether protein contains tandem repeated domains |
| `max_domain_length` | Length of largest domain |
| `min_domain_length` | Length of smallest domain |
| `interdomain_linker_total` | Total residues in linker regions |
| `frac_unassigned` | Fraction of residues not in any domain |
| `topology_string` | Ordered CATH superfamily string (e.g., "3.40-2.60-3.40") |
| `is_single_domain` | Boolean |
| `is_multi_domain` | Boolean |
| `fold_complexity` | CATH topology count (unique T-level codes) |

**Output:**
- `data/processed/domain_architecture_groel.tsv`
- `data/processed/domain_architecture_hsp60.tsv`
- `data/processed/domain_architecture_backgrounds.tsv`
- `results/domains/domain_feature_summary.tsv`

### Tools
- CATH API or FTP data files
- Chainsaw (Python package or clone from GitHub)
- Merizo (Python package or clone from GitHub)
- `biopython`, `pandas`, `numpy`

### Compute
- Pilot: Local. CATH API queries are lightweight. Chainsaw/Merizo are CPU-only and run in seconds per structure.
- Full-scale: Batch processing on HPC for full proteome Chainsaw/Merizo runs.

### Dependencies
Module D (requires AlphaFold structures for Chainsaw/Merizo).
Module B (requires background protein lists).

---

## Module F: N-domain vs C-region Stability Analysis

### Purpose
Test whether chaperonin substrates have less stable N-terminal regions compared to C-terminal regions. Measure stability using multiple orthogonal metrics (NOT pLDDT alone).

### Task F1: Compute per-residue and per-region stability features

**Script:** `scripts/compute_stability_features.py`
**Input:**
- AlphaFold PDB files (`data/raw/alphafold/`)
- Domain boundaries (`data/processed/cath_domain_boundaries.tsv`, `data/processed/fallback_domain_boundaries.tsv`)
- DSSP files (`data/interim/dssp/`)

**Metrics computed per residue:**

| Metric | Tool/Method | Notes |
|--------|-------------|-------|
| pLDDT | Extract from AlphaFold PDB B-factor | Confidence/disorder proxy only. NOT stability. |
| Contact order (absolute) | Custom: mean sequence separation of contacting residues (C-beta < 8 A) | Higher = slower folding |
| Contact order (relative) | Absolute contact order / chain length | Normalized for size |
| Local packing density | Count C-beta atoms within 10 A sphere per residue | Core quality |
| Hydrophobic burial | Fraction of hydrophobic residues with <20% solvent accessibility (from DSSP) | Core quality |
| Secondary structure content | From DSSP: fraction helix, sheet, coil per region | Fold character |

**Aggregation into three regions (for multi-domain proteins):**

| Region | Definition | Metrics aggregated |
|--------|------------|--------------------|
| Pre-domain tail | Residues 1 to (first_domain_start - 1) | Mean pLDDT, fraction disordered, length |
| First domain (N-domain) | First CATH domain boundaries | All 6 metrics above |
| Remainder (C-region) | Everything after first domain end | All 6 metrics above |

**For single-domain proteins:**
- Split at midpoint (or quartiles): compare N-terminal quarter vs C-terminal quarter.
- Also compute whole-domain metrics for comparison.

**Output:**
- `data/processed/stability_features_groel.tsv`
- `data/processed/stability_features_hsp60.tsv`
- `data/processed/stability_features_backgrounds.tsv`

Columns: `uniprot_id`, `n_domains`, `pre_domain_length`, `pre_domain_mean_plddt`, `pre_domain_frac_disordered`, `n_domain_mean_plddt`, `n_domain_contact_order_abs`, `n_domain_contact_order_rel`, `n_domain_packing_density`, `n_domain_hydrophobic_burial`, `n_domain_helix_frac`, `n_domain_sheet_frac`, `c_region_mean_plddt`, `c_region_contact_order_abs`, `c_region_contact_order_rel`, `c_region_packing_density`, `c_region_hydrophobic_burial`, `c_region_helix_frac`, `c_region_sheet_frac`, `n_minus_c_plddt`, `n_minus_c_contact_order`, `n_minus_c_packing`

### Task F2: FoldX stability estimation (Phase 2)

**Deferred to Phase 2.** FoldX requires license and is computationally intensive.

**Method (when executed):**
1. Repair AlphaFold structures with `FoldX --command=RepairPDB`.
2. Compute `FoldX --command=Stability` for per-domain free energy.
3. Compute `FoldX --command=AnalyseComplex` for domain interface energy in multi-domain proteins.

**Output:** `data/processed/foldx_stability_features.tsv`

### Tools
- `biopython` (PDB parsing, distance calculations)
- `mkdssp` (solvent accessibility, secondary structure)
- `numpy`, `scipy` (contact order calculation)
- FoldX (Phase 2 only, requires license)

### Compute
- Pilot: Local. Contact order and packing density calculations are O(N^2) per protein but N < 2,000 residues typically.
- Full-scale: HPC for FoldX runs on full proteome.

### Dependencies
Module D (AlphaFold structures, DSSP files).
Module E (domain boundaries for region decomposition).

---

## Module G: Mitochondrial Targeting Analysis

### Purpose
Characterize MTS properties of HSP60 substrates. Test whether conserved chaperonin substrates (GroEL in E. coli -> HSP60 in human) acquired mitochondrial targeting. Analyze MTS relationship to structural domains.

### Task G1: Extract MTS annotations from MitoCarta 3.0

**Script:** `scripts/compute_mts_features.py`
**Input:**
- `data/processed/human_mito_proteome.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- AlphaFold structures (for MTS structural context)
- Domain boundaries (to check MTS/domain overlap)

**Method:**
1. For all MitoCarta proteins: MTS is confirmed. Use UniProt "TRANSIT" annotation to get cleavage site position.
2. Fetch transit peptide annotations from UniProt for all HSP60 Tier 1 proteins.
3. If UniProt transit peptide annotation is missing, use TargetP 2.0 prediction.

**Features computed per protein:**

| Feature | Description |
|---------|-------------|
| `has_mts` | Boolean: transit peptide annotated or predicted |
| `mts_length` | Length of MTS in residues |
| `mts_charge` | Net charge of MTS sequence |
| `mts_hydrophobicity` | Mean Kyte-Doolittle hydrophobicity |
| `mts_amphipathic_helix` | Predicted amphipathic helix content (from sequence) |
| `mts_cleavage_site` | Position of R-2 motif or predicted cleavage |
| `mts_overlaps_domain` | Whether MTS extends into first structural domain |
| `first_domain_start` | Start position of first CATH domain |
| `mts_to_domain_gap` | Distance between MTS cleavage and first domain start |

### Task G2: MTS prediction for non-MitoCarta proteins

**Script:** Part of `scripts/compute_mts_features.py`
**Input:** FASTA sequences for proteins not in MitoCarta
**Tools:**
- **TargetP 2.0**: Primary MTS predictor. Use command-line version if available, or web server API.
- **SignalP 6.0**: Distinguish MTS from signal peptides.
- **DeepMito**: Sub-mitochondrial localization prediction.
- **MitoFates**: Optional (server may be offline). If available, use for MPP cleavage site prediction.

**Output:**
- `results/mts/targetp_predictions.tsv`
- `results/mts/signalp_predictions.tsv`
- `data/processed/mts_features_hsp60.tsv`

### Task G3: Cross-organism MTS analysis

**Script:** Part of `scripts/compute_mts_features.py`
**Input:**
- `results/homology/conserved_substrate_pairs.tsv` (from Module C)
- `data/processed/mts_features_hsp60.tsv`
- `data/processed/groel_substrates_standardized.tsv`

**Analysis:**
1. For each conserved substrate pair (GroEL substrate in E. coli, HSP60 substrate in human):
   - Does the human ortholog have an MTS?
   - How long is the MTS relative to the mature protein?
   - Does the MTS overlap with any structural domain?
2. For GroEL substrates without HSP60 orthologs:
   - Are the human orthologs (if any) not mitochondrial? Or are they mitochondrial but not HSP60 substrates?
3. Summary statistics: What fraction of conserved substrate pairs have MTS?

**Output:**
- `results/mts/conserved_substrate_mts_analysis.tsv`
- `results/mts/mts_analysis_summary.txt`

### Tools
- `biopython` (sequence feature extraction)
- TargetP 2.0 (command-line or web)
- SignalP 6.0 (command-line or web)
- DeepMito (web or local)
- `pandas`, `numpy`

### Compute
- Pilot: Local. Web-based predictions for ~266 proteins are feasible.
- Full-scale: Install command-line TargetP/SignalP on HPC for full proteome.

### Dependencies
Module C (ortholog pairs for cross-organism analysis).
Module E (domain boundaries for MTS/domain overlap).
Module D (AlphaFold structures for MTS structural context).

---

## Module H: Comparative Statistics

### Purpose
Perform all statistical tests. Compare substrate distributions against backgrounds. Test primary hypotheses with pre-registered framework.

### Task H1: Pre-registered primary hypotheses

| # | Hypothesis | Test type | Substrate set | Background | Features |
|---|-----------|-----------|---------------|------------|----------|
| H1 | Chaperonin substrates are enriched for multi-domain proteins | Logistic regression | GroEL (252) | E. coli cytoplasmic | `n_domains`, `is_multi_domain` |
| H2 | Same as H1 in human | Logistic regression | HSP60 Tier 1 (266) | Matrix (525) | `n_domains`, `is_multi_domain` |
| H3 | N-terminal domains of substrates are less stable than C-terminal regions | Paired test within substrates | GroEL + HSP60 | Self (paired N vs C) | `n_minus_c_contact_order`, `n_minus_c_packing` |
| H4 | N-terminal instability is greater in substrates than backgrounds | Interaction test | GroEL vs cytoplasmic | Background | `n_minus_c` delta features |
| H5 | Conserved substrates (orthologs in both sets) have extreme domain/stability features | Comparison within substrates | Conserved vs non-conserved substrates | — | All domain + stability features |
| H6 | HSP60 substrates have distinct MTS properties vs matrix background | Logistic regression | HSP60 Tier 1 | Matrix non-substrates | MTS features |

### Task H2: Multivariate analyses (primary)

**Script:** `scripts/run_statistics.py`
**Input:**
- `data/processed/domain_architecture_*.tsv`
- `data/processed/stability_features_*.tsv`
- `data/processed/mts_features_hsp60.tsv`
- Size-matched control draws

**Methods:**
1. **PCA** on the combined feature matrix (domain + stability features) for each organism. Visualize substrate vs background separation.
2. **Regularized logistic regression** (L1/elastic net) to classify substrate vs background. Use cross-validation (5-fold) to estimate accuracy and identify most discriminative features.
3. **Permutation importance** to rank features.

**Output:**
- `results/stats/pca_groel.tsv`, `results/stats/pca_hsp60.tsv`
- `results/stats/logistic_regression_groel.tsv`, `results/stats/logistic_regression_hsp60.tsv`
- `results/stats/feature_importance.tsv`

### Task H3: Univariate analyses (secondary)

**Methods:**
1. **Mann-Whitney U** for each feature, substrate vs background.
2. **Permutation tests** (10,000 permutations) for features where distributional assumptions are unclear.
3. **Westfall-Young step-down** for multiple testing correction (accounts for feature correlation).
4. **Effect sizes**: Cohen's d (for approximately normal features) or rank-biserial correlation (for non-parametric).

**Output:**
- `results/stats/univariate_tests_groel.tsv`
- `results/stats/univariate_tests_hsp60.tsv`
- `results/stats/multiple_testing_correction.tsv`

### Task H4: Cross-organism conservation analysis

Compare feature distributions between:
- Conserved substrates (in both GroEL and HSP60 sets)
- GroEL-only substrates
- HSP60-only substrates
- Non-substrate orthologs

**Output:**
- `results/stats/cross_organism_comparison.tsv`

### Tools
- `scipy.stats` (Mann-Whitney, permutation tests)
- `scikit-learn` (PCA, logistic regression, cross-validation)
- `statsmodels` (for detailed regression output, if needed)
- `numpy` (permutation-based correction)

### Compute
- Pilot: Local. All statistical tests are lightweight.
- Full-scale: Same. Statistics do not require HPC.

### Dependencies
Module E (domain architecture features).
Module F (stability features).
Module G (MTS features).
Module C (ortholog pairs for cross-organism analysis).

---

## Module I: Visualization & Figures

### Purpose
Generate publication-quality figures for all analyses.

### Planned figures

| # | Figure | Module source | Type |
|---|--------|---------------|------|
| 1 | Dataset overview: protein counts, overlaps, Venn diagrams | B | Summary |
| 2 | Domain count distributions: substrates vs backgrounds (E. coli and human, side by side) | E | Violin/box |
| 3 | CATH class distribution: pie/bar charts for substrates vs backgrounds | E | Bar chart |
| 4 | Domain architecture topology strings: heatmap of most common topologies | E | Heatmap |
| 5 | N-terminal vs C-terminal stability: paired scatter/delta histograms | F | Scatter + histogram |
| 6 | pLDDT profiles along protein length: average profile for substrates vs backgrounds | F | Line plot |
| 7 | Contact order by region: box plots for pre-domain, N-domain, C-region | F | Box plot |
| 8 | PCA biplots: substrates vs backgrounds in PC1-PC2 space | H | Scatter |
| 9 | Feature importance: bar chart from regularized logistic regression | H | Bar chart |
| 10 | MTS length distribution: HSP60 substrates vs matrix background | G | Histogram |
| 11 | MTS-to-domain gap: distribution and relationship to first domain stability | G | Scatter + histogram |
| 12 | Cross-organism comparison: conserved vs species-specific substrate features | H | Grouped violin |
| 13 | Orthology Sankey diagram: flow from GroEL substrates to human orthologs | C | Sankey |

### Style specifications
- Backend: `matplotlib` + `seaborn`
- Figure size: single column (3.5 inch) or double column (7.0 inch)
- DPI: 300 for raster, vector (PDF/SVG) preferred
- Color palette: colorblind-safe (use `seaborn` "colorblind" or custom palette)
- Font: Arial or Helvetica, 8-10 pt for axis labels

### Output
- `results/figures/fig{N}_{name}.pdf` (vector)
- `results/figures/fig{N}_{name}.png` (raster, 300 DPI)

### Tools
- `matplotlib`, `seaborn`
- `adjustText` (for label placement, if needed)

### Compute
Local only.

### Dependencies
All analysis modules (E, F, G, H) must be complete.

---

## Dependency Graph & Implementation Order

### ASCII dependency graph

```
Module A (Data Acquisition)         STATUS: MOSTLY DONE
  |
  v
Module B (Dataset Construction)     STATUS: TO DO
  |
  +---> Module C (Orthology)        STATUS: TO DO
  |       |
  +---> Module D (Structures)       STATUS: TO DO
  |       |
  |       v
  |     Module E (Domains)          STATUS: TO DO
  |       |
  |       v
  |     Module F (Stability)        STATUS: TO DO
  |       |
  |       +-----+
  |              |
  v              v
Module G (MTS Analysis)             STATUS: TO DO
  |              |
  v              v
Module H (Statistics)               STATUS: TO DO
  |
  v
Module I (Figures)                  STATUS: TO DO
```

### Detailed dependencies

| Module | Depends on | Can start when |
|--------|-----------|----------------|
| A | None | Now |
| B | A | A is complete (READY — A is mostly done) |
| C | B | B1, B2 complete (proteome FASTAs downloaded) |
| D | A, B | A complete (for UniProt IDs), B3 complete (for background lists) |
| E | D | D1 complete (AlphaFold structures), plus CATH data fetched |
| F | D, E | D1 + D2 (structures + DSSP) and E1-E2 (domain boundaries) |
| G | C, D, E | C (ortholog pairs), D (structures), E (domain boundaries) |
| H | E, F, G, C | All feature computation modules complete |
| I | H | All statistics complete |

### Parallelization opportunities

The following can run **in parallel** once Module B is complete:

```
B complete
  |
  +---> C (Orthology)        [independent track]
  |
  +---> D (Structures)       [independent track]
         |
         +---> E (Domains)   [depends on D]
         |
         +---> D2 (DSSP)     [depends on D1]
                |
                +---> F (Stability) [depends on D + E]
```

Module C and Modules D/E/F are fully independent and can execute simultaneously.

Module G requires outputs from C, D, and E, so it is the convergence point.

### Recommended execution order

```
WEEK 1:
  B1: Download E. coli proteome         (30 min)
  B2: Download human proteome           (30 min)
  B3: Build E. coli cytoplasmic bg      (1 hr)
  D1: Download AlphaFold pilot structs  (2-4 hrs, rate limited)

WEEK 2:
  [parallel track 1]
  C1: OrthoFinder                       (20-60 min)
  C2: MMseqs2 RBH                       (5 min)

  [parallel track 2]
  D2: DSSP on all structures            (1-2 hrs)
  E1: Fetch CATH assignments            (1-2 hrs)
  E2: Chainsaw/Merizo fallback          (1-2 hrs)
  E3: Compute domain architecture       (30 min)

WEEK 3:
  F1: Compute stability features        (2-4 hrs)
  B4: Build size-matched controls       (30 min)

WEEK 4:
  G1: Extract MTS annotations           (1 hr)
  G2: MTS prediction (TargetP, etc.)    (2-4 hrs, web submissions)
  G3: Cross-organism MTS analysis       (30 min)

WEEK 5:
  H1-H4: All statistical analyses      (2-4 hrs)

WEEK 6:
  I: Figure generation                  (1-2 days)
  Documentation, review, revision
```

---

## Pilot vs Full-Scale Comparison

| Aspect | Phase 1 (Pilot) | Phase 2 (Full-scale) |
|--------|-----------------|----------------------|
| **Proteins analyzed** | ~900 (252 GroEL + 266 HSP60 + ~400 background) | ~25,000 (full E. coli + human proteomes) |
| **AlphaFold structures** | ~4,500 (substrates + backgrounds) | ~25,000 |
| **Compute location** | All local (M1, 8 GB) | HPC for heavy tasks |
| **OrthoFinder** | Full proteomes locally (DIAMOND) | Same, HPC if needed |
| **Foldseek** | Not used | All-vs-all structural clustering on HPC |
| **InterProScan** | Not used (CATH only) | Full proteome InterPro as complementary layer on HPC |
| **FoldX** | Deferred | Run on HPC with license |
| **TargetP/SignalP** | Web server for ~266 proteins | Command-line on HPC |
| **Statistical power** | Sufficient for primary hypotheses (N > 200 per group) | Enables exploratory sub-group analyses |
| **Disk usage estimate** | ~3-5 GB | ~15-20 GB |
| **Time estimate** | 4-6 weeks | Additional 4-6 weeks |

### Phase 1 success criteria (before proceeding to Phase 2)

1. Domain architecture features computed for all pilot proteins with > 90% coverage.
2. At least 3 stability metrics computed per protein (pLDDT, contact order, packing density).
3. Ortholog pairs identified with concordant OrthoFinder + MMseqs2 results.
4. Primary hypotheses H1-H4 tested with pilot data. If effect sizes are negligible (Cohen's d < 0.2), reassess before scaling up.
5. MTS features extracted for > 90% of HSP60 Tier 1 substrates.

---

## Risk Register

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| AlphaFold structures missing for some proteins | Reduces coverage | Low (~1-3% expected) | Use UniProt accession mapping. For truly missing structures, exclude from structural analyses but include in sequence-only analyses. |
| CATH coverage < 50% for pilot proteins | Undermines domain analysis | Medium | Chainsaw/Merizo fallback is designed for this. Also consider Gene3D (HMM-based CATH) which has broader coverage than curated CATH. |
| OrthoFinder exceeds local RAM | Blocks orthology on local machine | Low-Medium | Fall back to MMseqs2 RBH only for pilot. Defer full OrthoFinder to HPC. |
| MitoFates server offline | Lose one MTS predictor | High | Already planned as optional. TargetP 2.0 is sufficient as primary predictor. |
| pLDDT misinterpreted as stability | Incorrect conclusions | N/A (mitigated by design) | Decision 2 explicitly separates pLDDT from stability. Multiple orthogonal metrics required. |
| MitoCarta 2.0 vs 3.0 discrepancies affect HSP60 substrate count | Changes in matrix assignment | Already quantified | 70 changes documented. All analysis uses 3.0. Sensitivity analysis with 2.0 annotations as supplementary check. |
| Small effect sizes make results inconclusive | Negative result | Medium | Pre-registered hypotheses prevent post-hoc fishing. Negative results are publishable if properly powered. Power analysis before Phase 2 expansion. |
| Disk space constraint (~18 GB free) | Cannot store all structures | Low-Medium | Monitor usage. AlphaFold PDBs compress well. Delete intermediate files after processing. |

---

## Appendix: Key References

| Reference | Use in project |
|-----------|---------------|
| Kerner MJ et al. (2005) Cell 122:209-220 | Source of GroEL substrate list (Dataset 4) |
| Bruderer R et al. (2020) Cell Stress Chaperones 25:1073-1085 | Source of HSP60 interactome (Dataset 5) |
| Rath S et al. (2021) Nucleic Acids Res 49:D1541-D1547 | MitoCarta 3.0 (Datasets 3, 7) |
| Sillitoe I et al. (2021) Nucleic Acids Res 49:D266-D271 | CATH domain database (Module E) |
| Jumper J et al. (2021) Nature 596:583-589 | AlphaFold (Module D) |
| Emms DM & Kelly S (2019) Genome Biol 20:238 | OrthoFinder (Module C) |
| Steinegger M & Soding J (2017) Nat Biotechnol 35:1026-1028 | MMseqs2 (Module C) |
| Armenteros JJA et al. (2019) Nat Biotechnol 37:420-423 | TargetP 2.0 (Module G) |
| Schymkowitz J et al. (2005) Nucleic Acids Res 33:W382-W388 | FoldX (Module F, Phase 2) |
| Wells JD et al. (2024) Bioinformatics 40:btae187 | Chainsaw domain parser (Module E) |
| Lau AM et al. (2023) Nat Commun 14:7370 | Merizo domain parser (Module E) |

---

*This document is the master reference for the Antah Asti Prarambh project. All implementation decisions should be consistent with the choices recorded here. Update this document when decisions change.*
