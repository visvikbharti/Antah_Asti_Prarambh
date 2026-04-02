# Antah Asti Prarambh

**"The End is the Beginning"** — Comparative structural proteomics of chaperonin substrates

[![Phase 1](https://img.shields.io/badge/Phase%201-Complete-brightgreen)]()
[![Phase 2](https://img.shields.io/badge/Phase%202-Complete-brightgreen)]()
[![Python](https://img.shields.io/badge/Python-3.9-blue)]()

## Overview

This project investigates the structural properties of chaperonin substrates across two organisms: **GroEL** in *Escherichia coli* (Kerner et al. 2005) and **HSP60/HSPD1** in human mitochondria (Bruderer et al. 2020). Using AlphaFold-predicted structures, CATH/Chainsaw domain assignments, and curated orthology, we analyze whether chaperonin substrates share conserved structural features and how mitochondrial targeting peptides relate to the first structural domain.

### Three Scientific Goals

| # | Goal | Question |
|---|------|----------|
| 1 | **Domain architecture** | Do chaperonin substrates have distinct domain folds compared to non-substrate backgrounds? Are patterns conserved between GroEL and HSP60? |
| 2 | **N-vs-C terminus stability** | Are N-terminal structural domains less stable than C-terminal regions? Is this asymmetry substrate-specific or universal? |
| 3 | **Mitochondrial targeting** | Is the MTS a separate pre-domain extension or does it overlap the first structural domain? |

### Key Findings

- **N-terminal contact order asymmetry is universal** — N-domains have higher relative contact order than C-regions across all multi-domain proteins (p < 1e-20), not just chaperonin substrates. This suggests chaperonins exploit a pre-existing structural property rather than driving substrate evolution.
- **GroEL substrates enriched in TIM barrels** (OR = 8.4) and Rossmann-like folds.
- **84.4% of mitochondrial transit peptides** are separate pre-domain extensions (median gap: 18 residues from cleavage site to first domain).
- **HSP60 substrates 3.3x enriched** for mitochondrial matrix localization.
- **Homolog pair conservation**: r = 0.82 for N-domain contact order between GroEL-HSP60 ortholog pairs.

---

## Pipeline Flowchart

```
                        ┌─────────────────────────────────────┐
                        │        PHASE 0: DATA CURATION       │
                        └─────────────────┬───────────────────┘
                                          │
              ┌───────────────────────────┼───────────────────────────┐
              ▼                           ▼                           ▼
    ┌─────────────────┐       ┌─────────────────┐       ┌─────────────────────┐
    │  GroEL (252)    │       │  HSP60 (266)    │       │  MitoCarta 3.0      │
    │  Kerner 2005    │       │  Bruderer 2020  │       │  1,136 mito / 525   │
    │  ID remapping   │       │  SILAC filter   │       │  matrix proteins    │
    └────────┬────────┘       └────────┬────────┘       └──────────┬──────────┘
             │                         │                           │
             └─────────────┬───────────┘                           │
                           ▼                                       │
              ┌────────────────────────┐                           │
              │  MODULE C: ORTHOLOGY   │                           │
              │  MMseqs2 RBH (40 pairs)│                           │
              │  Orthogroups (422)     │                           │
              │  Dataset 6 (69 pairs)  │                           │
              └────────────┬───────────┘                           │
                           │                                       │
                           ▼                                       │
              ┌────────────────────────┐                           │
              │  MODULE D: STRUCTURES  │                           │
              │  AlphaFold (1,382 CIF) │                           │
              │  DSSP (sec. structure) │                           │
              │  Quality validation    │                           │
              └────────────┬───────────┘                           │
                           │                                       │
         ┌─────────────────┼─────────────────┐                     │
         ▼                 ▼                 ▼                     │
┌──────────────┐  ┌──────────────┐  ┌──────────────┐              │
│ MODULE E:    │  │ Foldseek     │  │ FoldX 5.1    │              │
│ CATH/Gene3D  │  │ structural   │  │ thermo-      │              │
│ + Chainsaw   │  │ clustering   │  │ dynamic ΔG   │              │
│ (99.8% cov)  │  │ (1,155 clust)│  │ (Phase 2)    │              │
└──────┬───────┘  └──────┬───────┘  └──────┬───────┘              │
       │                 │                 │                       │
       └─────────────────┼─────────────────┘                       │
                         ▼                                         │
              ┌────────────────────────┐                           │
              │  MODULE F: N-vs-C      │                           │
              │  STABILITY ANALYSIS    │                           │
              │  3-region decomposition│                           │
              │  Contact order, pLDDT  │                           │
              │  Paired comparisons    │◄──────────────────────────┘
              └────────────┬───────────┘
                           │
              ┌────────────┼────────────┐
              ▼                         ▼
┌────────────────────────┐  ┌────────────────────────┐
│  MODULE G: MTS         │  │  MODULE H: STATISTICS  │
│  TARGETING ANALYSIS    │  │  9 hypotheses          │
│  Transit peptides      │  │  3 families            │
│  MTS-domain overlap    │  │  Hierarchical BH       │
│  Matrix enrichment     │  │  Effect sizes          │
└────────────┬───────────┘  └────────────┬───────────┘
             │                           │
             └─────────────┬─────────────┘
                           ▼
              ┌────────────────────────┐
              │  MODULE I: FIGURES     │
              │  6 publication figures │
              │  PDF + PNG (300 DPI)   │
              └────────────────────────┘
```

### Phase 2 HPC Dependency Graph

```
  download_alphafold (22 GB)
          │
    ┌─────┼─────────────────┐
    │     │                 │
    ▼     ▼                 ▼
 Foldseek  Chainsaw       FoldX (501 array jobs)
    │     │                 │
    ▼     ▼                 ▼
 Analyze  Unified         Collect ΔG
    │     domains           │
    └──┬──┘     ┌───────────┘
       │        │
       ▼        ▼
    Module F (stability)
       │
       ▼
    Module H (stats) ──► Module I (figures)
```

---

## Quick Start

```bash
# Clone
git clone https://github.com/visvikbharti/Antah_Asti_Prarambh.git
cd Antah_Asti_Prarambh

# Setup (conda environment + external data)
make setup

# Run Phase 1
make phase1

# Check results
make validate
```

See [`docs/INSTALLATION.md`](docs/INSTALLATION.md) for detailed setup, including external tools (Chainsaw, FoldX).

---

## Datasets

| # | Dataset | N | Source |
|---|---------|--:|--------|
| 1 | *E. coli* K-12 proteome | 4,403 | UniProt UP000000625 |
| 2 | Human reference proteome | 20,416 | UniProt UP000005640 |
| 3 | Human mito proteome | 1,136 | MitoCarta 3.0 |
| 4 | GroEL substrates | 252 | Kerner et al. 2005 |
| 5 | HSP60 interactome (Tier 1) | 266 | Bruderer et al. 2020 |
| 6 | Homolog pairs | 69 | RBH + orthogroup intersection |
| 7 | Mito matrix-only | 525 | MitoCarta 3.0 sub-localization |

---

## Methods

### Tools and Software

| Tool | Version | Purpose |
|------|---------|---------|
| AlphaFold DB | v4/v6 | Predicted protein structures |
| DSSP (mkdssp) | 2.2.1 | Secondary structure assignment |
| CATH/Gene3D | via InterPro API | Structural domain boundaries |
| Chainsaw | v3 | ML-based domain segmentation (gap-fill for CATH) |
| MMseqs2 | 18.8cc5c | Reciprocal best hits + orthology |
| Foldseek | 10.941cd33 | Structural clustering |
| FoldX | 5.1 | Thermodynamic stability (DeltaG) |
| Contact order | Plaxco 1998 | Folding complexity (8A cutoff, min_sep=6) |

### Statistical Framework

- **Hierarchical BH correction** across 3 goal families (domain architecture, stability asymmetry, matrix targeting)
- **9 pre-registered hypotheses** with family-level Simes gates
- **Compartment-matched AND size-matched controls** (1,000 permutations, +/-20% length bins)
- **Effect sizes**: Cohen's d for paired comparisons, odds ratios for enrichment

---

## Project Structure

```
Antah_Asti_Prarambh/
├── Makefile                             # Master build system (make help)
├── environment.yml                      # Conda environment definition
├── requirements.txt                     # Pip dependencies
├── scripts/                             # Phase 0: data cleaning and curation
│   ├── download_external_data.sh        # Download all external datasets
│   ├── validate_uniprot_accessions.py   # GroEL ID remapping
│   ├── filter_hsp60_interactome.py      # HSP60 SILAC filtering
│   ├── module_c_extract_fasta.py        # FASTA extraction
│   └── module_c_analyze_rbh.py          # RBH analysis
├── workflow/
│   ├── scripts/                         # Phase 1: pilot analysis (16 scripts)
│   │   ├── download_alphafold_pilot.py  # AlphaFold structure download
│   │   ├── run_dssp.py                  # Secondary structure
│   │   ├── get_cath_domains.py          # CATH/Gene3D API
│   │   ├── run_chainsaw_e2.py           # ML domain segmentation
│   │   ├── module_f_n_vs_c_analysis.py  # N-vs-C paired analysis
│   │   ├── module_g_mts_analysis.py     # MTS targeting
│   │   ├── module_h_comparative_stats.py # Hierarchical statistics
│   │   ├── generate_figures.py          # Publication figures
│   │   └── ...
│   └── phase2/                          # Phase 2: full-scale HPC pipeline
│       ├── Snakefile                    # Workflow dependency graph
│       ├── config.yaml                  # HPC configuration
│       ├── config.example.yaml          # Template (edit for your system)
│       ├── scripts/                     # Full-scale modules (F, H, I)
│       └── slurm_jobs/                  # 19 SLURM job scripts
├── data/
│   ├── raw/                             # Source data (downloaded or provided)
│   │   ├── custom/                      # Kerner 2005 + Bruderer 2020 tables
│   │   ├── uniprot/                     # Proteome FASTA + metadata
│   │   └── mitocarta/                   # MitoCarta 3.0
│   └── processed/                       # Curated datasets (7 files)
├── results/
│   ├── domains/                         # CATH, Chainsaw, Foldseek, unified
│   ├── homology/                        # RBH, orthogroups, substrate pairs
│   ├── structures/                      # DSSP summaries, quality validation
│   ├── termini/                         # N-vs-C paired metrics, contact order
│   ├── mts/                             # Targeting, MTS-domain overlap
│   ├── stats/                           # Hypothesis tests, corrected p-values
│   ├── figures/                         # Phase 1 figures (6 x PDF+PNG)
│   └── phase2/                          # Full-scale results (25,007 proteins)
│       ├── domains/                     # Chainsaw full-scale, unified
│       ├── foldseek/                    # 16,193 clusters
│       ├── stability/                   # N-vs-C with 2,648 pairs
│       ├── stats/                       # 56 tests, 25 significant
│       └── figures/                     # Phase 2 figures (6 x PDF+PNG)
└── docs/                                # Documentation
```

---

## Reproducing the Analysis

### Prerequisites

- Conda (Miniconda or Anaconda)
- 8 GB RAM (16 GB recommended)
- ~1 GB disk for Phase 1; ~50 GB for Phase 2

### Phase 1 (local, ~1,390 pilot proteins)

```bash
# Option A: Using Makefile (recommended)
make setup     # Creates conda env + downloads data
make phase1    # Runs complete pipeline

# Option B: Step by step
conda env create -f environment.yml
conda activate proteomics
bash scripts/download_external_data.sh

make phase0        # Data cleaning
make orthology     # Module C: RBH + orthogroups
make structures    # Module D: AlphaFold + DSSP + quality
make domains       # Module E: CATH + Chainsaw + Foldseek
make stability     # Module F: N-vs-C paired analysis
make targeting     # Module G: MTS targeting
make statistics    # Module H: hierarchical tests
make figures       # Module I: 6 publication figures
```

### Phase 2 (HPC, ~25,000 proteins)

```bash
# 1. Transfer to HPC
rsync -avz Antah_Asti_Prarambh/ user@hpc:/path/to/project/

# 2. Configure
cd workflow/phase2
cp config.example.yaml config.yaml
# Edit config.yaml: set project_dir, chainsaw path, foldx path

# 3. Submit pipeline
sbatch slurm_jobs/00_setup.sh              # Environment + validation
bash slurm_jobs/submit_pipeline.sh          # All jobs with dependencies

# Or via Snakemake:
snakemake --configfile config.yaml --profile slurm -j 50
```

See [`docs/HPC_DEPLOYMENT_GUIDE.md`](docs/HPC_DEPLOYMENT_GUIDE.md) for detailed HPC instructions.

---

## Data Availability

External data is **not included** in the repository due to size. The download script fetches everything automatically:

```bash
bash scripts/download_external_data.sh
```

| Dataset | Source | Size | Auto-download |
|---------|--------|------|:---:|
| *E. coli* K-12 proteome | [UniProt UP000000625](https://www.uniprot.org/proteomes/UP000000625) | 1.8 MB | Yes |
| Human proteome | [UniProt UP000005640](https://www.uniprot.org/proteomes/UP000005640) | 13 MB | Yes |
| MitoCarta 3.0 | [Broad Institute](https://www.broadinstitute.org/mitocarta) | 9.7 MB | Yes |
| GroEL substrates | Kerner et al. 2005, *Cell* 122(2):209-220 | Included | N/A |
| HSP60 interactome | Bruderer et al. 2020, *Mol Cell Proteomics* | Manual | No* |
| AlphaFold structures | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | 22 GB | Phase 2 only |

*\*Bruderer 2020 supplement may require institutional access. See download script for URL.*

---

## Documentation

| Document | Description |
|----------|-------------|
| **Setup & Guides** | |
| [`docs/INSTALLATION.md`](docs/INSTALLATION.md) | Installation guide for all tools |
| [`docs/HPC_DEPLOYMENT_GUIDE.md`](docs/HPC_DEPLOYMENT_GUIDE.md) | HPC deployment with SLURM |
| [`docs/METHODS_AND_PROTOCOLS.md`](docs/METHODS_AND_PROTOCOLS.md) | Reproducibility guide with exact commands |
| **Science** | |
| [`docs/PRIMARY_HYPOTHESES.md`](docs/PRIMARY_HYPOTHESES.md) | 9 pre-registered hypotheses |
| [`docs/RESULTS_NARRATIVE.md`](docs/RESULTS_NARRATIVE.md) | Manuscript-style results with all statistics |
| [`docs/DOCUMENTATION.md`](docs/DOCUMENTATION.md) | Master project documentation (~1,255 lines) |
| [`docs/COMPREHENSIVE_PROJECT_DOCUMENT.md`](docs/COMPREHENSIVE_PROJECT_DOCUMENT.md) | Complete reference (~1,618 lines) |
| **Verification** | |
| [`docs/PHASE1_VERIFICATION.md`](docs/PHASE1_VERIFICATION.md) | Phase 1 success criteria (all 5 PASS) |
| [`docs/PHASE2_RESULTS_REPORT.md`](docs/PHASE2_RESULTS_REPORT.md) | Full-scale Phase 2 results |
| **Collaboration** | |
| [`docs/COLLABORATOR_PRESENTATION.md`](docs/COLLABORATOR_PRESENTATION.md) | Meeting-ready presentation materials |
| [`docs/DATA_HANDOVER_INDEX.md`](docs/DATA_HANDOVER_INDEX.md) | File-by-file result guide |

---

## Key Results Summary

### Goal 1: Domain Architecture

| Metric | Value |
|--------|-------|
| CATH coverage | 82.8% (Gene3D); 99.8% unified (+ Chainsaw) |
| GroEL TIM barrel enrichment | OR = 8.4, p_BH = 2.3e-6 |
| Cross-organism fold conservation | 79.7% of homolog pairs share top superfamily |
| Foldseek clusters | 16,193 (full-scale) |

### Goal 2: N-vs-C Stability

| Metric | Value |
|--------|-------|
| Contact order asymmetry (N > C) | p = 1.05e-20 (universal) |
| Substrate-specific? | **No** — background shows same pattern (all p > 0.14) |
| GroEL class gradient? | **No** — KW p = 0.77 |
| Homolog pair RCO correlation | r = 0.82 |

### Goal 3: Mitochondrial Targeting

| Metric | Value |
|--------|-------|
| HSP60 matrix enrichment | OR = 3.29, p = 1.6e-16 |
| MTS as pre-domain extension | 84.4% (binomial p = 3.4e-51) |
| Median MTS-to-domain gap | 18 residues |
| Non-canonical matrix import | 21.1% of HSP60 substrates |

---

## Important Scientific Notes

1. **pLDDT is NOT thermodynamic stability.** It is AlphaFold's per-residue model confidence. FoldX DeltaG is the correct stability metric.

2. **The N-vs-C asymmetry is universal (key negative result).** Higher N-terminal contact order is a general property of multi-domain proteins, not a chaperonin-specific adaptation. Chaperonins may exploit a pre-existing structural asymmetry.

3. **FoldX was parameterized on experimental structures, not AlphaFold models.** This is a known limitation that should be caveated in publications.

4. **NDIC in HSP60 data = "Not Detected In Control" = very high enrichment**, not missing data.

---

## Status

- **Phase 1 (pilot)**: Complete and verified (all 5 success criteria PASS)
- **Phase 2 (full-scale)**: **COMPLETE** including FoldX thermodynamic stability (25,007 proteins, 0 failures)
- **FoldX key finding**: GroEL substrates have slightly lower total energy (median -38.6 vs -15.2 for *E. coli* bg, p=2.9e-3, d=-0.07; compartment-matched)
- **Statistics**: 60 tests across 3 families, 28 significant after hierarchical BH correction
- **Next**: Manuscript preparation

---

## Citation
And the key data sources:
- Kerner et al. (2005) Proteome-wide analysis of chaperonin-dependent protein folding in *Escherichia coli*. *Cell* 122(2):209-220.
- Bruderer et al. (2020) The HSPD1 interactome. *Mol Cell Proteomics*.
- Jumper et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583-589.

## License

This project is shared for academic and research purposes. Please contact the authors before reuse.

## Authors

Vishal Bharti — CSIR-Institute of Genomics and Integrative Biology (IGIB), New Delhi
