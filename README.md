# Antah Asti Prarambh

**"The End is the Beginning"** — Comparative structural proteomics of chaperonin substrates

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

## Project Structure

```
Antah_Asti_Prarambh/
├── scripts/                          # Phase 0: data cleaning and curation
│   ├── validate_uniprot_accessions.py
│   ├── filter_hsp60_interactome.py
│   ├── module_c_extract_fasta.py
│   └── module_c_analyze_rbh.py
├── workflow/
│   ├── scripts/                      # Phase 1: pilot analysis (16 scripts)
│   │   ├── download_alphafold_pilot.py
│   │   ├── run_dssp.py
│   │   ├── get_cath_domains.py
│   │   ├── run_chainsaw_e2.py
│   │   ├── module_f_n_vs_c_analysis.py
│   │   ├── module_g_mts_analysis.py
│   │   ├── module_h_comparative_stats.py
│   │   ├── generate_figures.py
│   │   └── ...
│   └── phase2/                       # Phase 2: full-scale HPC pipeline
│       ├── Snakefile                 # 10-rule dependency graph
│       ├── config.yaml               # Central configuration
│       ├── scripts/                   # Full-scale analysis modules (F, H, I)
│       └── slurm_jobs/               # 19 SLURM job scripts
├── data/
│   ├── raw/                          # Source data (custom tables, UniProt metadata)
│   └── processed/                    # Curated datasets (7 TSV/FASTA files)
├── results/
│   ├── domains/                      # CATH, Chainsaw, Foldseek, unified assignments
│   ├── homology/                     # RBH, orthogroups, substrate pairs
│   ├── structures/                   # DSSP summaries, quality validation
│   ├── termini/                      # N-vs-C paired metrics, contact order
│   ├── mts/                          # Targeting classification, MTS-domain overlap
│   ├── stats/                        # Hypothesis tests, corrected p-values
│   ├── figures/                      # Phase 1 publication figures (6 fig x PDF+PNG)
│   └── phase2/                       # Full-scale results (25,007 proteins)
└── docs/                             # Documentation, methods, session logs
```

## Reproducing the Analysis

### Phase 1 (local, ~1,390 pilot proteins)

```bash
# 1. Set up conda environment
conda create -n proteomics python=3.9
conda activate proteomics
pip install biopython pandas scipy matplotlib seaborn openpyxl
conda install -c bioconda mmseqs2 foldseek

# 2. Run modules sequentially
python scripts/validate_uniprot_accessions.py
python scripts/filter_hsp60_interactome.py
python workflow/scripts/download_alphafold_pilot.py
python workflow/scripts/run_dssp.py
python workflow/scripts/get_cath_domains.py
python workflow/scripts/run_chainsaw_e2.py
# ... (see docs/METHODS_AND_PROTOCOLS.md for full sequence)
```

### Phase 2 (HPC, ~25,000 proteins)

```bash
# Deploy to HPC, edit config.yaml for local paths
# Submit via SLURM pipeline:
bash workflow/phase2/slurm_jobs/submit_pipeline.sh
# Or use Snakemake:
snakemake --profile slurm -j 20
```

See [`docs/HPC_DEPLOYMENT_GUIDE.md`](docs/HPC_DEPLOYMENT_GUIDE.md) for full HPC instructions.

## Data Availability

- **AlphaFold structures**: [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) (not included due to size; ~22 GB for full proteomes)
- **MitoCarta 3.0**: [Broad Institute](https://www.broadinstitute.org/mitocarta)
- **GroEL substrates**: Kerner et al. 2005, *Cell* 122(2):209-220 — Table S3
- **HSP60 interactome**: Bruderer et al. 2020, *Mol Cell Proteomics* — Supplementary Table 4

## Documentation

| Document | Description |
|----------|-------------|
| [`docs/DOCUMENTATION.md`](docs/DOCUMENTATION.md) | Master project documentation (~1,255 lines) |
| [`docs/METHODS_AND_PROTOCOLS.md`](docs/METHODS_AND_PROTOCOLS.md) | Reproducibility guide with exact commands |
| [`docs/RESULTS_NARRATIVE.md`](docs/RESULTS_NARRATIVE.md) | Manuscript-style results with all statistics |
| [`docs/PRIMARY_HYPOTHESES.md`](docs/PRIMARY_HYPOTHESES.md) | 9 pre-registered hypotheses |
| [`docs/PHASE1_VERIFICATION.md`](docs/PHASE1_VERIFICATION.md) | Phase 1 success criteria (all 5 PASS) |
| [`docs/PHASE2_RESULTS_REPORT.md`](docs/PHASE2_RESULTS_REPORT.md) | Full-scale Phase 2 results |
| [`docs/COMPREHENSIVE_PROJECT_DOCUMENT.md`](docs/COMPREHENSIVE_PROJECT_DOCUMENT.md) | Complete reference (~1,618 lines) |

## Status

- **Phase 1 (pilot)**: Complete and verified
- **Phase 2 (full-scale)**: Complete except FoldX thermodynamic stability (~48% done on HPC)
- **Next**: Integrate FoldX DeltaG results, regenerate figures, prepare manuscript

## Citation

*Manuscript in preparation.*

## License

This project is shared for academic and research purposes. Please contact the authors before reuse.

## Authors

Vishal Bharti — CSIR-Institute of Genomics and Integrative Biology (IGIB), New Delhi
