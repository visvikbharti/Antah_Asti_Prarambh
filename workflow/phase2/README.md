# Phase 2: Full-Scale HPC Analysis

## Overview

Phase 2 scales the Antah Asti Prarambh pilot analysis (Phase 1: ~900 proteins) to the
full E. coli K-12 (4,403 proteins) and Human (20,416 proteins) proteomes. These analyses
require HPC resources due to memory (32-64 GB), compute time, and disk (~22 GB for
structures alone) requirements that exceed the local machine (8 GB RAM, 15 GB disk).

## Prerequisites

1. **Phase 1 complete** -- all pilot analyses validated (see `docs/TODO.md`)
2. **HPC account** with SLURM scheduler
3. **Software installed on HPC:**
   - Python 3.9+ with: pandas, numpy, scipy, statsmodels, pyyaml, requests, matplotlib, seaborn
   - Foldseek (v10+): `conda install -c bioconda foldseek`
   - Chainsaw: clone from `https://github.com/JudeWells/chainsaw`
   - FoldX 5.0 (requires license from https://foldxsuite.crg.eu/)
   - gemmi (for CIF-to-PDB conversion): `pip install gemmi`
4. **Disk space:** ~50 GB on scratch filesystem (structures + temp databases)

## Directory Structure

```
workflow/phase2/
  config.yaml                  # All paths, parameters, SLURM settings
  Snakefile                    # Snakemake workflow (dependency graph)
  download_alphafold_full.py   # Download full AlphaFold proteomes
  run_foldseek_full.py         # Full-scale Foldseek clustering
  run_foldx.py                 # FoldX stability calculations
  README.md                    # This file
```

## Quick Start

### 1. Transfer Phase 1 data to HPC

```bash
# On local machine:
rsync -avz /Users/vishalbharti/Downloads/Antah_Asti_Prarambh/ \
  user@hpc:/scratch/$USER/antah_asti_prarambh/phase1/
```

### 2. Edit config.yaml

Update paths in `config.yaml` to match your HPC filesystem:

```yaml
project_dir: "/scratch/$USER/antah_asti_prarambh"
phase1_dir: "/scratch/$USER/antah_asti_prarambh/phase1"
```

Update SLURM settings:

```yaml
slurm:
  partition: "your_partition"
  account: "your_allocation"
```

Update tool paths:

```yaml
foldseek:
  binary: "/path/to/foldseek"  # or just "foldseek" if in PATH

foldx:
  binary: "/path/to/foldx"
  rotabase: "/path/to/rotabase.txt"

chainsaw:
  install_dir: "/path/to/chainsaw"
```

### 3. Run with Snakemake (recommended)

```bash
# Load modules (example)
module load python/3.9
module load conda
conda activate proteomics

# Navigate to workflow directory
cd /scratch/$USER/antah_asti_prarambh/phase1/workflow/phase2

# Dry run -- see what will execute
snakemake -n --configfile config.yaml

# Run on SLURM cluster
snakemake --configfile config.yaml \
  --cluster "sbatch --partition={resources.partition} \
             --mem={resources.mem_mb}M \
             --cpus-per-task={threads} \
             --time={resources.runtime} \
             --output=logs/slurm_%j.out" \
  --jobs 50 \
  --latency-wait 60 \
  --rerun-incomplete
```

### 4. Run individual scripts manually (alternative)

If you prefer manual control or need to debug:

```bash
# Step 1: Download AlphaFold structures (~22 GB, 1-4 hours)
python download_alphafold_full.py --config config.yaml --organism both

# Step 2: Foldseek clustering (64 GB RAM, 16 CPUs, 6-24 hours)
python run_foldseek_full.py --config config.yaml --threads 16 --memory 64G

# Step 3: FoldX stability calculations (submit as array job)
python run_foldx.py --config config.yaml --generate-slurm
sbatch foldx_array.slurm
# After array jobs complete:
python run_foldx.py --config config.yaml --collect
```

## Resource Requirements

| Task | RAM | CPUs | Time | Disk |
|------|-----|------|------|------|
| AlphaFold download | 4 GB | 1 | 1-4 hrs | 22 GB |
| Foldseek createdb | 16 GB | 4 | 30 min | 5 GB |
| Foldseek search | 64 GB | 16 | 6-24 hrs | 30 GB temp |
| Foldseek cluster | 64 GB | 16 | 2-4 hrs | 5 GB |
| Chainsaw domains | 16 GB | 4 | 4-8 hrs | 1 GB |
| FoldX (per protein) | 2 GB | 1 | 1-3 min | <1 GB |
| FoldX (all, serial) | 4 GB | 1 | ~48 hrs | 5 GB |
| FoldX (array, 50/job) | 4 GB | 1 | 4 hrs/job | 5 GB |
| Module E-H analysis | 16 GB | 4 | 1-2 hrs | <1 GB |
| Figure generation | 8 GB | 2 | 30 min | <1 GB |

**Total disk estimate:** ~50 GB (including temporary files)

## Pipeline Dependency Graph

```
download_alphafold
       |
 +-----+-----+-------------------+
 |           |                   |
 v           v                   v
foldseek  chainsaw_domains    foldx_stability
 |           |                   |
 v           v                   v
foldseek  domain_metrics      stability_metrics
_analyze     |                   |
 |           +---+---+-----------+
 |               |   |
 v               v   v
module_e    module_f_analysis
 |               |
 +---+---+---+---+
     |       |
     v       v
module_g  module_h_stats
             |
             v
        module_i_figures
```

## Resume and Checkpointing

All scripts support resume:

- **download_alphafold_full.py**: JSON checkpoint file tracks bulk/individual downloads.
  Re-running skips completed downloads automatically.
- **run_foldseek_full.py**: Checks for existing database files and search results.
  Re-running skips completed steps.
- **run_foldx.py**: Per-protein JSON results in `per_protein/` directory.
  Re-running skips proteins with existing results.
- **Snakemake**: Automatically tracks which rules have completed outputs.

## FoldX Array Jobs (SLURM)

For 25,000 proteins at 50 proteins/job = 500 array tasks:

```bash
# Generate the SLURM script
python run_foldx.py --config config.yaml --generate-slurm

# Submit array job
JOB_ID=$(sbatch foldx_array.slurm | awk '{print $4}')
echo "Submitted array job: $JOB_ID"

# Submit collection job (runs after array completes)
sbatch --dependency=afterok:$JOB_ID collect_results.slurm

# Monitor progress
squeue -u $USER
sacct -j $JOB_ID --format=JobID,State,Elapsed,MaxRSS
```

## Troubleshooting

### Foldseek out of memory

Increase `--split-memory-limit` in config.yaml or pass `--memory 128G`.
The search step is the most memory-intensive.

### FoldX "rotabase not found"

Set the environment variable:
```bash
export FOLDX_ROTABASE=/path/to/rotabase.txt
```

Or update `foldx.rotabase` in config.yaml.

### Chainsaw import errors

Ensure Chainsaw dependencies are installed:
```bash
pip install torch einops
```

### Download failures

The download script supports resume. Simply re-run the same command.
For persistent failures, check:
- Network connectivity to EBI FTP
- Disk space on scratch filesystem
- File permissions

### Snakemake "incomplete files"

```bash
snakemake --configfile config.yaml --rerun-incomplete
```

## Output Files

After completion, key Phase 2 outputs will be in:

```
results/phase2/
  structures/
    structure_index_full.tsv       # Index of all 25K structures
    download_checkpoint.json       # Download progress
  foldseek/
    analysis/
      foldseek_clusters_full.tsv   # Cluster assignments for all proteins
      foldseek_full_summary.txt    # Clustering summary report
  domains/
    chainsaw_full_predictions.tsv  # Domain boundaries for all proteins
    domain_distribution_full.tsv   # Domain architecture statistics
  foldx/
    foldx_stability_all.tsv        # DeltaG for all proteins
  stats/
    stability_comparisons_full.tsv # N-vs-C comparison with FoldX
    statistics_summary_full.txt    # Full statistical analysis report
  figures/
    fig1-fig7 (PDF + PNG at 300 DPI)
```
