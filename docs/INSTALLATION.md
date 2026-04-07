# Installation Guide

Complete guide for setting up all software dependencies for the Antah Asti Prarambh pipeline.

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/visvikbharti/Antah_Asti_Prarambh.git
cd Antah_Asti_Prarambh

# 2. Create conda environment (includes Python + all bioinformatics tools)
conda env create -f environment.yml
conda activate proteomics

# 3. Download external datasets
bash scripts/download_external_data.sh

# 4. Run the analysis pipeline
make analysis
```

## Prerequisites

- **Operating system**: Linux or macOS (Apple Silicon works via Rosetta 2)
- **Conda**: Miniconda or Anaconda ([install guide](https://docs.conda.io/en/latest/miniconda.html))
- **Disk space**: ~1 GB for core analysis, ~50 GB for full-scale HPC pipeline (AlphaFold structures)
- **RAM**: 8 GB minimum, 16 GB recommended

## Core Environment (conda)

The `environment.yml` file defines the full conda environment:

```bash
conda env create -f environment.yml
conda activate proteomics
```

This installs:

| Package | Version | Purpose |
|---------|---------|---------|
| Python | 3.9 | Runtime |
| pandas | 2.2.2 | Data manipulation |
| numpy | latest | Numerical computing |
| scipy | 1.9.2 | Statistical tests |
| matplotlib | latest | Plotting |
| seaborn | latest | Statistical visualization |
| biopython | 1.78 | Sequence/structure parsing |
| mmseqs2 | 18.8cc5c | Reciprocal best hits, orthology |
| foldseek | 10.941cd33 | Structural clustering |
| dssp | latest | Secondary structure (mkdssp) |
| gemmi | latest | CIF/PDB file handling |
| snakemake | 7.32.4 | Workflow management |
| statsmodels | latest | Multiple testing correction |
| python-pptx | latest | Presentation generation |

### Updating the environment

```bash
conda env update -f environment.yml --prune
```

### Apple Silicon (M1/M2/M3) notes

MMseqs2 and Foldseek are x86_64 binaries on bioconda. They run via Rosetta 2:

```bash
# Rosetta 2 should auto-install on first use; if not:
softwareupdate --install-rosetta
```

## External Tools

These tools are **not installable via conda** and require separate setup.

### Chainsaw (ML domain segmentation)

Required for: Module E2 domain segmentation, full-scale domain assignment

```bash
# Install from GitHub
git clone https://github.com/JudeWells/chainsaw.git
cd chainsaw
pip install -e .

# Verify
python -c "from chainsaw import get_predictions; print('Chainsaw OK')"
```

**HPC deployment**: Install in your software directory and update `config.yaml`:
```yaml
chainsaw:
  install_dir: "/path/to/chainsaw"
```

### FoldX (thermodynamic stability)

Required for: Full-scale DeltaG stability calculations (HPC)

FoldX requires an **academic license** (free for non-commercial use):

1. Register at [https://foldxsuite.crg.eu/](https://foldxsuite.crg.eu/)
2. Download FoldX 5.x for Linux
3. Extract and note the binary path

```bash
# Extract
unzip foldx5_Linux.zip -d /path/to/foldx5/

# Verify
/path/to/foldx5/foldx --version
# Expected: FoldX 5.1 or similar

# Update config.yaml
# foldx:
#   binary: "/path/to/foldx5/foldx"
```

**Notes:**
- FoldX 5.1 does NOT require `rotabase.txt`
- FoldX was parameterized on experimental structures, not AlphaFold models — caveat this in publications
- License is typically valid for 1 year

### DSSP (secondary structure assignment)

Usually installed via conda (included in `environment.yml`). If you need it standalone:

```bash
# Via conda (recommended)
conda install -c conda-forge dssp

# Verify
mkdssp --version
```

### TargetP 2.0 and SignalP 6.0 (optional)

These DTU tools are used for transit peptide prediction but are **not required** — the pipeline uses UniProt transit peptide annotations instead.

If you want to run them:
1. Register at [https://services.healthtech.dtu.dk/](https://services.healthtech.dtu.dk/)
2. Download TargetP 2.0 and/or SignalP 6.0
3. Follow the DTU installation instructions

### IUPred2a (optional)

For intrinsic disorder prediction. Not required for the core pipeline.

1. Register at [https://iupred2a.elte.hu/](https://iupred2a.elte.hu/)
2. Download and follow installation instructions

## HPC Pipeline Setup

For running the full-scale pipeline (25,007 proteins) on an HPC cluster with SLURM:

### 1. Transfer project to HPC

```bash
rsync -avz /path/to/Antah_Asti_Prarambh/ \
    user@hpc:/path/to/Antah_Asti_Prarambh_hpc/
```

### 2. Configure paths

```bash
cd /path/to/Antah_Asti_Prarambh_hpc/workflow/phase2/
cp config.example.yaml config.yaml
# Edit config.yaml: set project_dir, phase1_dir, chainsaw.install_dir, foldx.binary
```

### 3. Run setup job

```bash
# Creates directories, installs conda env, validates files
sbatch workflow/phase2/slurm_jobs/00_setup.sh
```

### 4. Submit pipeline

```bash
# Submit all jobs with dependencies
bash workflow/phase2/slurm_jobs/submit_pipeline.sh
```

See [`docs/HPC_DEPLOYMENT_GUIDE.md`](HPC_DEPLOYMENT_GUIDE.md) for detailed HPC instructions.

## Verification

After installation, verify everything works:

```bash
# Activate environment
conda activate proteomics

# Check core tools
python3 -c "import pandas, scipy, matplotlib, Bio, statsmodels; print('Python packages OK')"
mmseqs version
foldseek version
mkdssp --version

# Check optional tools
python3 -c "import gemmi; print('gemmi OK')"
which chainsaw 2>/dev/null || echo "Chainsaw: install separately (see above)"

# Run validation
make validate
```

## Troubleshooting

### `libarchive.so.20: cannot open shared object file`
This is a conda-libmamba-solver issue. Safe to ignore — it doesn't affect package functionality.

### MMseqs2/Foldseek segfault on Apple Silicon
Ensure Rosetta 2 is installed: `softwareupdate --install-rosetta`

### `mkdssp` not found
Install via conda: `conda install -c conda-forge dssp`

### CATH/Gene3D API timeout
The `get_cath_domains.py` script uses the InterPro API with checkpointing. If it times out, re-run — it will resume from the last checkpoint.

### FoldX "License expired"
FoldX academic licenses are valid for ~1 year. Re-register at [foldxsuite.crg.eu](https://foldxsuite.crg.eu/) for a new license.

### SLURM `QOSMaxJobsPerUserLimit`
This means your HPC QOS limits concurrent jobs. The pipeline handles this — jobs queue and run as slots open. No action needed.

### `LD_LIBRARY_PATH` issues on HPC
Some HPC systems need this for gcc/numpy compatibility. Add to SLURM scripts:
```bash
export LD_LIBRARY_PATH=/path/to/gcc-8.4/lib64:$LD_LIBRARY_PATH
```
