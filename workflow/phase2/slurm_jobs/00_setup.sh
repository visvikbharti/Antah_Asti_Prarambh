#!/bin/bash
#SBATCH --job-name=aap_setup
#SBATCH --output=logs/00_setup_%j.out
#SBATCH --error=logs/00_setup_%j.err
#SBATCH --partition=compute
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

# =============================================================================
# Antah Asti Prarambh — Phase 2 Setup
# Creates conda environment, directory structure, validates data
# =============================================================================

set -euo pipefail

# ---- PATHS (hardcoded, no placeholders) ------------------------------------
export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
CONDA_BASE="/home/vishal.bharti/miniconda3"
# ----------------------------------------------------------------------------

echo "============================================================"
echo "  Phase 2 Setup — $(date)"
echo "  Project: ${PROJECT_DIR}"
echo "============================================================"

# Initialize conda
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Create proteomics conda env if it doesn't exist
if ! conda env list | grep -q "^${CONDA_ENV} "; then
    echo ""
    echo "Creating conda environment: ${CONDA_ENV}"
    conda create -n "${CONDA_ENV}" python=3.11 -y
    conda activate "${CONDA_ENV}"

    # Install packages from conda-forge and bioconda
    conda install -c conda-forge -c bioconda \
        pandas numpy scipy matplotlib seaborn biopython pyyaml \
        openpyxl requests tqdm gemmi foldseek -y

    # Install additional pip packages
    pip install statsmodels einops

    echo "Conda environment ${CONDA_ENV} created successfully."
else
    echo "Conda environment ${CONDA_ENV} already exists."
    conda activate "${CONDA_ENV}"
fi

echo ""
echo "Python: $(which python3) ($(python3 --version))"
echo "Conda env: ${CONDA_ENV}"

# Create directory structure (all on /lustre)
echo ""
echo "Creating directory structure..."
mkdir -p "${PROJECT_DIR}/data/raw/alphafold/full/ecoli"
mkdir -p "${PROJECT_DIR}/data/raw/alphafold/full/human"
mkdir -p "${PROJECT_DIR}/results/phase2/structures"
mkdir -p "${PROJECT_DIR}/results/phase2/domains"
mkdir -p "${PROJECT_DIR}/results/phase2/foldseek/databases"
mkdir -p "${PROJECT_DIR}/results/phase2/foldseek/search"
mkdir -p "${PROJECT_DIR}/results/phase2/foldseek/analysis"
mkdir -p "${PROJECT_DIR}/results/phase2/foldx/per_protein"
mkdir -p "${PROJECT_DIR}/results/phase2/foldx/work"
mkdir -p "${PROJECT_DIR}/results/phase2/stability"
mkdir -p "${PROJECT_DIR}/results/phase2/stats"
mkdir -p "${PROJECT_DIR}/results/phase2/figures"
mkdir -p "${PROJECT_DIR}/tmp/foldseek_tmp"
mkdir -p "${PROJECT_DIR}/logs"
echo "Directory structure created."

# Validate Phase 1 data (rsync'd from local Mac)
echo ""
echo "Checking Phase 1 files..."
MISSING=0
for f in \
  "${PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv" \
  "${PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv" \
  "${PROJECT_DIR}/data/processed/human_mito_proteome.tsv" \
  "${PROJECT_DIR}/data/processed/human_matrix_proteome.tsv" \
  "${PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv" \
  "${PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.fasta" \
  "${PROJECT_DIR}/data/raw/uniprot/human_proteome.fasta" \
  "${PROJECT_DIR}/results/domains/cath_domain_assignments.tsv" \
  "${PROJECT_DIR}/results/domains/cath_protein_summary.tsv" \
  "${PROJECT_DIR}/results/structures/structure_index.tsv" \
  "${PROJECT_DIR}/results/structures/dssp_summary.tsv" \
  "${PROJECT_DIR}/results/termini/n_vs_c_paired.tsv" \
  "${PROJECT_DIR}/results/mts/combined_targeting.tsv"; do
  if [ -f "$f" ]; then
    echo "  [OK]  $(basename $f)"
  else
    echo "  [MISSING]  $f"
    MISSING=$((MISSING + 1))
  fi
done

if [ $MISSING -gt 0 ]; then
  echo ""
  echo "WARNING: ${MISSING} Phase 1 files missing."
  echo "Transfer from local Mac:"
  echo "  rsync -avz /Users/vishalbharti/Downloads/Antah_Asti_Prarambh/ \\"
  echo "    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/"
else
  echo ""
  echo "All Phase 1 files present."
fi

# Validate tools
echo ""
echo "Checking software..."
for cmd in python3 foldseek; do
  if command -v $cmd &>/dev/null; then
    echo "  [OK]  $cmd: $(command -v $cmd)"
  else
    echo "  [MISSING]  $cmd — install via: conda install -c bioconda foldseek"
  fi
done

# Check Chainsaw
CHAINSAW_DIR="/lustre/vishal.bharti/software/chainsaw"
if [ -f "${CHAINSAW_DIR}/get_predictions.py" ]; then
  echo "  [OK]  Chainsaw: ${CHAINSAW_DIR}"
else
  echo "  [MISSING]  Chainsaw — install via:"
  echo "    cd /lustre/vishal.bharti/software"
  echo "    git clone https://github.com/JudeWells/chainsaw.git"
fi

# Check FoldX
FOLDX_BIN="/lustre/vishal.bharti/software/foldx5/foldx"
if [ -f "${FOLDX_BIN}" ]; then
  echo "  [OK]  FoldX: ${FOLDX_BIN}"
else
  echo "  [MISSING]  FoldX — requires license from https://foldxsuite.crg.eu/"
  echo "    After obtaining, extract to /lustre/vishal.bharti/software/foldx5/"
fi

# Check Python packages
python3 -c "
import sys
pkgs = ['pandas', 'numpy', 'scipy', 'yaml', 'matplotlib', 'seaborn',
        'Bio', 'statsmodels', 'gemmi']
missing = []
for p in pkgs:
    try:
        __import__(p)
    except ImportError:
        missing.append(p)
if missing:
    print(f'  [MISSING]  Python packages: {missing}')
    print('  Install via: pip install ' + ' '.join(missing))
else:
    print('  [OK]  All Python packages available')
"

# Check disk space
echo ""
AVAIL=$(df -BG "${PROJECT_DIR}" | tail -1 | awk '{print $4}' | tr -d 'G')
echo "Available disk: ${AVAIL} GB (need ~50 GB)"
if [ "${AVAIL}" -lt 50 ]; then
  echo "WARNING: Insufficient disk space. Need at least 50 GB."
fi

echo ""
echo "Setup complete. Review output above for any [MISSING] items."
echo "============================================================"
