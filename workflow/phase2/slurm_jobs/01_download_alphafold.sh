#!/bin/bash
#SBATCH --job-name=aap_download
#SBATCH --output=logs/01_download_%j.out
#SBATCH --error=logs/01_download_%j.err
#SBATCH --partition=compute
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

# =============================================================================
# Antah Asti Prarambh — Step 1: Download full AlphaFold proteomes
# ~2 GB E. coli + ~20 GB Human = ~22 GB total
# Resume-safe: re-run skips completed downloads
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 1: Download AlphaFold — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

cd "${SCRIPT_DIR}"

python3 download_alphafold_full.py \
    --config "${CONFIG}" \
    --organism both \
    --workers 8

echo ""
echo "Download complete — $(date)"
echo "============================================================"
