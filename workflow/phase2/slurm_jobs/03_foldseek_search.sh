#!/bin/bash
#SBATCH --job-name=aap_fs_search
#SBATCH --output=logs/03_foldseek_search_%j.out
#SBATCH --error=logs/03_foldseek_search_%j.err
#SBATCH --partition=compute
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00

# =============================================================================
# Antah Asti Prarambh — Step 3: Foldseek all-vs-all search
# Most resource-intensive step: 64 GB RAM, 16 CPUs, up to 24 hours
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 3: Foldseek search — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

cd "${SCRIPT_DIR}"

python3 run_foldseek_full.py \
    --config "${CONFIG}" \
    --step search \
    --threads ${SLURM_CPUS_PER_TASK:-16} \
    --memory 64G

echo ""
echo "Foldseek search complete — $(date)"
echo "============================================================"
