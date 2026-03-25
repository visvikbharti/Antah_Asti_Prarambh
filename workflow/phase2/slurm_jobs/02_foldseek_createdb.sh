#!/bin/bash
#SBATCH --job-name=aap_fs_createdb
#SBATCH --output=logs/02_foldseek_createdb_%j.out
#SBATCH --error=logs/02_foldseek_createdb_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

# =============================================================================
# Antah Asti Prarambh — Step 2: Foldseek create databases
# Creates structure databases from downloaded CIF files
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 2: Foldseek createdb — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

cd "${SCRIPT_DIR}"

python3 run_foldseek_full.py \
    --config "${CONFIG}" \
    --step createdb \
    --threads ${SLURM_CPUS_PER_TASK:-4}

echo ""
echo "Foldseek createdb complete — $(date)"
echo "============================================================"
