#!/bin/bash
#SBATCH --job-name=aap_fs_cluster
#SBATCH --output=logs/04_foldseek_cluster_%j.out
#SBATCH --error=logs/04_foldseek_cluster_%j.err
#SBATCH --partition=compute
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00

# =============================================================================
# Antah Asti Prarambh — Step 4: Foldseek cluster + export + analyze
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 4: Foldseek cluster — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

cd "${SCRIPT_DIR}"

# Cluster
python3 run_foldseek_full.py \
    --config "${CONFIG}" \
    --step cluster \
    --threads ${SLURM_CPUS_PER_TASK:-16}

# Export results
python3 run_foldseek_full.py \
    --config "${CONFIG}" \
    --step export

# Analyze
python3 run_foldseek_full.py \
    --config "${CONFIG}" \
    --step analyze-only

echo ""
echo "Foldseek cluster + analysis complete — $(date)"
echo "============================================================"
