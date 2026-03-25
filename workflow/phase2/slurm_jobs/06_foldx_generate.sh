#!/bin/bash
#SBATCH --job-name=aap_foldx_gen
#SBATCH --output=logs/06_foldx_generate_%j.out
#SBATCH --error=logs/06_foldx_generate_%j.err
#SBATCH --partition=compute
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

# =============================================================================
# Antah Asti Prarambh — Step 6a: Generate FoldX array job script
# This creates the SLURM array job file, then you submit it separately
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 6a: Generate FoldX SLURM array script — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Check FoldX exists
if [ ! -f "/lustre/vishal.bharti/software/foldx5/foldx" ]; then
    echo "ERROR: FoldX not found at /lustre/vishal.bharti/software/foldx5/foldx"
    echo "FoldX requires a license from https://foldxsuite.crg.eu/"
    echo "After obtaining, extract to /lustre/vishal.bharti/software/foldx5/"
    exit 1
fi

cd "${SCRIPT_DIR}"

python3 run_foldx.py \
    --config "${CONFIG}" \
    --generate-slurm

FOLDX_DIR="${PROJECT_DIR}/results/phase2/foldx"
echo ""
echo "Generated SLURM array script at: ${FOLDX_DIR}/foldx_array.slurm"
echo ""
echo "Next steps:"
echo "  1. Review: cat ${FOLDX_DIR}/foldx_array.slurm"
echo "  2. Submit: sbatch ${FOLDX_DIR}/foldx_array.slurm"
echo "  3. After completion, submit: sbatch workflow/phase2/slurm_jobs/07_foldx_collect.sh"
echo "============================================================"
