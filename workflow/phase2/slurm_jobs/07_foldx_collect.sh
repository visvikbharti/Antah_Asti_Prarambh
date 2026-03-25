#!/bin/bash
#SBATCH --job-name=aap_foldx_collect
#SBATCH --output=logs/07_foldx_collect_%j.out
#SBATCH --error=logs/07_foldx_collect_%j.err
#SBATCH --partition=compute
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

# =============================================================================
# Antah Asti Prarambh — Step 6b: Collect FoldX array job results
# Run AFTER all FoldX array tasks complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
SCRIPT_DIR="${PROJECT_DIR}/workflow/phase2"
CONFIG="${SCRIPT_DIR}/config.yaml"

echo "============================================================"
echo "  Step 6b: Collect FoldX results — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

cd "${SCRIPT_DIR}"

python3 run_foldx.py \
    --config "${CONFIG}" \
    --collect

FOLDX_OUT="${PROJECT_DIR}/results/phase2/foldx/foldx_stability_all.tsv"
if [ -f "${FOLDX_OUT}" ]; then
    NPROTEINS=$(tail -n +2 "${FOLDX_OUT}" | wc -l)
    echo "Collected results for ${NPROTEINS} proteins"
else
    echo "WARNING: Output file not found at ${FOLDX_OUT}"
fi

echo ""
echo "FoldX collection complete — $(date)"
echo "============================================================"
