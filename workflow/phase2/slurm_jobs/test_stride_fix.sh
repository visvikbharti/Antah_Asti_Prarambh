#!/bin/bash
#SBATCH --job-name=test_stride
#SBATCH --output=logs/test_stride_%j.out
#SBATCH --error=logs/test_stride_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00

set -euo pipefail
source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CHAINSAW_DIR="/lustre/vishal.bharti/software/chainsaw"
TESTDIR="${PROJECT_DIR}/tmp/test_single"

rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"

# Test with 3 known multi-domain proteins
for acc in P0ABB0 P00934 P04693; do
    ln -sf ${PROJECT_DIR}/data/raw/alphafold/full/ecoli/AF-${acc}-F1-model_v6.cif ${TESTDIR}/
done

echo "Testing Chainsaw with fixed STRIDE..."
python3 ${CHAINSAW_DIR}/get_predictions.py \
    --structure_directory ${TESTDIR} \
    --output ${PROJECT_DIR}/tmp/test_stride_out.tsv \
    --use_first_chain 2>&1

echo ""
echo "=== Results ==="
cat ${PROJECT_DIR}/tmp/test_stride_out.tsv
