#!/bin/bash
#SBATCH --job-name=aap_chainsaw
#SBATCH --output=logs/05_chainsaw_%j.out
#SBATCH --error=logs/05_chainsaw_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=3-00:00:00

# =============================================================================
# Antah Asti Prarambh — Step 5: Chainsaw domain prediction (v2 — STRIDE fixed)
# Runs on all ~25,000 structures at ~7s each = ~49 hours (3-day walltime)
# v2 fix: correct STRIDE binary (heiniglab/stride, not bioconda StriDe)
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
CHAINSAW_DIR="/lustre/vishal.bharti/software/chainsaw"

# Directories for structures
AF_ECOLI="${PROJECT_DIR}/data/raw/alphafold/full/ecoli"
AF_HUMAN="${PROJECT_DIR}/data/raw/alphafold/full/human"
OUTPUT="${PROJECT_DIR}/results/phase2/domains/chainsaw_full_predictions.tsv"

echo "============================================================"
echo "  Step 5: Chainsaw domain prediction — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Fix: prioritize conda's libstdc++ over system gcc-8.4 (which lacks GLIBCXX_3.4.29)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Validate Chainsaw installation
if [ ! -f "${CHAINSAW_DIR}/get_predictions.py" ]; then
    echo "ERROR: Chainsaw not found at ${CHAINSAW_DIR}"
    echo "Install: cd /lustre/vishal.bharti/software && git clone https://github.com/JudeWells/chainsaw.git"
    exit 1
fi

mkdir -p "$(dirname ${OUTPUT})"

# Clean old batch files from previous (broken) run
rm -f "${PROJECT_DIR}/results/phase2/domains/chainsaw_batch_"*.tsv

# Collect all CIF files
TMPLIST=$(mktemp -p "${PROJECT_DIR}/tmp")
find "${AF_ECOLI}" "${AF_HUMAN}" -name "AF-*-F1-model_v*.cif" 2>/dev/null > "${TMPLIST}"
TOTAL=$(wc -l < "${TMPLIST}")
echo "Total structures to process: ${TOTAL}"

# Process in batches of 500
BATCH_SIZE=500
BATCH_ID=0
BATCH_DIR="${PROJECT_DIR}/results/phase2/domains"
BATCH_FILES=()

while IFS= read -r line; do
    BATCH_FILES+=("$line")
    if [ ${#BATCH_FILES[@]} -ge ${BATCH_SIZE} ]; then
        echo "Processing batch ${BATCH_ID} (${#BATCH_FILES[@]} files)... $(date)"

        # Create temp directory with symlinks (on Lustre, not /tmp)
        BATCH_TMPDIR=$(mktemp -d -p "${PROJECT_DIR}/tmp" "chainsaw_batch_${BATCH_ID}_XXXX")
        for cif in "${BATCH_FILES[@]}"; do
            ln -sf "$cif" "${BATCH_TMPDIR}/$(basename $cif)"
        done

        BATCH_OUT="${BATCH_DIR}/chainsaw_batch_$(printf '%04d' ${BATCH_ID}).tsv"
        python3 "${CHAINSAW_DIR}/get_predictions.py" \
            --structure_directory "${BATCH_TMPDIR}" \
            --output "${BATCH_OUT}" \
            --use_first_chain || echo "Batch ${BATCH_ID} had errors"

        rm -rf "${BATCH_TMPDIR}"
        BATCH_FILES=()
        BATCH_ID=$((BATCH_ID + 1))
    fi
done < "${TMPLIST}"

# Process remaining files
if [ ${#BATCH_FILES[@]} -gt 0 ]; then
    echo "Processing final batch ${BATCH_ID} (${#BATCH_FILES[@]} files)... $(date)"
    BATCH_TMPDIR=$(mktemp -d -p "${PROJECT_DIR}/tmp" "chainsaw_batch_${BATCH_ID}_XXXX")
    for cif in "${BATCH_FILES[@]}"; do
        ln -sf "$cif" "${BATCH_TMPDIR}/$(basename $cif)"
    done

    BATCH_OUT="${BATCH_DIR}/chainsaw_batch_$(printf '%04d' ${BATCH_ID}).tsv"
    python3 "${CHAINSAW_DIR}/get_predictions.py" \
        --structure_directory "${BATCH_TMPDIR}" \
        --output "${BATCH_OUT}" \
        --use_first_chain || echo "Final batch had errors"

    rm -rf "${BATCH_TMPDIR}"
fi

rm -f "${TMPLIST}"

# Merge all batch results
echo "Merging batch results..."
HEADER_DONE=0
for batch_file in "${BATCH_DIR}"/chainsaw_batch_*.tsv; do
    if [ -f "$batch_file" ]; then
        if [ ${HEADER_DONE} -eq 0 ]; then
            cat "$batch_file" > "${OUTPUT}"
            HEADER_DONE=1
        else
            tail -n +2 "$batch_file" >> "${OUTPUT}"
        fi
    fi
done

NPROTEINS=$(tail -n +2 "${OUTPUT}" | wc -l)
echo "Chainsaw complete: ${NPROTEINS} proteins in ${OUTPUT}"
echo "Finished: $(date)"
echo "============================================================"
