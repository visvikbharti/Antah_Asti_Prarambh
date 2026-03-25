#!/bin/bash
#SBATCH --job-name=aap_mod_f
#SBATCH --output=logs/09_module_f_%j.out
#SBATCH --error=logs/09_module_f_%j.err
#SBATCH --partition=compute
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00

# =============================================================================
# Antah Asti Prarambh — Step 9: Module F (N-vs-C structural analysis)
# Computes contact order, pLDDT, secondary structure per region
# Uses Stride for SS, gemmi for CIF parsing, Plaxco formula for CO
# NO FoldX dependency — contact order is the primary folding kinetics proxy
# Requires: Step 8 (Module E domains) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"

echo "============================================================"
echo "  Step 9: Module F — N-vs-C structural analysis — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Check inputs
for f in \
  "${PROJECT_DIR}/results/phase2/domains/unified_domain_assignments_full.tsv" \
  "${PROJECT_DIR}/results/phase2/domains/chainsaw_full_predictions.tsv"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing input: $f"
    exit 1
  fi
done

# Ensure tmp directory exists
mkdir -p "${PROJECT_DIR}/tmp"

python3 -u "${PROJECT_DIR}/workflow/phase2/scripts/module_f_full.py"

echo ""
echo "Module F complete — $(date)"
echo "============================================================"
