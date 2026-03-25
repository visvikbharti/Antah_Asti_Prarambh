#!/bin/bash
#SBATCH --job-name=aap_mod_i
#SBATCH --output=logs/11_module_i_%j.out
#SBATCH --error=logs/11_module_i_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

# =============================================================================
# Antah Asti Prarambh — Step 11: Module I (Publication figures)
# 6 publication-quality figures from full-scale analysis
# Requires: Steps 9 (Module F) + 10 (Module H) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
export MPLBACKEND=Agg

echo "============================================================"
echo "  Step 11: Module I — Figures — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

python3 -u "${PROJECT_DIR}/workflow/phase2/scripts/module_i_full.py"

echo ""
echo "Module I complete — $(date)"
echo "============================================================"
