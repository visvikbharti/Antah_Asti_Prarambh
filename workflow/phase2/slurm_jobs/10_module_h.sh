#!/bin/bash
#SBATCH --job-name=aap_mod_h
#SBATCH --output=logs/10_module_h_%j.out
#SBATCH --error=logs/10_module_h_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

# =============================================================================
# Antah Asti Prarambh — Step 10: Module H (Comparative statistics)
# Hierarchical BH correction across 3 hypothesis families
# Requires: Steps 8 (Module E) + 9 (Module F) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"

echo "============================================================"
echo "  Step 10: Module H — Statistics — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

for f in \
  "${PROJECT_DIR}/results/phase2/domains/unified_domain_assignments_full.tsv" \
  "${PROJECT_DIR}/results/phase2/stability/n_vs_c_paired_full.tsv"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing input: $f"
    exit 1
  fi
done

python3 -u "${PROJECT_DIR}/workflow/phase2/scripts/module_h_full.py"

echo ""
echo "Module H complete — $(date)"
echo "============================================================"
