#!/bin/bash
# =============================================================================
# Antah Asti Prarambh — Phase 2: Submit analysis jobs (Modules F, H, I)
# Run after FoldX collection and Module E are complete.
#
# Usage:
#   bash submit_analysis.sh <FOLDX_COLLECT_JOB_ID>
#
# If FoldX collect and Module E are already done, use:
#   bash submit_analysis.sh done
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SLURM_DIR}"

if [ $# -lt 1 ]; then
    echo "Usage: bash submit_analysis.sh <FOLDX_COLLECT_JOB_ID>"
    echo "       bash submit_analysis.sh done  (if upstream jobs already complete)"
    exit 1
fi

FOLDX_JOB="$1"

echo "============================================================"
echo "  Submitting analysis pipeline (Modules F, H, I)"
echo "============================================================"

if [ "${FOLDX_JOB}" = "done" ]; then
    # No dependency — submit immediately
    echo "Submitting Step 9: Module F (stability)..."
    JOB_MOD_F=$(sbatch --parsable 09_module_f_stability.sh)

    echo "Submitting Step 10: Module H (statistics)..."
    JOB_MOD_H=$(sbatch --parsable --dependency=afterok:${JOB_MOD_F} 10_module_h_stats.sh)

    echo "Submitting Step 11: Module I (figures)..."
    JOB_MOD_I=$(sbatch --parsable --dependency=afterok:${JOB_MOD_H} 11_module_i_figures.sh)
else
    # Module F depends on FoldX collect
    echo "Submitting Step 9: Module F (stability, after FoldX job ${FOLDX_JOB})..."
    JOB_MOD_F=$(sbatch --parsable --dependency=afterok:${FOLDX_JOB} 09_module_f_stability.sh)

    echo "Submitting Step 10: Module H (statistics, after Module F)..."
    JOB_MOD_H=$(sbatch --parsable --dependency=afterok:${JOB_MOD_F} 10_module_h_stats.sh)

    echo "Submitting Step 11: Module I (figures, after Module H)..."
    JOB_MOD_I=$(sbatch --parsable --dependency=afterok:${JOB_MOD_H} 11_module_i_figures.sh)
fi

echo ""
echo "  ${JOB_MOD_F} (Module F: stability)"
echo "    └─> ${JOB_MOD_H} (Module H: statistics)"
echo "          └─> ${JOB_MOD_I} (Module I: figures)"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "============================================================"
