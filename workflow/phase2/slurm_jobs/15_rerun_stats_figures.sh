#!/bin/bash
#SBATCH --job-name=aap_rerun
#SBATCH --output=logs/15_rerun_%j.out
#SBATCH --error=logs/15_rerun_%j.err
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

# =============================================================================
# Antah Asti Prarambh — Step 15: Re-run Module H (stats) + Module I (figures)
# After steps 12-14 fill the full-proteome gaps
# Requires: Steps 12 (DSSP), 13 (Contact Order), 14 (CATH) complete
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 15: Re-run stats + figures with full-proteome data — $(date)"
echo "============================================================"

echo "This step should be submitted AFTER steps 12, 13, 14 complete."
echo "It will re-integrate all full-scale data and re-run Module H + I."
echo ""
echo "Check that these files exist:"
echo "  ${RESULTS}/structures/dssp_summary_full.tsv"
echo "  ${RESULTS}/stability/contact_order_full.tsv"
echo "  ${RESULTS}/stability/n_vs_c_rco_paired_full.tsv"
echo "  ${RESULTS}/domains/cath_domain_assignments_full.tsv"
echo ""

for f in \
  "${RESULTS}/stability/n_vs_c_rco_paired_full.tsv" \
  "${RESULTS}/foldx/foldx_stability_all.tsv"; do
  if [ ! -f "$f" ]; then
    echo "WARNING: Missing: $f"
  else
    echo "OK: $f ($(wc -l < "$f") lines)"
  fi
done

echo ""
echo "Ready for Module H + I re-run. Submit manually after verification."
echo "============================================================"
