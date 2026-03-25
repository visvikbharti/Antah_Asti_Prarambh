#!/bin/bash
# =============================================================================
# Antah Asti Prarambh — Phase 2: Master SLURM Pipeline Submission
# Submits all jobs with proper dependency chains
#
# Usage:
#   cd /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs
#   bash submit_pipeline.sh
#
# QOS: common (max 300G mem, 5-day walltime, max 5 concurrent jobs)
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SLURM_DIR}"

# Create logs directory (on Lustre)
mkdir -p logs

echo "============================================================"
echo "  Antah Asti Prarambh — Phase 2 Pipeline Submission"
echo "  The End is the Beginning"
echo "  Submitting from: ${SLURM_DIR}"
echo "============================================================"
echo ""

# Step 0: Setup and validation
echo "Submitting Step 0: Setup..."
JOB_SETUP=$(sbatch --parsable 00_setup.sh)
echo "  Job ${JOB_SETUP}: setup"

# Step 1: Download AlphaFold (depends on setup)
echo "Submitting Step 1: Download AlphaFold..."
JOB_DOWNLOAD=$(sbatch --parsable --dependency=afterok:${JOB_SETUP} 01_download_alphafold.sh)
echo "  Job ${JOB_DOWNLOAD}: download_alphafold (after ${JOB_SETUP})"

# Steps 2-5 run in parallel after download completes
# (max 3 concurrent: Foldseek createdb + Chainsaw + FoldX generate → within QOS limit of 5)

# --- Foldseek branch ---
echo "Submitting Step 2: Foldseek createdb..."
JOB_FS_DB=$(sbatch --parsable --dependency=afterok:${JOB_DOWNLOAD} 02_foldseek_createdb.sh)
echo "  Job ${JOB_FS_DB}: foldseek_createdb (after ${JOB_DOWNLOAD})"

echo "Submitting Step 3: Foldseek search..."
JOB_FS_SEARCH=$(sbatch --parsable --dependency=afterok:${JOB_FS_DB} 03_foldseek_search.sh)
echo "  Job ${JOB_FS_SEARCH}: foldseek_search (after ${JOB_FS_DB})"

echo "Submitting Step 4: Foldseek cluster..."
JOB_FS_CLUSTER=$(sbatch --parsable --dependency=afterok:${JOB_FS_SEARCH} 04_foldseek_cluster.sh)
echo "  Job ${JOB_FS_CLUSTER}: foldseek_cluster (after ${JOB_FS_SEARCH})"

# --- Chainsaw branch (parallel with Foldseek) ---
echo "Submitting Step 5: Chainsaw..."
JOB_CHAINSAW=$(sbatch --parsable --dependency=afterok:${JOB_DOWNLOAD} 05_chainsaw.sh)
echo "  Job ${JOB_CHAINSAW}: chainsaw (after ${JOB_DOWNLOAD})"

# --- FoldX branch (parallel with Foldseek and Chainsaw) ---
echo "Submitting Step 6a: FoldX generate array script..."
JOB_FX_GEN=$(sbatch --parsable --dependency=afterok:${JOB_DOWNLOAD} 06_foldx_generate.sh)
echo "  Job ${JOB_FX_GEN}: foldx_generate (after ${JOB_DOWNLOAD})"

# NOTE: FoldX array job must be submitted manually after 06 completes!
echo ""
echo "  *** MANUAL STEP REQUIRED ***"
echo "  After job ${JOB_FX_GEN} completes, manually submit the FoldX array job:"
echo "    cat /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/foldx_array.slurm"
echo "    ARRAY_JOB=\$(sbatch --parsable /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/foldx_array.slurm)"
echo "    COLLECT_JOB=\$(sbatch --parsable --dependency=afterok:\${ARRAY_JOB} 07_foldx_collect.sh)"
echo "    bash submit_analysis.sh \${COLLECT_JOB}"
echo ""

# --- Module E: needs Foldseek clusters + Chainsaw ---
echo "Submitting Step 8: Module E (domains)..."
JOB_MOD_E=$(sbatch --parsable --dependency=afterok:${JOB_FS_CLUSTER}:${JOB_CHAINSAW} 08_module_e_domains.sh)
echo "  Job ${JOB_MOD_E}: module_e (after foldseek_cluster + chainsaw)"

echo ""
echo "============================================================"
echo "  Jobs submitted. Pipeline dependency graph:"
echo ""
echo "  ${JOB_SETUP} (setup)"
echo "    └─> ${JOB_DOWNLOAD} (download)"
echo "          ├─> ${JOB_FS_DB} -> ${JOB_FS_SEARCH} -> ${JOB_FS_CLUSTER} (foldseek)"
echo "          ├─> ${JOB_CHAINSAW} (chainsaw)"
echo "          └─> ${JOB_FX_GEN} (foldx generate) -> [MANUAL: array + collect]"
echo ""
echo "          foldseek_cluster + chainsaw -> ${JOB_MOD_E} (module E)"
echo ""
echo "  Remaining analysis (Module F, H, I) submitted via:"
echo "    bash submit_analysis.sh <FOLDX_COLLECT_JOB_ID>"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Logs:    ls ${SLURM_DIR}/logs/"
echo "============================================================"
