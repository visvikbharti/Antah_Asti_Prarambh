#!/bin/bash
#SBATCH --job-name=aap_foldx
#SBATCH --partition=compute
#SBATCH --array=0-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/logs/foldx_%A_%a.out
#SBATCH --error=/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/logs/foldx_%A_%a.err

# ============================================================================
# FoldX Stability Calculations — SLURM Array Job
# Antah Asti Prarambh Phase 2
#
# Total proteins: 25007
# Chunk size: 50
# Total chunks: 501
# ============================================================================

echo "=== FoldX Array Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Run chunk
python3 -u /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/run_foldx.py \
    --config /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/config.yaml \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --chunk-size 50

echo ""
echo "End: $(date)"
echo "Exit code: $?"
