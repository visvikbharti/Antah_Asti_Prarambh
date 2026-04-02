#!/bin/bash
#SBATCH --job-name=aap_foldx_col
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/logs/foldx_collect_%j.out
#SBATCH --dependency=afterok:${1:-0}

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Run after all array jobs complete
python3 -u /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/run_foldx.py \
    --config /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/config.yaml \
    --collect

echo "Collection complete: $(date)"
