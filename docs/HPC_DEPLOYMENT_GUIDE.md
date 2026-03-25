# Antah Asti Prarambh -- Phase 2 HPC Deployment Guide

**Purpose:** Step-by-step instructions to deploy and run the full-scale proteomics
pipeline (4,403 E. coli + 20,416 Human proteins) on HPC.

**HPC:** IGIB HPC (tejas.igib.res.in), CentOS 7, SLURM, Lustre filesystem
**Date updated:** 2026-03-17
**Estimated wall-clock time:** 2-4 days
**Estimated disk usage:** ~50 GB on Lustre

---

## Quick Reference

| Setting | Value |
|---------|-------|
| Username | `vishal.bharti` |
| Login node | `hn1` via `ssh vishal.bharti@tejas.igib.res.in` |
| Project dir | `/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc` |
| Partition | `compute` (40 CPUs, 380G RAM/node) |
| QOS | `common` (max 300G mem, 5-day walltime, max 5 concurrent jobs) |
| Conda | `/home/vishal.bharti/miniconda3` |
| Conda env | `proteomics` (created by 00_setup.sh) |
| Chainsaw | `/lustre/vishal.bharti/software/chainsaw` |
| FoldX | `/lustre/vishal.bharti/software/foldx5` (needs license) |

---

## Step 1: Transfer Project to HPC

From your local Mac:

```bash
rsync -avz --progress \
    /Users/vishalbharti/Downloads/Antah_Asti_Prarambh/ \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/
```

This transfers all Phase 1 data, results, scripts, and Phase 2 pipeline scripts.
Estimated transfer: ~700 MB (excluding AlphaFold CIF files from pilot, which
are not needed on HPC since Phase 2 downloads fresh structures).

To exclude the pilot CIF files (saves ~466 MB):

```bash
rsync -avz --progress \
    --exclude='data/raw/alphafold/pilot/' \
    /Users/vishalbharti/Downloads/Antah_Asti_Prarambh/ \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/
```

---

## Step 2: Install Software on HPC

SSH into the HPC:

```bash
ssh vishal.bharti@tejas.igib.res.in
```

### 2a. Install Chainsaw (domain prediction)

```bash
cd /lustre/vishal.bharti/software
git clone https://github.com/JudeWells/chainsaw.git
cd chainsaw
# Chainsaw needs torch + einops (torch is in base env, einops will be in proteomics env)
```

### 2b. Install FoldX (thermodynamic stability)

FoldX requires a free academic license:
1. Register at https://foldxsuite.crg.eu/
2. Download FoldX 5 for Linux
3. Extract to `/lustre/vishal.bharti/software/foldx5/`

```bash
mkdir -p /lustre/vishal.bharti/software/foldx5
cd /lustre/vishal.bharti/software/foldx5
# Copy the downloaded foldx archive here and extract:
# tar -xzf foldx5Linux.tar.gz
# Verify:
ls -la foldx rotabase.txt
```

**Note:** If you don't have FoldX yet, Steps 0-5 and 8 (Module E) can still run.
FoldX is only needed for Steps 6-7 and Module F. You can deploy and start the
pipeline today while waiting for the FoldX license.

---

## Step 3: Verify Transfer and Submit Pipeline

```bash
cd /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs

# Verify scripts are there
ls -la *.sh

# Submit the full pipeline
bash submit_pipeline.sh
```

This submits all jobs with proper dependency chains:

```
Job 12345 (setup)
  └─> Job 12346 (download AlphaFold ~22 GB)
        ├─> Job 12347 -> 12348 -> 12349 (Foldseek: createdb -> search -> cluster)
        ├─> Job 12350 (Chainsaw domain prediction)
        └─> Job 12351 (FoldX generate array) -> [MANUAL: array + collect]

        foldseek_cluster + chainsaw -> Job 12352 (Module E: domains)
```

---

## Step 4: Monitor Progress

```bash
# See your running/pending jobs
squeue -u $USER

# See completed jobs with resource usage
sacct --format=JobID,JobName,State,Elapsed,MaxRSS -j <JOBID>

# Check specific job output
cat /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/logs/01_download_*.out

# Watch for job completion
watch -n 60 squeue -u $USER
```

---

## Step 5: FoldX (Manual Step)

After job `06_foldx_generate` completes (check with `squeue`):

```bash
cd /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs

# Review the generated array script
cat /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/foldx_array.slurm

# Submit the array job
ARRAY_JOB=$(sbatch --parsable /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/foldx_array.slurm)
echo "FoldX array job: $ARRAY_JOB"

# Submit the collect job (runs after all array tasks finish)
COLLECT_JOB=$(sbatch --parsable --dependency=afterok:${ARRAY_JOB} 07_foldx_collect.sh)
echo "Collect job: $COLLECT_JOB"
```

---

## Step 6: Submit Analysis Jobs

After FoldX collect completes:

```bash
cd /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs

# Submit Module F -> Module H -> Module I chain
bash submit_analysis.sh ${COLLECT_JOB}

# OR if FoldX collect is already done:
bash submit_analysis.sh done
```

This submits:
- Module F: N-vs-C stability analysis with FoldX DeltaG
- Module H: Hierarchical statistics
- Module I: Publication figures

---

## Step 7: Check Final Results

After all jobs complete:

```bash
ls /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/

# Expected:
#   structures/structure_index_full.tsv
#   foldseek/analysis/foldseek_clusters_full.tsv
#   domains/chainsaw_full_predictions.tsv
#   domains/unified_domain_assignments_full.tsv
#   domains/domain_distribution_full.tsv
#   foldx/foldx_stability_all.tsv
#   stability/n_vs_c_paired_full.tsv
#   stats/corrected_pvalues_full.tsv
#   stats/statistics_summary_full.txt
#   figures/*.pdf, *.png
```

### Transfer results back to Mac

```bash
# On your local Mac:
rsync -avz --progress \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/ \
    /Users/vishalbharti/Downloads/Antah_Asti_Prarambh/results/phase2/
```

---

## Troubleshooting

### Job failed — check logs
```bash
cat /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/logs/<jobname>_<jobid>.err
```

### Out of memory
Increase `--mem` in the specific SLURM script. Max is 300G on `compute` partition.

### Out of time (5-day QOS limit)
The `common` QOS allows up to 5 days. If Foldseek search exceeds this, reduce
sensitivity in `config.yaml` (from 7.5 to 5.5) or use `newcompute` partition.

### Conda environment issues
```bash
source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
# If env doesn't exist, re-run setup:
sbatch 00_setup.sh
```

### Foldseek not found
```bash
conda activate proteomics
conda install -c bioconda foldseek -y
which foldseek
```

### Chainsaw errors (CUDA/GPU)
Chainsaw uses PyTorch but runs on CPU by default. If you see CUDA errors:
```bash
export CUDA_VISIBLE_DEVICES=""
```

### Dependency job cancelled (due to upstream failure)
```bash
# Check what failed
sacct --format=JobID,JobName,State -j <JOBID>
# Fix the issue, then resubmit the failed job and downstream:
sbatch <failed_script>.sh
```

---

## Pipeline Summary

| When | What to do |
|------|------------|
| Once | Step 1: rsync project to HPC |
| Once | Step 2: Install Chainsaw (+ FoldX when license arrives) |
| Once | Step 3: `bash submit_pipeline.sh` |
| Wait | Steps 0-5 + 8 run automatically via SLURM dependencies |
| Manual | Step 5: Submit FoldX array + collect job |
| Manual | Step 6: `bash submit_analysis.sh <COLLECT_JOB_ID>` |
| Done | Step 7: Check results, transfer back to Mac |

**Estimated compute:** download (1-4 hrs) + Foldseek (8-24 hrs) + Chainsaw (~24-48 hrs)
+ FoldX (~4 hrs as array) + analysis (~2 hrs) = **~2-4 days total**

---

*All paths are hardcoded to your IGIB HPC environment. No placeholders to replace.*
