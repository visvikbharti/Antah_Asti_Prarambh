# Antah Asti Prarambh — Session 4 Documentation

**Date:** 2026-03-18
**Session:** 4 (continuation of Phase 2 HPC deployment)
**Previous sessions:** 1-2 (Phase 0+1 local pilot), 3 (Phase 2 HPC initial deployment)

---

## Table of Contents

1. [Session 4 Objectives](#1-session-4-objectives)
2. [What We Found — Three Critical Bugs](#2-what-we-found--three-critical-bugs)
3. [Debugging Process — Step by Step](#3-debugging-process--step-by-step)
4. [Fixes Applied](#4-fixes-applied)
5. [New Scripts Written](#5-new-scripts-written)
6. [Current HPC Job Status](#6-current-hpc-job-status)
7. [Complete Directory Structure — HPC](#7-complete-directory-structure--hpc)
8. [Complete Directory Structure — Local Mac](#8-complete-directory-structure--local-mac)
9. [File Transfer Commands](#9-file-transfer-commands)
10. [Column Name Reference (Critical)](#10-column-name-reference-critical)
11. [HPC Environment Reference](#11-hpc-environment-reference)
12. [TODO List for Session 5](#12-todo-list-for-session-5)
13. [Overall Pipeline Status](#13-overall-pipeline-status)
14. [Prompt for Next Session](#14-prompt-for-next-session)

---

## 1. Session 4 Objectives

Session 3 left two SLURM jobs running on HPC:
- **Chainsaw (job 93689):** Domain prediction on 4,457 structures (~8 hrs)
- **Module E (job 93700):** Domain integration (dependency on Chainsaw)

Session 4 goals:
1. Check if both jobs completed successfully
2. Review output files for correctness
3. Help with FoldX setup (license pending)
4. Submit remaining analysis chain (Modules F→H→I)
5. Transfer results back to Mac

**Actual outcome:** Both jobs completed but outputs were fundamentally broken. Three critical bugs were discovered, fixed, and the pipeline was re-submitted.

---

## 2. What We Found — Three Critical Bugs

### Bug 1: Wrong STRIDE Binary (CRITICAL — Root Cause)

| Aspect | Detail |
|--------|--------|
| **Symptom** | ALL 4,457 Chainsaw predictions had `ndom=0`, `chopping=NULL` |
| **Diagnosis** | `grep -c "Stride failed"` in error log returned **4,457** — 100% failure rate |
| **Root cause** | `conda install -c bioconda stride` installed "StriDe" v0.0.1 (genomics string tool by Yao-Ting Huang), NOT the protein secondary structure STRIDE (Frishman & Argos, 1995) |
| **Evidence** | Running `stride/stride 2>&1 | head -3` showed "Program: StriDe / Version: 0.0.1 / Contact: ythuang@cs.ccu.edu.tw" |
| **Fix** | Compiled correct STRIDE from `github.com/heiniglab/stride`: `git clone`, `make` |
| **Verification** | Test job 93750: P0ABB0 got 3 domains (28-96, 97-375, 384-510), P04693 got 2 domains, P00934 got 2 domains. Stride failures: **0**. |
| **Impact** | 8 hours of compute wasted (session 3 Chainsaw run was entirely useless) |

### Bug 2: Human AlphaFold CIFs Still Compressed

| Aspect | Detail |
|--------|--------|
| **Symptom** | Only 4,457 structures found (should be ~25,007) |
| **Diagnosis** | `ls human/*.cif.gz | head -5` showed all human files were `.cif.gz` (gzipped) |
| **Root cause** | AlphaFold bulk download produces `.cif.gz` + `.pdb.gz` files. E. coli was extracted (4,371 `.cif`), human was not. The `find *.cif` command in Chainsaw script missed all compressed files. |
| **Fix** | SLURM job `05a_decompress_human.sh`: `gunzip *.cif.gz`, `rm *.pdb.gz` |
| **Result** | 23,672 human CIF files extracted (includes F1+F2+F3 fragments; pipeline filters to F1 via `AF-*-F1-model_v*.cif` pattern). Total: 25,007 F1 structures. |

### Bug 3: Module E Column Name Mismatches (6 Wrong Lookups)

| Aspect | Detail |
|--------|--------|
| **Symptom** | Module E log showed: "Proteins with CATH: 0", "Substrate counts: GroEL=0, HSP60=0". All accessions empty. |
| **Diagnosis** | Script used `"accession"` for all files, but actual column names differ |
| **Root cause** | Phase 1 files use inconsistent naming conventions |

Column name mapping (what the script expected vs reality):

| File | Script Used | Actual Column |
|------|------------|---------------|
| `cath_domain_assignments.tsv` | `accession` | `uniprot_accession` |
| `cath_protein_summary.tsv` | `accession` | `uniprot_accession` |
| `groel_substrates_standardized.tsv` | `accession` | `current_accession` |
| `hsp60_tier1_substrates.tsv` | `accession` | `uniprot_id` |
| `human_matrix_proteome.tsv` | `accession` | `uniprot_id` |
| `human_mito_proteome.tsv` | `accession` | `uniprot_id` |

**Fix:** Rewrote Module E with correct column names and fallback detection patterns.

---

## 3. Debugging Process — Step by Step

### Phase 1: Verify job completion
```bash
sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 93689,93700
# Both COMPLETED with ExitCode 0:0
# Chainsaw: 07:44:44, Module E: 00:01:07
```

### Phase 2: Check outputs
```bash
# Chainsaw: all predictions are ndom=0
awk -F'\t' 'NR>1{print $4}' chainsaw_full_predictions.tsv | sort | uniq -c
# Output: 4457 0   ← ALL ZERO

# Known multi-domain proteins also ndom=0:
grep "P0ABB0" chainsaw_full_predictions.tsv
# P0ABB0 (AtpA, 513 res, 3 known domains): ndom=0, chopping=NULL ← WRONG

# Module E: empty accessions
head -5 unified_domain_assignments_full.tsv
# accession column is blank for all rows
```

### Phase 3: Identify root causes
```bash
# Stride failure count = total proteins = 100% failure
grep -c "Stride failed" logs/05_chainsaw_93689.err
# Output: 4457

# Check Stride binary identity
stride/stride 2>&1 | head -3
# "Program: StriDe / Version: 0.0.1" ← WRONG PROGRAM

# Check human directory
ls data/raw/alphafold/full/human/ | head -5
# AF-A0A024R1R8-F1-model_v6.cif.gz ← Still compressed!

# Module E log
cat logs/08_module_e_93700.out
# "Proteins with CATH: 0" ← Column name mismatch
```

### Phase 4: Fix and verify
```bash
# Install correct STRIDE
git clone https://github.com/heiniglab/stride.git stride_real
cd stride_real && make
./stride 2>&1 | head -3
# "Action: secondary structure assignment" ← CORRECT

# Test Chainsaw with fixed STRIDE (SLURM test job 93750)
# Result: P0ABB0=3 domains, P04693=2 domains, P00934=2 domains ← CORRECT

# Submit fix chain
JOB1=$(sbatch --parsable 05a_decompress_human.sh)  # 93762
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} 05_chainsaw.sh)  # 93763
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 08_module_e_domains.sh)  # 93764
```

### Phase 5: Verify fix is working
```bash
# Decompression completed: 23,672 human CIF files
tail -5 logs/05a_decompress_93762.out
# "After: 23672 .cif files"

# Chainsaw found all structures:
head -3 logs/05_chainsaw_93763.out
# "Total structures to process: 25007" ← CORRECT (was 4457 before)

# Zero Stride failures:
grep -c "Stride failed" logs/05_chainsaw_93763.err
# 0 ← PERFECT (was 4457 before)
```

---

## 4. Fixes Applied

### Scripts modified on HPC:

| Script | Change | Version |
|--------|--------|---------|
| `05_chainsaw.sh` | Walltime 2d→3d, cleans old batch files | v2 |
| `08_module_e_domains.sh` | All 6 column names fixed, fallback detection | v2 |

### New scripts created:

| Script | Purpose |
|--------|---------|
| `05a_decompress_human.sh` | Gunzip human CIFs, remove PDB.gz |
| `test_stride_fix.sh` | Verify STRIDE fix on 3 known multi-domain proteins |

### Software installed on HPC:

| Software | Path | Source |
|----------|------|--------|
| STRIDE (correct) | `/lustre/vishal.bharti/software/stride_real/stride` | `github.com/heiniglab/stride` |
| Symlink | `/lustre/vishal.bharti/software/chainsaw/stride/stride` → above | |

---

## 5. New Scripts Written

### Complete rewrites of Modules F, H, I (locally, not yet on HPC)

The old versions (09_module_f_stability.sh, 10_module_h_stats.sh, 11_module_i_figures.sh) were:
- Entirely FoldX-dependent (FoldX not available)
- Missing core computations (no CIF reading, no contact order, no Stride SS)
- Had the same column name bugs as Module E
- Only 4 basic figures instead of 6

**New Module F** (`workflow/phase2/scripts/module_f_full.py` + `slurm_jobs/09_module_f.sh`):
- Reads CIF files with gemmi for coordinates + pLDDT
- Runs Stride for per-residue secondary structure
- Computes contact order (Plaxco et al. 1998 formula: CA-CA < 8Å, seq_sep ≥ 6)
- Parses BOTH CATH and Chainsaw domain boundaries (correct chopping format)
- Defines three regions: pre-domain tail, N-domain, C-region
- Computes per-region: sequence composition, SS fractions, pLDDT stats, contact order
- Creates paired N-vs-C comparison table
- Processes substrates + all E. coli + all human mito (~6,000 proteins, ~4-8 hours)
- Walltime: 1 day, 32 GB RAM

**New Module H** (`workflow/phase2/scripts/module_h_full.py` + `slurm_jobs/10_module_h.sh`):
- Family 1: Domain architecture (Fisher's exact per superfamily, chi-squared for CATH class)
- Family 2: N-vs-C asymmetry (Wilcoxon signed-rank paired, Mann-Whitney U between groups, Kruskal-Wallis for GroEL class)
- Family 3: MTS targeting (Fisher's exact, binomial for pre-domain dominance)
- Hierarchical BH correction (within-family + between-family)
- Effect sizes: rank-biserial r, odds ratios with 95% CI, Cramér's V, eta-squared
- Correct column names throughout
- Walltime: 4 hours, 16 GB RAM

**New Module I** (`workflow/phase2/scripts/module_i_full.py` + `slurm_jobs/11_module_i.sh`):
- Fig 1: Domain architecture (CATH class stacked bar, top superfamilies, domain count)
- Fig 2: N-vs-C stability (contact order violins, pLDDT violins, N-C difference heatmap)
- Fig 3: GroEL class effects (CO by class, pLDDT by class with box+strip plots)
- Fig 4: MTS targeting (classification bars, gap histogram, MTS vs domain scatter)
- Fig 5: Orthology (orthogroup overlap, N-domain RCO conservation scatter)
- Fig 6: Summary (mean N-C RCO bars, significance overview, MTS pie)
- All colorblind-friendly, PDF+PNG at 300 DPI
- Walltime: 2 hours, 8 GB RAM

---

## 6. Current HPC Job Status

As of 2026-03-18 ~14:00 IST:

| Job ID | Name | Status | Started | Est. Duration | Est. Completion |
|--------|------|--------|---------|---------------|-----------------|
| 93762 | aap_gunzip (decompress) | **COMPLETED** | ~12:45 IST | 24 min | Done |
| 93763 | aap_chainsaw (Chainsaw v2) | **RUNNING** | ~13:09 IST | ~46 hours | ~Mar 20 11:00 AM |
| 93764 | aap_mod_e (Module E v2) | PENDING (dep) | After 93763 | ~2 min | ~Mar 20 11:02 AM |

### Monitoring commands:
```bash
# Check job status
squeue -u vishal.bharti
sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 93763,93764

# Check Chainsaw progress (batch numbers advance every ~55 min)
tail -5 /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/logs/05_chainsaw_93763.out

# Verify zero Stride failures
grep -c "Stride failed" /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/logs/05_chainsaw_93763.err

# Check total structures (should be 25007)
head -3 /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/logs/05_chainsaw_93763.out
```

---

## 7. Complete Directory Structure — HPC

```
/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/
├── data/
│   ├── raw/
│   │   ├── alphafold/
│   │   │   ├── full/
│   │   │   │   ├── ecoli/          # 4,371 .cif files (uncompressed)
│   │   │   │   └── human/          # 23,672 .cif files (freshly decompressed)
│   │   │   └── pilot/              # ~1,382 .cif files (Phase 1)
│   │   ├── uniprot/
│   │   │   ├── ecoli_k12_proteome.fasta + .tsv
│   │   │   └── human_proteome.fasta + .tsv
│   │   └── mitocarta/
│   │       └── Human.MitoCarta3.0.xls
│   └── processed/
│       ├── groel_substrates_standardized.tsv    # col: current_accession
│       ├── hsp60_tier1_substrates.tsv           # col: uniprot_id
│       ├── human_mito_proteome.tsv              # col: uniprot_id
│       ├── human_matrix_proteome.tsv            # col: uniprot_id
│       └── groel_hsp60_homologs.tsv
├── results/
│   ├── domains/                   # Phase 1 pilot CATH results
│   │   ├── cath_domain_assignments.tsv    # col: uniprot_accession
│   │   └── cath_protein_summary.tsv       # col: uniprot_accession
│   ├── homology/                  # Phase 1 RBH + orthogroups
│   ├── structures/                # Phase 1 DSSP + structure index
│   ├── termini/                   # Phase 1 N-vs-C results
│   ├── mts/                       # Phase 1 MTS targeting
│   ├── stats/                     # Phase 1 statistics
│   ├── figures/                   # Phase 1 figures
│   └── phase2/                    # ★ Phase 2 full-scale results
│       ├── domains/
│       │   ├── chainsaw_full_predictions.tsv          # ← being generated by job 93763
│       │   ├── chainsaw_batch_XXXX.tsv                # ← intermediate batch files
│       │   ├── unified_domain_assignments_full.tsv    # ← will be generated by job 93764
│       │   └── domain_distribution_full.tsv           # ← will be generated by job 93764
│       ├── foldseek/
│       │   ├── db/                                    # Foldseek database files
│       │   ├── search/                                # Foldseek search results
│       │   └── analysis/
│       │       └── foldseek_clusters_full.tsv         # ✓ Complete (16,193 clusters)
│       ├── stability/             # ← Module F output (not yet run)
│       ├── stats/                 # ← Module H output (not yet run)
│       └── figures/               # ← Module I output (not yet run)
├── workflow/
│   └── phase2/
│       ├── slurm_jobs/
│       │   ├── logs/              # All SLURM output/error logs
│       │   ├── 05_chainsaw.sh     # v2: 3-day walltime, cleans old batches
│       │   ├── 05a_decompress_human.sh  # NEW: gunzip human CIFs
│       │   ├── 08_module_e_domains.sh   # v2: column names fixed
│       │   ├── 09_module_f.sh     # NEW: calls module_f_full.py
│       │   ├── 10_module_h.sh     # NEW: calls module_h_full.py
│       │   ├── 11_module_i.sh     # NEW: calls module_i_full.py
│       │   └── test_stride_fix.sh # Test script for STRIDE verification
│       └── scripts/               # ← needs to be created on HPC + files transferred
│           ├── module_f_full.py   # ← TO TRANSFER
│           ├── module_h_full.py   # ← TO TRANSFER
│           └── module_i_full.py   # ← TO TRANSFER
├── tmp/                           # Temporary files (Chainsaw batches, Stride)
└── docs/                          # Documentation
```

### Software locations on HPC:
```
/lustre/vishal.bharti/software/
├── chainsaw/                      # Chainsaw domain predictor
│   ├── get_predictions.py
│   ├── saved_models/
│   └── stride/
│       └── stride → ../../stride_real/stride   # SYMLINK to correct STRIDE
├── stride_real/                   # Correct STRIDE (heiniglab/stride, compiled)
│   └── stride                     # Binary
├── foldseek/
│   └── bin/foldseek              # Standalone binary (v10-941cd33)
└── foldx5/                       # FoldX placeholder (needs license)
```

---

## 8. Complete Directory Structure — Local Mac

```
~/Downloads/Antah_Asti_Prarambh/
├── workflow/phase2/
│   ├── slurm_jobs/
│   │   ├── 05_chainsaw.sh          # v2 (updated)
│   │   ├── 05a_decompress_human.sh  # NEW
│   │   ├── 08_module_e_domains.sh   # v2 (updated)
│   │   ├── 09_module_f.sh          # NEW (wrapper)
│   │   ├── 09_module_f_stability.sh # OLD (FoldX-dependent, obsolete)
│   │   ├── 10_module_h.sh          # NEW (wrapper)
│   │   ├── 10_module_h_stats.sh     # OLD (obsolete)
│   │   ├── 11_module_i.sh          # NEW (wrapper)
│   │   ├── 11_module_i_figures.sh   # OLD (obsolete)
│   │   └── test_stride_fix.sh
│   └── scripts/
│       ├── module_f_full.py         # ★ NEW — full N-vs-C structural analysis
│       ├── module_h_full.py         # ★ NEW — comparative statistics
│       └── module_i_full.py         # ★ NEW — 6 publication figures
├── docs/
│   ├── SESSION_CONTINUITY.md        # Sessions 1-2 handoff
│   ├── SESSION4_DOCUMENTATION.md    # ★ THIS FILE
│   ├── HPC_DEPLOYMENT_GUIDE.md
│   ├── PHASE1_VERIFICATION.md
│   └── ...
└── (all Phase 1 data and results)
```

---

## 9. File Transfer Commands

### Transfer Module F/H/I scripts to HPC (required before submitting):
```bash
# From Mac terminal:
scp ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/scripts/module_f_full.py \
    ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/scripts/module_h_full.py \
    ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/scripts/module_i_full.py \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/scripts/

scp ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/slurm_jobs/09_module_f.sh \
    ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/slurm_jobs/10_module_h.sh \
    ~/Downloads/Antah_Asti_Prarambh/workflow/phase2/slurm_jobs/11_module_i.sh \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/
```

### Submit analysis chain (after Module E completes):
```bash
cd /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs
mkdir -p /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/scripts

JOB_F=$(sbatch --parsable 09_module_f.sh)
JOB_H=$(sbatch --parsable --dependency=afterok:${JOB_F} 10_module_h.sh)
JOB_I=$(sbatch --parsable --dependency=afterok:${JOB_H} 11_module_i.sh)
echo "F=$JOB_F  H=$JOB_H  I=$JOB_I"
```

### Transfer results back to Mac (after everything completes):
```bash
# From Mac terminal:
rsync -avz --progress \
    vishal.bharti@tejas.igib.res.in:/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/ \
    ~/Downloads/Antah_Asti_Prarambh/results/phase2/
```

---

## 10. Column Name Reference (Critical)

This table is essential for any future script that integrates multiple data files:

| File | Accession Column | Other Key Columns |
|------|-----------------|-------------------|
| `cath_domain_assignments.tsv` | `uniprot_accession` | `domain_start`, `domain_end`, `cath_superfamily`, `cath_class` |
| `cath_protein_summary.tsv` | `uniprot_accession` | `n_domains`, `protein_length`, `has_cath_assignment` |
| `groel_substrates_standardized.tsv` | `current_accession` | `groel_class`, `length`, `location_category` |
| `hsp60_tier1_substrates.tsv` | `uniprot_id` | `gene_name`, `median_silac_ratio`, `percent_occupancy` |
| `human_mito_proteome.tsv` | `uniprot_id` | `sub_mito_localization`, `is_matrix` |
| `human_matrix_proteome.tsv` | `uniprot_id` | `sub_mito_localization`, `is_matrix` |
| `ecoli_k12_proteome.tsv` | `Entry` | `Length`, `Subcellular location [CC]` |
| `human_proteome.tsv` | `Entry` | `Length` |
| `groel_hsp60_homologs.tsv` | (paired) | `groel_accession`, `hsp60_accession` |
| `structure_index.tsv` | `uniprot_accession` | `source_dataset`, `residues_modeled` |
| `combined_targeting.tsv` | `uniprot_accession` | `targeting_classification`, `mitocarta_is_matrix` |
| `mts_domain_relationship.tsv` | `uniprot_accession` | `transit_peptide_end`, `first_domain_start`, `mts_is_pre_domain` |
| `chainsaw_full_predictions.tsv` | `chain_id` (extract w/ regex) | `ndom`, `chopping`, `confidence` |

---

## 11. HPC Environment Reference

| Setting | Value |
|---------|-------|
| **Hostname** | `tejas.igib.res.in` |
| **Username** | `vishal.bharti` |
| **Project dir** | `/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc` |
| **Partition** | `compute` (40 CPUs, 380G RAM/node) |
| **QOS** | `common` (max 300G, 5-day walltime, max 5 concurrent) |
| **Conda** | `/home/vishal.bharti/miniconda3` (init with `source .../conda.sh`) |
| **Conda env** | `proteomics` (Python 3.11, gemmi 0.7.5, torch, pandas, scipy, statsmodels) |
| **LD_LIBRARY_PATH fix** | `export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"` |
| **STRIDE** | `/lustre/vishal.bharti/software/stride_real/stride` |
| **Chainsaw** | `/lustre/vishal.bharti/software/chainsaw` |
| **Foldseek** | `/lustre/vishal.bharti/software/foldseek/bin/foldseek` |
| **FoldX** | NOT INSTALLED (needs license from foldxsuite.crg.eu) |
| **libmamba-solver** | BROKEN (use `--solver=classic` for conda installs) |

---

## 12. TODO List for Session 5

### Immediate (when Chainsaw+Module E finish):
- [ ] Verify Chainsaw job 93763 completed (exit code 0)
- [ ] Check Stride failure count (should be 0 or near-zero)
- [ ] Check total proteins processed (should be 25,007)
- [ ] Verify Module E job 93764 completed
- [ ] Review Module E outputs: unified_domain_assignments_full.tsv, domain_distribution_full.tsv
- [ ] Check that accession columns are populated (not empty!)
- [ ] Check substrate counts are correct (GroEL=252, HSP60=266)
- [ ] Transfer Module F/H/I scripts to HPC (see Section 9)
- [ ] Create scripts directory on HPC: `mkdir -p workflow/phase2/scripts`
- [ ] Submit analysis chain: F → H → I (see Section 9)

### After Module F/H/I complete (~5-9 hours after submission):
- [ ] Review Module F outputs: region_boundaries, contact_order, n_vs_c_paired
- [ ] Review Module H outputs: corrected_pvalues, statistics_summary
- [ ] Review Module I outputs: 6 figures (PDF+PNG)
- [ ] Transfer all Phase 2 results to Mac (see Section 9)

### FoldX (if license obtained):
- [ ] Download FoldX binary to `/lustre/vishal.bharti/software/foldx5/`
- [ ] Place license file (`rotabase.txt`) in foldx5 directory
- [ ] Run FoldX jobs (06_foldx_generate.sh + 07_foldx_collect.sh)
- [ ] Re-run Module F with FoldX DeltaG integration

### Longer-term:
- [ ] Write comprehensive DOCUMENTATION.md v3.0 with Phase 2 results
- [ ] Compare Phase 1 pilot vs Phase 2 full-scale results
- [ ] Prepare manuscript figures from Phase 2 outputs
- [ ] Consider targetp/signalp for expanded MTS prediction (need DTU license)

---

## 13. Overall Pipeline Status

```
Phase 2 HPC Pipeline Progress
==============================

[✓] Step 01: AlphaFold download       — 25,007 structures (job 92733)
[✓] Step 02: Foldseek createdb        — 6 min (job 92734)
[✓] Step 03: Foldseek search          — 16 min (job 92735)
[✓] Step 04: Foldseek cluster         — 3 min (job 92737) → 16,193 clusters
[✓] Step 05a: Decompress human CIFs   — 24 min (job 93762) → 23,672 files
[⏳] Step 05: Chainsaw domain prediction — ~46 hrs (job 93763) → RUNNING
[⏳] Step 08: Module E (domain integration) — ~2 min (job 93764) → PENDING
[—] Step 06: FoldX generate           — SKIPPED (no license)
[—] Step 07: FoldX collect            — SKIPPED (no license)
[○] Step 09: Module F (N-vs-C analysis) — ~4-8 hrs — scripts ready, not submitted
[○] Step 10: Module H (statistics)     — ~10 min — scripts ready, not submitted
[○] Step 11: Module I (figures)        — ~5 min — scripts ready, not submitted

Legend: [✓] Complete  [⏳] Running/Pending  [—] Skipped  [○] Ready to submit
```

---

## 14. Prompt for Next Session

Use this prompt when starting session 5:

---

> Continuing Antah Asti Prarambh Phase 2 HPC deployment (session 5). Please check your memory files for full context, especially `project_session4_progress.md` and the updated `MEMORY.md`.
>
> **Session 4 summary:** Found and fixed 3 critical bugs (wrong STRIDE binary, compressed human CIFs, Module E column mismatches). Re-submitted Chainsaw (job 93763, 25,007 structures) + Module E (job 93764) as dependency chain. Also wrote complete rewrites of Modules F, H, I (locally, not yet on HPC).
>
> **Current status:** I'll paste the `sacct` output for jobs 93763 and 93764 below. Module F/H/I scripts are on my Mac at `workflow/phase2/scripts/` and `workflow/phase2/slurm_jobs/` — not yet transferred to HPC.
>
> Please help me:
> 1. Verify Chainsaw and Module E completed successfully (I'll paste logs)
> 2. Review the output files for correctness
> 3. Transfer Module F/H/I scripts to HPC and submit the analysis chain
> 4. Once F/H/I complete, review results and transfer to Mac
> 5. FoldX license status: [HAVE IT / STILL PENDING]
>
> Key HPC paths: project at `/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc`, SLURM scripts at `workflow/phase2/slurm_jobs/`, logs in `logs/` subdirectory.
>
> [PASTE sacct output here]
> [PASTE tail of Chainsaw log here]
> [PASTE tail of Module E log here]

---
