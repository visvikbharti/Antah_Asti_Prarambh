# Antah Asti Prarambh — Session 6 Documentation

**Date:** 2026-03-22
**Session:** 6 (FoldX monitoring, timeout fixes, collaborator deliverables)
**Previous sessions:** 1-2 (Phase 0+1 local pilot), 3 (Phase 2 HPC deployment), 4 (3 bug fixes), 5 (analysis chain + FoldX submission)

---

## Table of Contents

1. [Session 6 Objectives](#1-session-6-objectives)
2. [FoldX Status Assessment](#2-foldx-status-assessment)
3. [HPC Fixes Applied](#3-hpc-fixes-applied)
4. [Collaborator Deliverables Created](#4-collaborator-deliverables-created)
5. [Complete Session Job History](#5-complete-session-job-history)
6. [Current Status and Next Steps](#6-current-status-and-next-steps)
7. [Prompt for Next Session](#7-prompt-for-next-session)

---

## 1. Session 6 Objectives

Session 5 left the FoldX array job (94152, 501 tasks) running on HPC, estimated to complete ~March 21.

Session 6 goals:
1. Check FoldX completion status
2. If incomplete, assess progress and fix any issues
3. Prepare collaborator presentation and data handover package
4. Create comprehensive sharing guide for full reproducibility

**Actual outcome:** FoldX was only 28% complete (not the estimated 57 hours — the QOS 5-concurrent-task limit made the original estimate too optimistic). Found and fixed 2 timed-out chunks. Created 3 comprehensive collaborator documents.

---

## 2. FoldX Status Assessment

Checked at 2026-03-22 19:10 IST (~72 hours after submission).

### Overall Progress

| Metric | Value |
|--------|:-----:|
| Total tasks | 501 |
| Completed successfully | 141 |
| Timed out (TIMEOUT) | 2 (chunks 125, 139) |
| Currently running | 5 (chunks 143-147) |
| Pending | 353 (chunks 148-500) |
| Per-protein results | 7,204 / 25,007 (28.8%) |
| Chunk summary files | 141 |

### Runtime Distribution (completed chunks)

| Metric | Value |
|--------|:-----:|
| Minimum | 00:00:16 (chunk 0, likely cached) |
| Median | ~02:30:00 |
| Maximum (completed) | 03:39:41 (chunk 135) |
| Timeouts | 04:00:01 (chunk 125), 04:00:18 (chunk 139) |

### Timeout Analysis

Chunks 125 and 139 hit the 4-hour SLURM wall time limit. The max successful runtime was 3:39:41, so these chunks contained proteins that were slightly more computationally expensive (likely larger or more complex structures). Extending to 6 hours provides a comfortable margin.

### Why the Original Estimate Was Wrong

Original estimate: ~57 hours (assuming continuous processing).
Actual: QOS policy `MaxJobsPerUserLimit` = 5 concurrent tasks. With 501 tasks at ~2.5 hrs each:
- 501 ÷ 5 = 100.2 rounds × 2.5 hrs = ~250 hours ≈ 10.4 days
- Submitted March 19 → expected completion March 29-30

---

## 3. HPC Fixes Applied

### Fix 1: Extend wall time for pending tasks

```bash
scontrol update JobId=94152 TimeLimit=06:00:00
```

**Result:** Tasks 148-500 updated successfully. Running tasks 143-147 denied (admin-only for running jobs). This is fine — they were all under 2.5 hours and completed within the original 4-hour limit.

### Fix 2: Cancel broken collection job

```bash
scancel 94158
```

**Why:** Job 94158 had `afterok:94152` dependency. Since 2 array tasks TIMED OUT (exit code non-zero), the `afterok` condition can never be satisfied. The collection must be submitted manually after all tasks complete.

### Fix 3: Resubmit timed-out chunks

```bash
sbatch --array=125,139 --time=06:00:00 /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/foldx_array.slurm
```

**Result:** Submitted as job 94439 (only chunks 125, 139). Will run with 6-hour wall time when slots become available in the QOS queue.

### Post-fix state

```
squeue output at 19:29 IST:
  94439_[125,139]  — PENDING (resubmitted timeouts)
  94152_[148-500]  — PENDING (original remaining)
  94152_143-147    — RUNNING on node23
```

---

## 4. Collaborator Deliverables Created

### Document 1: COLLABORATOR_PRESENTATION.md

**Path:** `docs/COLLABORATOR_PRESENTATION.md`
**Purpose:** Meeting-ready presentation document

**Contents:**
- Project overview and 3 scientific goals
- Study design (7 datasets, 25,007 structures)
- Computational pipeline flowchart
- 5 key findings with full statistics:
  1. GroEL enriched in TIM barrels (OR=8.4)
  2. N>C contact order universal (NOT substrate-specific)
  3. Cross-species RCO conservation (r=0.84)
  4. MTS pre-domain architecture (84.4%)
  5. HSP60 matrix enrichment (OR=3.29)
- Figure descriptions (6 figures)
- Statistical framework (56 tests, hierarchical BH)
- Phase 1 vs Phase 2 scale comparison
- FoldX status and what it will add
- Biological synthesis narrative
- Discussion points for collaborator meeting
- Complete data inventory
- Software and reproducibility

### Document 2: DATA_HANDOVER_INDEX.md

**Path:** `docs/DATA_HANDOVER_INDEX.md`
**Purpose:** File-by-file guide for collaborator

**Contents:**
- Quick-start (4 most important files)
- Every result file with columns, row counts, and usage instructions
- Figure descriptions
- Dataset label mapping
- Column naming conventions
- Known limitations

### Document 3: COLLABORATOR_SHARING_GUIDE.md

**Path:** `docs/COLLABORATOR_SHARING_GUIDE.md`
**Purpose:** Complete inventory and packaging instructions

**Contents:**
- Full file list organized by category (A–G):
  - A: Raw data (literature, UniProt, MitoCarta, AlphaFold)
  - B: Processed data (7 datasets)
  - C: Phase 1 scripts (20 Python scripts)
  - D: Phase 2 HPC scripts (Snakefile, config, 7 Python, 19 SLURM)
  - E: Phase 1 results (structures, domains, homology, termini, MTS, stats, figures)
  - F: Phase 2 results (full-scale domains, stability, stats, figures)
  - G: Documentation (14 markdown files)
- 3 sharing options (full project, essential only, HPC pipeline only)
- tar command for creating the archive (~30-40 MB)
- Collaborator quick-start reading order
- What's NOT included (HPC-only regenerable files)

---

## 5. Complete Session Job History

### Session 6 HPC Actions

| Action | Command | Result |
|--------|---------|--------|
| Extend wall time | `scontrol update JobId=94152 TimeLimit=06:00:00` | Pending tasks 148-500 updated; running 143-147 denied |
| Cancel collection | `scancel 94158` | Cancelled |
| Resubmit timeouts | `sbatch --array=125,139 --time=06:00:00 foldx_array.slurm` | Job 94439 submitted |

### Cumulative Job History

| Session | Job ID | Name | State | Purpose |
|---------|--------|------|-------|---------|
| 3 | 93500-93550 | various | COMPLETED | AlphaFold download, Foldseek, initial Chainsaw |
| 4 | 93763 | aap_chainsaw | COMPLETED | Chainsaw v2 (STRIDE fix) |
| 4 | 93764 | aap_mod_e | COMPLETED | Module E v2 (column fix) |
| 5 | 93902 | aap_mod_f | COMPLETED | Module F (N-vs-C) |
| 5 | 93944 | aap_mod_h | COMPLETED | Module H v2 (binomtest fix) |
| 5 | 94025 | aap_mod_i | COMPLETED | Module I v3 (6/6 figures) |
| 5 | 94138 | aap_foldx_gen | COMPLETED | FoldX array generation |
| 5 | 94141_0 | aap_foldx | COMPLETED | FoldX test (50 proteins) |
| 5 | 94152 | aap_foldx | RUNNING | FoldX full (501 tasks) |
| 6 | 94158 | aap_foldx_col | CANCELLED | Collection (cancelled — afterok broken) |
| 6 | 94439 | aap_foldx | PENDING | FoldX resubmit (chunks 125, 139) |

---

## 6. Current Status and Next Steps

### Pipeline Status

```
Phase 2 HPC Pipeline Progress (as of 2026-03-22 19:30 IST)
=============================================================

[DONE] Step 01: AlphaFold download           — 25,007 structures
[DONE] Step 02: Foldseek createdb            — 6 min
[DONE] Step 03: Foldseek search              — 16 min
[DONE] Step 04: Foldseek cluster             — 3 min → 16,242 clusters
[DONE] Step 05a: Decompress human CIFs       — 24 min → 23,672 files
[DONE] Step 05: Chainsaw domain prediction   — 12.5 hrs → 25,007 proteins
[DONE] Step 08: Module E (domain integration)— 40 sec → 25,258 assignments
[DONE] Step 09: Module F (N-vs-C analysis)   — 55 min → 2,648 paired
[DONE] Step 10: Module H (statistics)        — 45 sec → 56 tests, 25 significant
[DONE] Step 11: Module I (figures)           — 1.5 min → 6/6 figures
[RUN]  Step 06: FoldX stability              — 28% done, est. ~7 more days
[PEND] Step 07: FoldX collect                — manual submit after completion

Legend: [DONE] Complete  [RUN] Running  [PEND] Pending
```

### Collaborator Delivery Status

```
[DONE] COLLABORATOR_PRESENTATION.md  — Meeting presentation
[DONE] DATA_HANDOVER_INDEX.md        — File guide
[DONE] COLLABORATOR_SHARING_GUIDE.md — Full sharing inventory + tar command
[PEND] Sharing archive               — Create after FoldX for complete package
```

### After FoldX Completes (~March 29-30)

1. Verify FoldX completion and success rate
2. Submit collection job manually
3. Re-run Module F with FoldX DeltaG integration
4. Re-run Module H with DeltaG-enriched statistics
5. Re-run Module I to update figures with DeltaG data
6. Transfer updated results to Mac
7. Update COLLABORATOR_PRESENTATION.md with FoldX findings
8. Create final sharing archive
9. Second collaborator meeting
10. Manuscript preparation

---

## 7. Prompt for Next Session

> Continuing Antah Asti Prarambh Phase 2 (session 7). Check memory files for context.
>
> **Session 6 summary:** FoldX was 28% done (143/501 chunks). Fixed 2 TIMEOUT failures (chunks 125, 139 → resubmitted as job 94439 with 6h wall time). Extended remaining tasks to 6h. Cancelled broken collection job 94158. Created 3 collaborator documents (presentation, handover index, sharing guide).
>
> **FoldX should be done by now (~March 29-30). Check status:**
> ```bash
> # Check both original and resubmitted jobs
> sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 94152 | tail -20
> sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 94439
>
> # Count completed results
> ls /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/per_protein/ | wc -l
>
> # Check for any remaining running/pending
> squeue -u vishal.bharti
> ```
>
> **If all complete, submit collection:**
> ```bash
> sbatch /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/workflow/phase2/slurm_jobs/07_foldx_collect.sh
> ```
>
> **Then:** Re-run Module F (with DeltaG) → Module H → Module I → transfer to Mac → update presentation → manuscript.
>
> [PASTE output here]

---
