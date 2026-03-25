# Antah Asti Prarambh — Session 5 Documentation

**Date:** 2026-03-19
**Session:** 5 (Phase 2 HPC — verification, analysis chain, FoldX deployment)
**Previous sessions:** 1-2 (Phase 0+1 local pilot), 3 (Phase 2 HPC initial deployment), 4 (3 critical bug fixes)

---

## Table of Contents

1. [Session 5 Objectives](#1-session-5-objectives)
2. [Verification of Chainsaw + Module E](#2-verification-of-chainsaw--module-e)
3. [Module F → H → I Analysis Chain](#3-module-f--h--i-analysis-chain)
4. [Bugs Found and Fixed](#4-bugs-found-and-fixed)
5. [Figure Polishing](#5-figure-polishing)
6. [FoldX Installation and Deployment](#6-foldx-installation-and-deployment)
7. [Phase 2 Results Summary](#7-phase-2-results-summary)
8. [Key Scientific Findings](#8-key-scientific-findings)
9. [Complete File Inventory — Phase 2 Results](#9-complete-file-inventory--phase-2-results)
10. [HPC Job History](#10-hpc-job-history)
11. [Current Status and Next Steps](#11-current-status-and-next-steps)
12. [Prompt for Next Session](#12-prompt-for-next-session)

---

## 1. Session 5 Objectives

Session 4 left two jobs running:
- Chainsaw v2 (job 93763): 25,007 proteins with correct STRIDE binary
- Module E v2 (job 93764): domain integration with fixed column names

Session 5 goals:
1. Verify both jobs completed successfully
2. Review output files for correctness
3. Transfer Module F/H/I scripts to HPC and run analysis chain
4. Review results and transfer to Mac
5. Set up FoldX (license obtained during session)

**Actual outcome:** All goals achieved. Three minor bugs fixed in Modules H and I. Figures polished. FoldX 5.1 installed and full array job submitted (25,007 proteins, ~2.5 days).

---

## 2. Verification of Chainsaw + Module E

### Chainsaw v2 (Job 93763) — VERIFIED

| Metric | Expected | Actual | Status |
|--------|----------|--------|--------|
| Exit code | 0:0 | 0:0 | PASS |
| Proteins processed | 25,007 | 25,007 | PASS |
| Stride failures | 0 | 0 | PASS |
| Runtime | ~46 hrs | ~12.5 hrs | Faster than expected |
| Start | Wed Mar 18 13:09 IST | — | — |
| End | Thu Mar 19 01:44 IST | — | — |

Domain count distribution from Chainsaw:

| ndom | Count | Percent |
|------|-------|---------|
| 0 | 1,602 | 6.4% |
| 1 | 9,472 | 37.9% |
| 2 | 7,011 | 28.0% |
| 3 | 3,748 | 15.0% |
| 4 | 1,623 | 6.5% |
| 5+ | 1,551 | 6.2% |

### Module E v2 (Job 93764) — VERIFIED

| Metric | Expected | Actual | Status |
|--------|----------|--------|--------|
| Exit code | 0:0 | 0:0 | PASS |
| Unified assignments | ~25,000 | 25,258 | PASS |
| Accession columns | Populated | Populated | PASS |
| GroEL substrate count | 252 | 252 | PASS |
| HSP60 substrate count | 266 | 266 | PASS |
| Runtime | ~2 min | 40 sec | — |

Source breakdown: CATH = 1,390 proteins, Chainsaw = 23,868 proteins.

Verification commands used:
```bash
sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 93763,93764
grep -c "Stride failed" logs/05_chainsaw_93763.err  # → 0
head -5 unified_domain_assignments_full.tsv  # accessions populated
grep -c "groel" unified_domain_assignments_full.tsv  # → 257 (includes overlap)
grep -c "hsp60" unified_domain_assignments_full.tsv  # → 291 (includes overlap)
```

---

## 3. Module F → H → I Analysis Chain

### Module F — N-vs-C Structural Analysis (Job 93902)

| Metric | Value |
|--------|-------|
| State | COMPLETED |
| Runtime | 55 minutes |
| Proteins analyzed | 5,322 (with domain boundaries) |
| Paired comparisons | 2,648 (multi-domain, N-vs-C) |
| Contact orders computed | 11,824 |

Output files:
- `stability/region_boundaries_full.tsv` — 5,322 records
- `stability/n_vs_c_paired_full.tsv` — 2,648 records
- `stability/contact_order_full.tsv` — 11,824 records
- `stability/structure_metrics_full.tsv` — 11,824 records
- `stability/sequence_metrics_full.tsv` — 11,824 records

### Module H — Comparative Statistics (Job 93944)

| Metric | Value |
|--------|-------|
| State | COMPLETED |
| Runtime | 45 seconds |
| Total tests | 56 |
| Significant (hierarchical BH) | 25 |

Results by hypothesis family:

| Family | Significant | Total | Key findings |
|--------|-------------|-------|-------------|
| Domain architecture (H1) | 9 | 24 | GroEL enriched in TIM barrels (OR=8.4, p=1.9e-8), CATH 1.10.10.10 (OR=50.9, p=3.5e-9) |
| N-vs-C asymmetry (H2) | 14 | 30 | N>C contact order universal (p=7.1e-18 in mito_bg). NOT substrate-specific. No GroEL class effect. |
| MTS targeting (H3) | 2 | 2 | HSP60 matrix enrichment OR=3.29 (p=1.6e-16). MTS pre-domain 84.4% (p=3.4e-51). |

### Module I — Publication Figures (Job 94025)

| Metric | Value |
|--------|-------|
| State | COMPLETED |
| Runtime | 1 min 32 sec |
| Figures generated | 6/6 |

---

## 4. Bugs Found and Fixed

### Bug 1: `scipy.stats.binom_test` removed (Module H)

| Aspect | Detail |
|--------|--------|
| **Script** | `module_h_full.py`, line 452 |
| **Symptom** | `AttributeError: module 'scipy.stats' has no attribute 'binom_test'` |
| **Cause** | `binom_test` deprecated and removed in newer scipy versions |
| **Fix** | `stats.binom_test(n, total, 0.5)` → `stats.binomtest(n, total, 0.5).pvalue` |
| **Impact** | Module H crashed at Family 3 (MTS targeting) |

### Bug 2: `savefig.bbox_inches` invalid rcParam (Module I)

| Aspect | Detail |
|--------|--------|
| **Script** | `module_i_full.py`, line 38 |
| **Symptom** | `KeyError: 'savefig.bbox_inches is not a valid rc parameter'` |
| **Cause** | Older matplotlib on HPC doesn't support this rcParam |
| **Fix** | Removed from `plt.rcParams.update()`, added `bbox_inches="tight"` to each `savefig()` call |
| **Impact** | Module I crashed immediately, no figures generated |

### Bug 3: Palette key mismatch (Module I, Figure 2)

| Aspect | Detail |
|--------|--------|
| **Script** | `module_i_full.py`, line 45 |
| **Symptom** | `Error: The palette dictionary is missing keys: {'N-domain', 'C-region'}` |
| **Cause** | Palette keys were `n_domain`/`c_region` but seaborn hue values were `N-domain`/`C-region` |
| **Fix** | Changed `NC_COLORS` keys to match: `{"N-domain": "#4ECDC4", "C-region": "#FF6B6B"}` |
| **Impact** | Figure 2 not generated (5/6 figures only) |

### Preventive fix: Python output buffering

| Aspect | Detail |
|--------|--------|
| **Scripts** | All three SLURM wrappers (`09_module_f.sh`, `10_module_h.sh`, `11_module_i.sh`) |
| **Symptom** | No log output visible for 16+ minutes despite active processing |
| **Cause** | Python buffers stdout when writing to SLURM log files (non-interactive) |
| **Fix** | Changed `python3` → `python3 -u` (unbuffered output) in all SLURM scripts |

---

## 5. Figure Polishing

After initial figures were generated on HPC, a polished version was created locally using `module_i_polished.py`:

**Improvements applied (all using real data only — no fabrication):**
- Colorblind-friendly palette (Wong 2011, Nature Methods)
- Bold panel labels (A, B, C)
- Clean spines (removed top/right borders)
- Sample sizes on every panel
- Real p-values and significance annotations from Module H output
- Kruskal-Wallis p-values on GroEL class comparison
- Matrix enrichment OR on MTS targeting figure
- Binomial p-value on MTS-domain pie chart
- Clean category names (no underscores, proper capitalization)
- Pearson correlation from actual data (r=0.84, p=5.3e-13, n=45)

---

## 6. FoldX Installation and Deployment

### License and Installation

| Step | Detail |
|------|--------|
| Registration | foldxsuite.crg.eu, user: VISHAL BHARTI, CSIR-IGIB |
| Version | FoldX 5.1 (Linux 64-bit), 30.55 MB |
| Binary | `foldx_20270131` (symlinked to `foldx`) |
| Location | `/lustre/vishal.bharti/software/foldx5/` |
| License expiry | 2026-12-31 |
| Rotabase | NOT required (FoldX 5.1 feature) |

### Pipeline Setup

1. `run_foldx.py` — orchestrator script: CIF→PDB conversion (gemmi), RepairPDB, Stability calculation
2. `06_foldx_generate.sh` — generates SLURM array job script
3. `07_foldx_collect.sh` — merges per-protein JSON results into TSV
4. `config.yaml` — updated with FoldX 5.1 paths and no rotabase

### Test Run

| Metric | Value |
|--------|-------|
| Test chunk | 50 proteins (chunk 0) |
| Success rate | 50/50 (100%) |
| Runtime | 34 minutes (~40 sec/protein) |
| Total energy | First protein (A0A385XJ53): 18.81 kcal/mol |
| Component energies | null (FoldX 5.1 output format — total energy sufficient) |

### Full Run Submitted

| Metric | Value |
|--------|-------|
| Array job ID | 94152 |
| Collection job ID | 94158 (afterok:94152) |
| Total tasks | 501 |
| Proteins per task | 50 |
| Total proteins | 25,007 |
| Concurrent tasks | 5 (QOS limit) |
| Est. completion | ~2026-03-21 evening IST (~57 hours) |

---

## 7. Phase 2 Results Summary

### Data Scale

| Dataset | Phase 1 (pilot) | Phase 2 (full-scale) |
|---------|-----------------|---------------------|
| Total structures | 1,382 | 25,007 |
| E. coli structures | — | 4,371 |
| Human structures | — | ~20,636 (F1 only) |
| Foldseek clusters | 1,155 | 16,242 |
| CATH domain assignments | 1,151 | 1,390 |
| Chainsaw domain predictions | ~200 | 25,007 |
| Unified assignments | 1,387 | 25,258 |
| Paired N-vs-C comparisons | ~400 | 2,648 |
| Contact orders | ~800 | 11,824 |
| Statistical tests | 281 | 56 (streamlined) |
| Significant results | 22 | 25 |

### Output Files on Mac

| Directory | Files | Size | Description |
|-----------|-------|------|-------------|
| `results/phase2/domains/` | 4 files | 5.3 MB | Domain predictions + unified assignments |
| `results/phase2/stability/` | 5 files | 4.2 MB | Contact order, N-vs-C paired, region boundaries |
| `results/phase2/stats/` | 3 files | 981 KB | Corrected p-values, statistics summary |
| `results/phase2/figures/` | 12 files | 1.6 MB | 6 figures (PDF + PNG, 300 DPI) |
| `results/phase2/foldseek/` | 3 files | 2.3 MB | Cluster membership, summary |

---

## 8. Key Scientific Findings

### Finding 1: N-terminal domains have higher contact order than C-terminal regions — UNIVERSALLY

This is the central finding and it's robust across all protein groups:

| Dataset | n | Median N-C diff | Wilcoxon p | Effect size (r) |
|---------|---|----------------|------------|----------------|
| GroEL substrates | 124 | 0.043 | 8.9e-05 | 0.41 |
| HSP60 substrates | 131 | 0.059 | 5.3e-06 | 0.46 |
| Matrix background | 251 | 0.069 | 2.4e-09 | 0.43 |
| Mito background | 425 | 0.064 | 7.1e-18 | 0.48 |

**Critical interpretation:** This asymmetry is NOT specific to chaperonin substrates. Mann-Whitney tests comparing substrates vs. background show no significant difference (GroEL vs E. coli bg: p=0.058; HSP60 vs mito bg: p=0.536). The N>C contact order asymmetry is a general property of multi-domain proteins, consistent with co-translational folding constraints.

### Finding 2: GroEL substrates enriched in specific CATH superfamilies

| Superfamily | Description | OR | p (corrected) |
|-------------|------------|-----|---------------|
| 3.20.20.70 | TIM barrel | 8.37 | 2.3e-07 |
| 1.10.10.10 | Arc repressor-like | 50.86 | 8.3e-08 |
| 3.30.420.40 | Nucleotidyltransferase | 6.01 | 5.8e-03 |
| 3.40.640.10 | Muconolactone isomerase | 2.56 | 3.98e-02 |

GroEL CATH class distribution significantly different from background (chi2=16.79, p=2.1e-03, Cramer's V=0.089).

### Finding 3: HSP60 substrates enriched for mitochondrial matrix

| Test | Result |
|------|--------|
| HSP60 matrix enrichment | OR=3.29 [2.46-4.40], p=1.6e-16 |
| MTS pre-domain dominance | 368/436 = 84.4%, binomial p=3.4e-51 |
| Median MTS-to-domain gap | 12 residues |

### Finding 4: No GroEL substrate class effect

| Test | Statistic | p-value |
|------|-----------|---------|
| Class effect on RCO | Kruskal-Wallis H=0.52 | 0.77 |
| Class effect on pLDDT | Kruskal-Wallis H=0.17 | 0.92 |

All three GroEL classes (I, II, III) show similar N>C contact order asymmetry.

### Finding 5: N-domain contact order conserved across species

Pearson correlation of N-domain RCO between GroEL↔HSP60 homolog pairs:
- r = 0.84, p = 5.3e-13, n = 45 pairs

This is strong evidence that the N-terminal folding complexity is evolutionarily conserved.

---

## 9. Complete File Inventory — Phase 2 Results

### On Mac (`~/Downloads/Antah_Asti_Prarambh/results/phase2/`)

```
results/phase2/
├── domains/
│   ├── chainsaw_full_predictions.tsv           (25,007 rows)
│   ├── chainsaw_full_predictions_annotated.tsv  (25,007 rows)
│   ├── unified_domain_assignments_full.tsv      (25,258 rows)
│   └── domain_distribution_full.tsv             (56 rows)
├── stability/
│   ├── region_boundaries_full.tsv               (5,322 rows)
│   ├── n_vs_c_paired_full.tsv                   (2,648 rows)
│   ├── contact_order_full.tsv                   (11,824 rows)
│   ├── structure_metrics_full.tsv               (11,824 rows)
│   └── sequence_metrics_full.tsv                (11,824 rows)
├── stats/
│   ├── corrected_pvalues_full.tsv               (56 rows)
│   ├── statistics_summary_full.txt              (534 lines)
│   └── stability_comparisons_full.tsv           (2,648 rows)
├── figures/
│   ├── fig1_domain_architecture.{pdf,png}
│   ├── fig2_n_vs_c_stability.{pdf,png}
│   ├── fig3_groel_class_comparison.{pdf,png}
│   ├── fig4_mts_targeting.{pdf,png}
│   ├── fig5_orthology.{pdf,png}
│   └── fig6_summary.{pdf,png}
├── foldseek/
│   └── analysis/
│       ├── combined_cluster_membership.tsv      (27,063 rows)
│       ├── foldseek_clusters_full.tsv           (16,242 clusters)
│       └── foldseek_full_summary.txt
└── foldx/                                       (PENDING — job 94152 running)
    ├── per_protein/{accession}.json             (50 so far, 25,007 expected)
    └── foldx_stability_all.tsv                  (after collection)
```

### On HPC (`/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/`)

Same as above, plus:
- Raw AlphaFold structures: `data/raw/alphafold/full/ecoli/` (4,371 CIF) and `data/raw/alphafold/full/human/` (23,672 CIF)
- FoldX working directories: `results/phase2/foldx/work/{accession}/`
- SLURM logs: `workflow/phase2/slurm_jobs/logs/`

---

## 10. HPC Job History

### Session 5 Jobs

| Job ID | Name | State | ExitCode | Elapsed | Purpose |
|--------|------|-------|----------|---------|---------|
| 93763 | aap_chainsaw | COMPLETED | 0:0 | 12:34:36 | Chainsaw v2 (25,007 proteins) |
| 93764 | aap_mod_e | COMPLETED | 0:0 | 00:00:40 | Module E v2 (domain integration) |
| 93902 | aap_mod_f | COMPLETED | 0:0 | 00:55:00 | Module F (N-vs-C analysis) |
| 93903 | aap_mod_h | FAILED | 1:0 | 00:01:07 | Module H v1 (binom_test bug) |
| 93905 | aap_mod_i | — | — | — | Module I (DependencyNeverSatisfied, cancelled) |
| 93944 | aap_mod_h | COMPLETED | 0:0 | 00:00:45 | Module H v2 (fixed) |
| 93945 | aap_mod_i | FAILED | 1:0 | 00:00:12 | Module I v1 (bbox_inches bug) |
| 93973 | aap_mod_i | COMPLETED | 0:0 | 00:00:19 | Module I v2 (5/6 figures, palette bug) |
| 94025 | aap_mod_i | COMPLETED | 0:0 | 00:01:32 | Module I v3 (6/6 figures) |
| 94138 | aap_foldx_gen | COMPLETED | 0:0 | 00:00:52 | FoldX generate array script |
| 94141_0 | aap_foldx | COMPLETED | 0:0 | 00:33:57 | FoldX test (50 proteins, 50/50 success) |
| 94152 | aap_foldx | RUNNING | — | — | FoldX full (501 tasks, 25,007 proteins) |
| 94158 | aap_foldx_col | PENDING | — | — | FoldX collect (afterok:94152) |

### Cumulative Session History

| Session | Date | Key Achievement |
|---------|------|----------------|
| 1-2 | 2026-02-xx | Phase 0+1 local pilot (1,390 proteins) |
| 3 | 2026-03-17 | Phase 2 HPC deployment (AlphaFold, Foldseek, initial Chainsaw) |
| 4 | 2026-03-18 | Found/fixed 3 critical bugs (STRIDE, decompression, column names) |
| 5 | 2026-03-19 | Verified all results, ran analysis chain, installed FoldX, submitted full run |

---

## 11. Current Status and Next Steps

### Pipeline Status

```
Phase 2 HPC Pipeline Progress (as of 2026-03-19 17:30 IST)
=============================================================

[DONE] Step 01: AlphaFold download           — 25,007 structures
[DONE] Step 02: Foldseek createdb            — 6 min
[DONE] Step 03: Foldseek search              — 16 min
[DONE] Step 04: Foldseek cluster             — 3 min → 16,242 clusters
[DONE] Step 05a: Decompress human CIFs       — 24 min → 23,672 files
[DONE] Step 05: Chainsaw domain prediction   — 12.5 hrs → 25,007 proteins, 0 Stride failures
[DONE] Step 08: Module E (domain integration)— 40 sec → 25,258 unified assignments
[DONE] Step 09: Module F (N-vs-C analysis)   — 55 min → 2,648 paired comparisons
[DONE] Step 10: Module H (statistics)        — 45 sec → 56 tests, 25 significant
[DONE] Step 11: Module I (figures)           — 1.5 min → 6/6 publication figures
[RUN]  Step 06: FoldX stability              — ~57 hrs → job 94152 (501 array tasks)
[PEND] Step 07: FoldX collect                — job 94158 (afterok:94152)

Legend: [DONE] Complete  [RUN] Running  [PEND] Pending
```

### After FoldX Completes (~2026-03-21)

1. Verify FoldX completion (success rate, check chunk summaries)
2. Re-run Module F with FoldX DeltaG integration
3. Re-run Module H with DeltaG-enriched statistical tests
4. Re-run Module I to update figures with DeltaG data
5. Transfer updated results to Mac
6. Prepare final collaborator delivery package

---

## 12. Prompt for Next Session

> Continuing Antah Asti Prarambh Phase 2 (session 6). Check memory files for context.
>
> **Session 5 summary:** Verified Chainsaw (25,007 proteins, 0 Stride failures) + Module E (25,258 assignments). Ran F→H→I chain (3 bug fixes: binomtest, bbox_inches, palette). Polished 6 figures. Installed FoldX 5.1 on HPC. Submitted full FoldX array job (94152, 501 tasks).
>
> **FoldX status:**
> ```bash
> sacct --format=JobID,JobName,State,ExitCode,Elapsed -j 94152,94158
> ls /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/chunks/ | wc -l
> ls /lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/results/phase2/foldx/per_protein/ | wc -l
> ```
>
> [PASTE output here]

---
