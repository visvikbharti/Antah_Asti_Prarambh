# Antah Asti Prarambh — Comprehensive Project Document

## "The End is the Beginning": A Comparative Structural Proteomics Study of Chaperonin Substrates

**Principal Investigator:** Vishal Bharti, CSIR-Institute of Genomics and Integrative Biology (IGIB), New Delhi
**Date:** April 1, 2026 (updated)
**Status:** Phase 2 COMPLETE including FoldX thermodynamic stability (25,007/25,007 proteins, 0 failures). 60 statistical tests, 28 significant. Manuscript in preparation.

---

## Table of Contents

1. [Project Overview and Motivation](#1-project-overview-and-motivation)
2. [The Three Scientific Goals](#2-the-three-scientific-goals)
3. [The Seven Datasets](#3-the-seven-datasets)
4. [Critical Scientific Decisions and Rationale](#4-critical-scientific-decisions-and-rationale)
5. [Pre-Registered Hypotheses](#5-pre-registered-hypotheses)
6. [Project Directory Structure](#6-project-directory-structure)
7. [Phase 1: Pilot-Scale Analysis (Local Mac)](#7-phase-1-pilot-scale-analysis-local-mac)
   - [Module A: Data Acquisition and Cleaning](#module-a-data-acquisition-and-cleaning)
   - [Module B: Dataset Construction](#module-b-dataset-construction)
   - [Module C: Orthology and Homology](#module-c-orthology-and-homology)
   - [Module D: Structure Acquisition and Quality Control](#module-d-structure-acquisition-and-quality-control)
   - [Module E: Structural Domain Assignment](#module-e-structural-domain-assignment)
   - [Module F: N-vs-C Terminus Stability Analysis](#module-f-n-vs-c-terminus-stability-analysis)
   - [Module G: Mitochondrial Targeting Signal Analysis](#module-g-mitochondrial-targeting-signal-analysis)
   - [Module H: Comparative Statistics](#module-h-comparative-statistics)
   - [Module I: Publication Figures](#module-i-publication-figures)
8. [Phase 2: Full-Scale HPC Pipeline](#8-phase-2-full-scale-hpc-pipeline)
   - [HPC Infrastructure](#hpc-infrastructure)
   - [Pipeline Architecture](#pipeline-architecture)
   - [Phase 2 Directory Structure (Two `phase2` Folders)](#phase-2-directory-structure-two-phase2-folders)
   - [Step-by-Step HPC Execution](#step-by-step-hpc-execution)
   - [FoldX Thermodynamic Stability](#foldx-thermodynamic-stability)
9. [Results and Interpretations](#9-results-and-interpretations)
   - [Goal 1: Domain Architecture](#goal-1-domain-architecture-results)
   - [Goal 2: N-vs-C Structural Asymmetry](#goal-2-n-vs-c-structural-asymmetry-results)
   - [Goal 3: Mitochondrial Targeting](#goal-3-mitochondrial-targeting-results)
   - [Cross-Species Conservation](#cross-species-conservation-results)
   - [Statistical Summary](#statistical-summary)
10. [Biological Synthesis and Interpretation](#10-biological-synthesis-and-interpretation)
11. [Limitations and Caveats](#11-limitations-and-caveats)
12. [Software and Reproducibility](#12-software-and-reproducibility)
13. [Session History](#13-session-history)
14. [Appendix: Complete File Inventory](#14-appendix-complete-file-inventory)

---

## 1. Project Overview and Motivation

### What is this project about?

Every living cell contains thousands of proteins that must fold into precise three-dimensional structures to function. Most proteins fold spontaneously, but a significant fraction — estimated at 10-15% of the proteome — requires assistance from molecular machines called **chaperonins**. These are large barrel-shaped protein complexes that encapsulate unfolded or misfolded proteins and provide a protected environment for them to fold correctly.

Two chaperonin systems are the focus of this study:

1. **GroEL/GroES** in *Escherichia coli* (bacteria): The best-studied chaperonin. GroEL is a ~800 kDa double-ring complex that works with its co-chaperonin GroES. It has been shown to interact with approximately 250 substrate proteins, classified by Kerner et al. (2005) into three dependency classes:
   - **Class I** (38 proteins): Can fold spontaneously but are accelerated by GroEL
   - **Class II** (126 proteins): Partially dependent — fold slowly without GroEL
   - **Class III** (84 proteins): Obligate substrates — cannot fold without GroEL and aggregate in its absence

2. **HSP60 (HSPD1)** in humans: The mitochondrial homolog of GroEL. It resides in the mitochondrial matrix and assists in folding proteins that are imported from the cytoplasm after translation on cytoplasmic ribosomes. HSP60 works with its co-chaperonin HSP10 (HSPE1).

### Why compare them?

GroEL and HSP60 are evolutionary homologs separated by approximately **2 billion years** of evolution — since the endosymbiotic event that gave rise to mitochondria from an alpha-proteobacterial ancestor. Despite this vast evolutionary distance, both perform the same fundamental function: assisting protein folding. This creates a natural experiment:

- **If substrate properties are conserved**: It implies strong evolutionary constraints on what makes a protein require chaperonin assistance — the structural features that make folding difficult are so fundamental that they have been maintained across 2 billion years.
- **If substrate properties have diverged**: It reveals how the chaperonin-substrate relationship has co-evolved with the transition from a free-living bacterium to an organelle within a eukaryotic cell.

### The project name

"Antah Asti Prarambh" is Sanskrit for "The End is the Beginning." This refers to the central scientific question: **Does the C-terminus (the "end" of translation) represent a new "beginning" for chaperonin-assisted folding?** Proteins are synthesized from N-terminus to C-terminus on the ribosome. The N-terminal region emerges first, folds first, and may already be structured before the C-terminus is even synthesized. This creates an inherent asymmetry: the C-terminal region is the last to emerge and the last to fold — potentially the region most in need of chaperonin assistance. Our project tests whether this asymmetry is a universal feature of protein architecture or specific to chaperonin substrates.

---

## 2. The Three Scientific Goals

### Goal 1: Structural Domain Distribution

**Question:** Do chaperonin substrates have distinctive structural domain architectures compared to the rest of the proteome?

**Rationale:** Not all protein folds are equally difficult to achieve. TIM barrels (beta/alpha-barrel topology), for example, have a complex topology that requires multiple strands and helices to come together in a specific arrangement. If chaperonins preferentially assist proteins with specific fold types, it reveals the structural basis of chaperonin substrate selectivity.

**Approach:** Assign structural domains to all proteins using CATH (Class, Architecture, Topology, Homologous superfamily) classification. Compare the distribution of domain superfamilies between chaperonin substrates and size-matched background proteomes using Fisher's exact test with odds ratios.

### Goal 2: N-terminal vs C-terminal Stability Asymmetry

**Question:** Are N-terminal structural domains more complex (higher contact order) than C-terminal regions? Is this asymmetry greater in chaperonin substrates?

**Rationale:** During co-translational folding, the N-terminal region of a protein emerges from the ribosome first. If N-terminal domains adopt complex folds (high contact order = many long-range contacts), they may be more difficult to fold and more likely to misfold during synthesis. Chaperonins could have evolved to specifically assist proteins where this N-terminal folding challenge is most severe.

**Key concept — Contact Order:** Contact order (Plaxco et al., 1998) measures the average sequence distance between residues that are in physical contact in the folded protein. High contact order means the fold requires residues that are far apart in sequence to come together — these folds are topologically complex and fold more slowly. Relative contact order (RCO) normalizes by protein length.

**Three-region decomposition:**
- **Pre-domain tail**: Residues 1 to the start of the first structural domain (includes signal/transit peptides if present)
- **N-domain**: The first structural domain (defined by CATH or Chainsaw boundaries)
- **C-region**: Everything after the first domain ends (may contain additional domains)

### Goal 3: Mitochondrial Matrix Targeting and MTS Architecture

**Question:** How do mitochondrial targeting signals (MTS) relate spatially to the first structural domain in HSP60 substrates?

**Rationale:** HSP60 substrates are nuclear-encoded proteins that must be imported into the mitochondrial matrix. Most carry an N-terminal mitochondrial targeting signal (MTS, also called transit peptide) that is cleaved upon import. After cleavage, the protein emerges into the matrix in an unfolded state and encounters HSP60. The spatial relationship between the MTS cleavage site and the first structural domain determines whether the protein arrives at HSP60 with its N-terminal domain already partially folded or completely unfolded.

**Key question:** Is the MTS a separate "pre-domain" extension (cleaved off before the first domain, leaving the first domain to emerge unfolded into the matrix), or does it overlap with the first structural domain (meaning domain integrity is disrupted by cleavage)?

---

## 3. The Seven Datasets

The project assembles seven carefully curated datasets that together enable all three scientific goals:

### Dataset 1: *E. coli* K-12 Proteome (Background for GroEL)

| Property | Value |
|----------|-------|
| **Source** | UniProt reference proteome UP000000625 |
| **Size** | 4,403 proteins |
| **Files** | `data/raw/uniprot/ecoli_k12_proteome.fasta`, `data/raw/uniprot/ecoli_k12_proteome.tsv` |
| **Purpose** | Background/control for GroEL substrate comparisons |

This is the complete set of reviewed proteins in *E. coli* strain K-12 (substrain MG1655). It serves as the denominator for all GroEL enrichment calculations: "Is superfamily X more common among GroEL substrates than expected by chance given its frequency in the full proteome?"

### Dataset 2: Human Proteome (Background for HSP60)

| Property | Value |
|----------|-------|
| **Source** | UniProt reference proteome UP000005640 |
| **Size** | 20,416 proteins |
| **Files** | `data/raw/uniprot/human_proteome.fasta`, `data/raw/uniprot/human_proteome.tsv` |
| **Purpose** | Background for orthology analysis; parent set for Datasets 3, 5, 7 |

The canonical human proteome (one protein per gene). Used for OrthoFinder-based orthology detection and as the outermost background layer for HSP60 comparisons.

### Dataset 3: Human Mitochondrial Proteome (Compartment-Matched Background)

| Property | Value |
|----------|-------|
| **Source** | MitoCarta 3.0 (Broad Institute) |
| **Size** | 1,136 proteins |
| **File** | `data/processed/human_mito_proteome.tsv` |
| **Purpose** | Compartment-matched background for HSP60 analysis |

**Why compartment-matched?** HSP60 resides in the mitochondrial matrix. Comparing HSP60 substrates to the full human proteome would conflate "properties of mitochondrial proteins" with "properties of chaperonin substrates." By using the mitochondrial proteome as background, we control for the general properties of mitochondrial proteins and isolate features specific to chaperonin dependence.

**MitoCarta 3.0 details:** MitoCarta is a curated inventory of the human mitochondrial proteome based on mass spectrometry of purified mitochondria, GFP localization, and literature curation. Version 3.0 includes sub-compartment localization: matrix, inner membrane (MIM), intermembrane space (IMS), and outer membrane (MOM). We use v3.0 specifically because v2.0→v3.0 involved 70 localization reclassifications (e.g., 52 respiratory chain subunits moved from "matrix" to "inner membrane"), which significantly affects our analyses.

### Dataset 4: GroEL Substrates

| Property | Value |
|----------|-------|
| **Source** | Kerner et al., Cell 2005 (Table S3) |
| **Size** | 252 proteins |
| **Files** | `data/processed/groel_substrates_standardized.tsv`, `data/processed/groel_substrates.fasta` |
| **Purpose** | Primary substrate set for GroEL analyses |

**Data cleaning challenges:** The Kerner 2005 dataset used UniProt accessions from 2005. In the intervening 21 years, 149 of 252 accessions (59%) had been "demerged" — a single accession was split into multiple strain-specific entries. Our validation script (`scripts/validate_uniprot_accessions.py`) queries the current UniProt API, identifies all demerged targets, and selects the K-12 (strain MG1655, taxon 83333) specific entry, preferring reviewed Swiss-Prot entries over unreviewed TrEMBL. Four accessions were plasmid-specific and not present in the K-12 reference proteome.

**GroEL dependency classes:**
- **Class I** (38 proteins): Spontaneous folders; interact with GroEL but do not require it. These are "opportunistic" substrates.
- **Class II** (126 proteins): Partially dependent; fold slowly or incompletely without GroEL.
- **Class III** (84 proteins): Obligate substrates; aggregate completely without GroEL. These are the proteins that *require* the chaperonin barrel.

### Dataset 5: HSP60 Tier-1 Substrates

| Property | Value |
|----------|-------|
| **Source** | Morten et al., Molecular Cell 2020 (co-immunoprecipitation + SILAC) |
| **Size** | 266 proteins (from 325 raw) |
| **Files** | `data/processed/hsp60_tier1_substrates.tsv`, `data/processed/hsp60_tier1_substrates.fasta` |
| **Purpose** | Primary substrate set for HSP60 analyses |

**Why careful filtering is critical:** The raw dataset comes from co-immunoprecipitation (co-IP) of HSP60 followed by mass spectrometry. Co-IP captures any protein that physically associates with HSP60 — not only true substrates but also co-chaperones, bystander proteins in the matrix, and contaminants. SILAC (Stable Isotope Labeling by Amino acids in Cell culture) quantification provides an enrichment ratio: how much more of each protein is found in the HSP60 pull-down versus the control.

**Filtering steps (script: `scripts/filter_hsp60_interactome.py`):**
1. **Exclude baits** (2): HSPD1 (HSP60 itself), HSPE1 (HSP10 co-chaperonin)
2. **Exclude co-chaperones** (4): TRAP1, HSPA9 (mortalin), GRPEL1, DNAJA3 — these are part of the chaperone machinery, not substrates
3. **Exclude contaminants** (4): Immunoglobulins, keratins, tubulins — common co-IP contaminants
4. **SILAC filtering**: Require median imputed SILAC ratio > 5 (stringent enrichment over control)
5. **MitoCarta confirmation**: Require MitoCarta 3.0 annotation for Tier 1

**NDIC handling:** "NDIC" (Not Detected In Control) means the protein was found only in the HSP60 pull-down and not in the control. This is paradoxically the *strongest* evidence of enrichment, not missing data. NDIC values are imputed as 2× the 95th percentile of observed SILAC ratios per replicate.

**Quality tiers:**
- **Tier 1** (266 proteins): MitoCarta-confirmed + median SILAC > 5 — high-confidence substrates
- **Tier 2**: Candidate substrates with median SILAC 2-5 — moderate confidence
- **Tier 3**: Candidate substrates with SILAC ≤ 2 — low confidence
- **Excluded**: Baits, co-chaperones, contaminants, non-candidate

### Dataset 6: Cross-Species Homolog Pairs

| Property | Value |
|----------|-------|
| **Source** | Computed from OrthoFinder + MMseqs2 RBH |
| **Size** | 69 pairs |
| **File** | `data/processed/groel_hsp60_homologs.tsv` |
| **Purpose** | Test conservation of substrate properties across species |

This dataset identifies GroEL substrates that have homologs among HSP60 substrates. A protein pair in Dataset 6 means: "This *E. coli* protein is a chaperonin substrate, its human ortholog is also a chaperonin substrate, and both are structural/functional counterparts."

**Two complementary methods:**
1. **Reciprocal Best Hit (RBH)**: MMseqs2 easy-rbh identifies 40 pairs where each protein is the other's best match. Simple but misses many-to-many relationships.
2. **OrthoFinder orthogroups**: Connected-component clustering on bidirectional MMseqs2 all-vs-all search identifies 422 orthogroups; 34 contain both GroEL and HSP60 substrates, yielding 62 substrate pairs.

**Merged result:** 69 unique pairs (33 found by both methods, 7 RBH-only, 29 orthogroup-only). The "evidence" column tracks how each pair was identified.

### Dataset 7: Mitochondrial Matrix-Only Proteome

| Property | Value |
|----------|-------|
| **Source** | MitoCarta 3.0, filtered to `is_matrix = 1` |
| **Size** | 525 proteins |
| **File** | `data/processed/human_matrix_proteome.tsv` |
| **Purpose** | Tightly compartment-matched background for HSP60 |

**Why a separate matrix set?** HSP60 physically resides in the matrix. Many mitochondrial proteins are in the inner membrane, intermembrane space, or outer membrane — they never encounter HSP60. Using only matrix proteins as background is the most stringent test: "Among proteins that are physically in the same compartment as HSP60, which ones actually require its assistance?"

---

## 4. Critical Scientific Decisions and Rationale

Nine key decisions were made at the project's outset, each with explicit reasoning documented in `analysis.md`:

### Decision 1: CATH/Chainsaw for Domain Boundaries (NOT InterPro)

**Decision:** Use CATH structural domains from Gene3D (primary) supplemented by Chainsaw ML predictions (secondary). Do NOT use InterPro as the primary source.

**Why:** InterPro provides sequence-derived domain boundaries based on HMM profiles (e.g., Pfam, SMART). These boundaries reflect sequence conservation, not structural domain boundaries. For a study about structural folding, we need boundaries that reflect where one independently folding unit ends and another begins. CATH domains are defined by structural classification of experimentally solved structures. Chainsaw (Wells et al., 2024) is a deep learning model trained on CATH domains that predicts structural boundaries from AlphaFold structures using STRIDE secondary structure assignments.

**Impact:** Using InterPro would have produced different (often incorrect) domain boundaries, leading to inaccurate N-domain vs C-region decompositions and potentially misleading stability comparisons.

### Decision 2: Contact Order + FoldX for Stability (NOT pLDDT Alone)

**Decision:** Use relative contact order (Plaxco et al., 1998) as the primary folding kinetics proxy, and FoldX computed ΔG (kcal/mol) as the thermodynamic stability measure. AlphaFold pLDDT is used as a model confidence metric only.

**Why this is CRITICAL:** pLDDT (predicted Local Distance Difference Test) is AlphaFold's per-residue confidence score. It measures how confident AlphaFold is about its prediction, NOT the thermodynamic stability of that protein region. A region with low pLDDT could be genuinely disordered OR it could be a stable domain that AlphaFold simply has insufficient data to model accurately. Conflating pLDDT with stability is a widespread error in the structural proteomics literature.

**Contact order** (CO) measures folding kinetics: high CO means the fold requires many long-range contacts (topologically complex, folds slowly). This is a structural property directly computed from 3D coordinates.

**FoldX** computes the empirical free energy of folding (ΔG) using a semi-empirical force field that accounts for van der Waals clashes, solvation, hydrogen bonding, electrostatics, and entropy. A more negative ΔG means a more thermodynamically stable fold.

Together, CO (folding difficulty) and ΔG (folded-state stability) provide a complete picture that pLDDT alone cannot.

### Decision 3: OrthoFinder on Full Proteomes (RBH as Supplementary)

**Decision:** Run orthology detection on the full *E. coli* (4,403) and human (20,416) proteomes, not just on the substrate sets. Use RBH as a supplementary cross-check.

**Why:** RBH between substrate sets (252 × 266) only finds one-to-one best matches. Many chaperonin substrates belong to large families (e.g., aminotransferases, dehydrogenases) where multiple paralogs exist. RBH misses many-to-many relationships. OrthoFinder uses graph-based clustering on the full proteome, correctly grouping paralogs and orthologs into orthogroups. This increased our substrate pair count from 40 (RBH) to 69 (merged).

### Decision 4: SILAC-Based HSP60 Filtering with MitoCarta Confirmation

**Decision:** Filter HSP60 interactome by SILAC enrichment ratio and MitoCarta membership. Do not use all co-IP hits as substrates.

**Why:** Co-immunoprecipitation captures physical interactions, not functional substrate relationships. Without SILAC quantification, abundant matrix proteins (e.g., ribosomal subunits, metabolic enzymes) that are not true substrates would contaminate the dataset. Requiring both SILAC enrichment > 5 AND MitoCarta annotation ensures high confidence. The percent_occupancy column from the original paper is retained as a continuous weight for sensitivity analyses.

### Decision 5: Compartment-Matched AND Size-Matched Controls

**Decision:** All enrichment tests use both compartment-matched AND size-matched background controls.

**Why compartment-matched:** Comparing HSP60 substrates to the full human proteome conflates "mitochondrial protein properties" with "chaperonin substrate properties." By comparing within the matrix proteome (Dataset 7), we isolate chaperonin-specific features.

**Why size-matched:** Larger proteins tend to have more domains, different amino acid compositions, and different fold distributions. Without size matching, enrichment of a domain superfamily among substrates could simply reflect the fact that substrates happen to be larger than average. We stratify by 10 kDa bins and draw size-matched controls from the appropriate proteome.

### Decision 6: Three-Region N-Terminal Decomposition

**Decision:** Decompose each protein into three regions: (1) pre-domain tail, (2) first structural domain (N-domain), (3) remainder (C-region). Do NOT conflate the literal N-terminus, the transit peptide, and the first structural domain.

**Why:** These three concepts are biologically distinct:
- The **pre-domain tail** (residues 1 to first_domain_start - 1) may include the transit peptide, unstructured linkers, or signal sequences
- The **N-domain** (first structural domain) is the first independently folding unit
- The **C-region** (everything after the first domain) contains all subsequent domains

Conflating these would mix transit peptide properties with domain stability measurements, producing biologically meaningless results.

### Decision 7: MitoCarta 3.0 as Ground Truth

**Decision:** Use MitoCarta 3.0 (not v2.0) as the authoritative source for mitochondrial localization. Use computational predictors (TargetP, DeepMito) only for proteins lacking MitoCarta annotation.

**Why:** MitoCarta is based on experimental evidence (mass spectrometry of purified mitochondria, GFP localization). Predictors have ~15-30% false positive/negative rates. Between v2.0 and v3.0, 70 proteins were reclassified — notably, 52 respiratory chain subunits moved from "matrix" to "inner membrane," which significantly changes our background set.

### Decision 8: Hierarchical Testing with Pre-Registered Hypotheses

**Decision:** Pre-register all hypotheses before analysis. Use hierarchical Benjamini-Hochberg correction: first correct within each hypothesis family, then apply a Simes test across families.

**Why:** With 56 tests in Phase 2 (281 in the Phase 1 pilot), uncorrected p-values would produce many false positives. Standard Bonferroni correction is too conservative. Hierarchical BH groups related tests into families (domain architecture, stability asymmetry, MTS targeting) and corrects within each family first, then adjusts across families. This preserves power for related tests while controlling the overall false discovery rate at 0.05.

### Decision 9: Structural Domain Boundaries from Chainsaw Use STRIDE (Not Bioconda "stride")

**Decision:** Use the heiniglab/stride protein secondary structure assignment program, NOT the bioconda package called "stride" (which is a genomic variant caller).

**Why:** This was a hard-won lesson. `conda install stride` installs a completely different tool (a genomics tool for structural variant calling). The STRIDE needed for Chainsaw is the protein secondary structure assignment program by Frishman and Argos. It must be compiled from source (heiniglab/stride on GitHub). Using the wrong tool produces zero output with no error message, wasting hours of debugging.

---

## 5. Pre-Registered Hypotheses

Nine hypotheses were pre-registered before Module H analysis, organized into three families:

### Family 1: Domain Architecture (H1)

| ID | Hypothesis | Test | Effect Size |
|----|-----------|------|-------------|
| H1.1 | GroEL substrates (especially Class III) enriched for specific CATH superfamilies (TIM barrels, Rossmann folds) | Fisher's exact test per superfamily | Odds ratio + 95% CI |
| H1.2 | HSP60 substrates show enrichment for specific fold architectures vs size-matched matrix background | Fisher's exact test per superfamily | Odds ratio + 95% CI |
| H1.3 | Structural fold enrichment conserved between GroEL and HSP60 (shared orthogroups show similar folds) | Chi-squared / Fisher's | Cramer's V |

### Family 2: N-vs-C Stability Asymmetry (H2)

| ID | Hypothesis | Test | Effect Size |
|----|-----------|------|-------------|
| H2.1 | N-terminal domains have different contact order than C-terminal regions (paired within-protein test) | Wilcoxon signed-rank | Rank-biserial r |
| H2.2 | N-vs-C asymmetry GREATER in chaperonin substrates than in background proteomes | Mann-Whitney U | Rank-biserial r |
| H2.3 | Class III obligate GroEL substrates show GREATER N-vs-C asymmetry than Class I spontaneous folders | Kruskal-Wallis H | Eta-squared |

### Family 3: MTS Targeting (H3)

| ID | Hypothesis | Test | Effect Size |
|----|-----------|------|-------------|
| H3.1 | HSP60 substrates enriched for matrix localization vs general mitochondrial proteome | Fisher's exact test | Odds ratio |
| H3.2 | MTS-bearing HSP60 substrates have distinct first-domain properties vs non-MTS substrates | Mann-Whitney U | Rank-biserial r |
| H3.3 | MTS is predominantly a pre-domain extension, NOT part of the first structural domain | Binomial test | Observed proportion |

---

## 6. Project Directory Structure

```
/Users/vishalbharti/Downloads/Antah_Asti_Prarambh/
│
├── config/                                # Reserved for configuration files
├── logs/                                  # Execution logs
│
├── data/                                  # All input data (~512 MB)
│   ├── raw/
│   │   ├── uniprot/                       # Downloaded proteomes
│   │   │   ├── ecoli_k12_proteome.fasta   # Dataset 1 (4,403 proteins)
│   │   │   ├── ecoli_k12_proteome.tsv
│   │   │   ├── human_proteome.fasta       # Dataset 2 (20,416 proteins)
│   │   │   └── human_proteome.tsv
│   │   ├── alphafold/
│   │   │   ├── pilot/                     # Phase 1: 1,382 CIF structures (~466 MB)
│   │   │   │   └── AF-{ACC}-F1-model_v{4|6}.cif
│   │   │   └── full/                      # Phase 2: on HPC only (~22 GB)
│   │   ├── mitocarta/
│   │   │   └── Human.MitoCarta3.0.xls     # MitoCarta 3.0 Excel
│   │   └── custom/                        # Original literature data
│   │
│   ├── interim/                           # Intermediate processing (empty)
│   │
│   └── processed/                         # Cleaned, standardized datasets
│       ├── groel_substrates_standardized.tsv     # Dataset 4 (252 proteins)
│       ├── groel_substrates.fasta                # Dataset 4 sequences
│       ├── hsp60_interactome_standardized.tsv    # Full HSP60 (all tiers)
│       ├── hsp60_tier1_substrates.tsv            # Dataset 5 (266 proteins)
│       ├── hsp60_tier1_substrates.fasta          # Dataset 5 sequences
│       ├── human_mito_proteome.tsv               # Dataset 3 (1,136 proteins)
│       ├── human_matrix_proteome.tsv             # Dataset 7 (525 proteins)
│       ├── groel_hsp60_homologs.tsv              # Dataset 6 (69 pairs)
│       ├── hsp60_filtering_report.txt
│       └── mitocarta_summary_report.txt
│
├── results/                               # All analysis outputs (~186 MB)
│   ├── homology/                          # Module C outputs
│   ├── structures/                        # Module D outputs
│   ├── domains/                           # Module E outputs
│   ├── termini/                           # Module F outputs
│   ├── mts/                               # Module G outputs
│   ├── stats/                             # Module H outputs (Phase 1)
│   ├── figures/                           # Module I outputs (Phase 1)
│   ├── dataset_build/                     # Intermediate dataset construction
│   │
│   └── phase2/                            # *** PHASE 2 RESULTS (transferred from HPC) ***
│       ├── structures/                    # Full-scale structure index
│       ├── domains/                       # Chainsaw + unified domains (25,258 proteins)
│       ├── stability/                     # N-vs-C paired, contact order (all proteins)
│       ├── stats/                         # Full-scale statistics (56 tests)
│       ├── figures/                       # Publication figures (6 × PDF+PNG)
│       ├── foldseek/                      # Structural clustering (16,193 clusters)
│       └── foldx/                         # FoldX results (COMPLETE — 25,007 proteins)
│
├── scripts/                               # Standalone analysis scripts (4 files)
│   ├── validate_uniprot_accessions.py     # GroEL accession demerging
│   ├── filter_hsp60_interactome.py        # HSP60 SILAC filtering
│   ├── module_c_extract_fasta.py          # FASTA extraction
│   └── module_c_analyze_rbh.py            # RBH analysis
│
├── workflow/                              # Pipeline code
│   ├── scripts/                           # Phase 1 module scripts (16 Python files)
│   │   ├── download_alphafold_pilot.py    # Module D
│   │   ├── run_dssp.py                    # Module D (DSSP)
│   │   ├── validate_structure_quality.py  # Module D (quality)
│   │   ├── get_cath_domains.py            # Module E (CATH)
│   │   ├── run_chainsaw_e2.py             # Module E (Chainsaw)
│   │   ├── compute_domain_structural_metrics.py  # Module E (metrics)
│   │   ├── analyze_foldseek.py            # Module E (clustering)
│   │   ├── domain_distribution_summary.py # Module E (summary)
│   │   ├── run_orthology.py               # Module C (orthogroups)
│   │   ├── build_dataset6_homologs.py     # Module C (Dataset 6)
│   │   ├── module_f_n_vs_c_analysis.py    # Module F (N-vs-C)
│   │   ├── module_f_extension_chainsaw.py # Module F (extend to Chainsaw)
│   │   ├── module_g_mts_analysis.py       # Module G (MTS)
│   │   ├── module_h_comparative_stats.py  # Module H (statistics)
│   │   ├── generate_figures.py            # Module I (figures)
│   │   └── parse_mitocarta.py             # Module A (MitoCarta)
│   │
│   └── phase2/                            # *** PHASE 2 PIPELINE CODE ***
│       ├── Snakefile                      # Snakemake orchestration
│       ├── config.yaml                    # Central configuration
│       ├── README.md                      # Deployment guide
│       ├── download_alphafold_full.py     # Full AlphaFold download
│       ├── run_foldseek_full.py           # Full Foldseek clustering
│       ├── run_foldx.py                   # FoldX stability wrapper
│       ├── scripts/                       # Analysis modules
│       │   ├── module_f_full.py           # Full-scale N-vs-C
│       │   ├── module_h_full.py           # Full-scale statistics
│       │   ├── module_i_full.py           # Full-scale figures
│       │   └── module_i_polished.py       # Polished figures
│       ├── slurm_jobs/                    # 19 SLURM job scripts
│       │   ├── 00_setup.sh → 11_module_i.sh
│       │   ├── submit_pipeline.sh
│       │   └── submit_analysis.sh
│       ├── envs/                          # Conda environments
│       └── rules/                         # Snakemake rules
│
├── docs/                                  # Documentation (16 files)
│   ├── COMPREHENSIVE_PROJECT_DOCUMENT.md  # This document
│   ├── COLLABORATOR_PRESENTATION.md
│   ├── DATA_HANDOVER_INDEX.md
│   ├── COLLABORATOR_SHARING_GUIDE.md
│   ├── PRIMARY_HYPOTHESES.md
│   ├── METHODS_AND_PROTOCOLS.md
│   ├── PROJECT_PLAN.md
│   ├── RESULTS_NARRATIVE.md
│   ├── PHASE1_VERIFICATION.md
│   ├── PHASE2_RESULTS_REPORT.md
│   ├── HPC_DEPLOYMENT_GUIDE.md
│   ├── DOCUMENTATION.md
│   ├── SESSION_CONTINUITY.md
│   ├── SESSION4_DOCUMENTATION.md
│   ├── SESSION5_DOCUMENTATION.md
│   ├── SESSION6_DOCUMENTATION.md
│   └── TODO.md
│
├── hsp60_interactome_clean.tsv            # Raw HSP60 data (root)
├── kerner_2005_groel_interactors_clean.csv       # Raw GroEL data
├── kerner_2005_groel_interactors_table_s3.csv    # GroEL Table S3
├── 12192_2020_1080_MOESM4_ESM.xlsx        # Supplementary data
├── foldx5_Linux.zip                       # FoldX 5.1 binary
└── Quesries-strategies.docx               # Initial planning dialogue
```

### The Two `phase2` Folders — Explained

A common point of confusion is that the project contains **two separate `phase2` directories** that serve different purposes:

1. **`workflow/phase2/`** — This is the **pipeline code**. It contains all the scripts, SLURM job definitions, configuration files, and Snakemake workflow that define *how* the Phase 2 analysis is executed on the HPC cluster. Think of this as the "recipe."
   - Contents: `Snakefile`, `config.yaml`, `download_alphafold_full.py`, `run_foldseek_full.py`, `run_foldx.py`, `scripts/` (4 analysis modules), `slurm_jobs/` (19 SLURM scripts)

2. **`results/phase2/`** — This is the **output data**. It contains all the results that were produced by running the Phase 2 pipeline on the HPC cluster and then transferred back to the local Mac. Think of this as the "finished dish."
   - Contents: `structures/` (index), `domains/` (assignments), `stability/` (N-vs-C), `stats/` (p-values), `figures/` (7 publication figures), `foldseek/` (clusters), `foldx/` (COMPLETE — 25,007 proteins)

---

## 7. Phase 1: Pilot-Scale Analysis (Local Mac)

Phase 1 was executed entirely on a local Apple M1 MacBook Air (8 GB RAM) and analyzed the core substrate sets: 252 GroEL substrates, 266 HSP60 substrates, 1,136 mitochondrial proteins, and 525 matrix proteins — totaling 1,390 unique proteins. This pilot served to:
1. Validate all computational methods before scaling to full proteomes
2. Produce preliminary results to verify biological plausibility
3. Debug data integration issues (column names, accession formats, file paths)

### Module A: Data Acquisition and Cleaning

**Purpose:** Standardize all input datasets to a common format with current UniProt accessions.

#### Step A1: GroEL Substrate Standardization

**Script:** `scripts/validate_uniprot_accessions.py`

**What it does:** Takes the 252 GroEL substrates from Kerner et al. (2005) and validates every accession against the current UniProt database. Since the paper is from 2005, many accessions have been "demerged" (one old accession split into multiple strain-specific entries).

**Algorithm:**
1. Load 252 accessions from the cleaned CSV
2. For each accession, query the UniProt REST API (`/uniprotkb/{accession}`)
3. If ACTIVE: record current metadata (gene name, length, organism, etc.)
4. If DEMERGED: fetch all target accessions, filter to *E. coli* K-12 (taxon 83333), prefer reviewed Swiss-Prot entries
5. If MERGED: follow to the new accession
6. If DELETED/NOT_FOUND: flag for manual review
7. Extract SCOP fold codes from Table S3 raw text using regex (`[a-g]\.\d+`)
8. Categorize subcellular locations (cytoplasmic, periplasmic, etc.)
9. Write standardized TSV with 16 columns

**Result:** 252 proteins standardized. 149 demerged accessions resolved. 4 plasmid-specific proteins flagged. All proteins mapped to current K-12 entries.

**Output columns:** `original_accession`, `current_accession`, `accession_status`, `entry_name`, `gene_symbol`, `protein_name`, `organism`, `length`, `reviewed`, `groel_class`, `mass_kDa`, `subcellular_location`, `location_category`, `scop_folds`, `description_clean`

#### Step A2: HSP60 Interactome Filtering

**Script:** `scripts/filter_hsp60_interactome.py`

**What it does:** Takes the raw HSP60 co-IP dataset (325 proteins) and applies SILAC-based quality filtering to identify high-confidence substrates.

**Algorithm:**
1. Parse 3 SILAC ratio columns (co-IP replicates 1, 2, 3)
2. Handle NDIC values: impute as 2× the 95th percentile of observed values in each replicate (NDIC = "Not Detected In Control" = very high enrichment)
3. Compute per-protein metrics: n_silac_replicates, n_ndic_replicates, median_silac_imputed
4. Flag proteins: bait (HSPD1/HSPE1), co-chaperones (TRAP1, HSPA9, GRPEL1, DNAJA3), contaminants (immunoglobulins, keratins, tubulins)
5. Assign quality tiers: Tier 1 (266), Tier 2, Tier 3, Excluded
6. Generate comprehensive filtering report

**Result:** 325 → 266 Tier-1 substrates. Median SILAC enrichment ratio: 22.2 (very high confidence).

#### Step A3: MitoCarta Parsing

**Script:** `workflow/scripts/parse_mitocarta.py`

**What it does:** Parses the MitoCarta 3.0 Excel file to extract the human mitochondrial proteome (Dataset 3) and the matrix-only subset (Dataset 7).

**Algorithm:**
1. Read sheet "A Human MitoCarta3.0" from the Excel file
2. Extract: UniProt ID, gene symbol, MitoCarta score, sub-mitochondrial localization text
3. Create binary flags: `is_matrix` (text contains "Matrix"), `is_im` (inner membrane), `is_ims` (intermembrane space), `is_om` (outer membrane)
4. Write full mito proteome (1,136 proteins) and matrix subset (525 proteins)
5. Cross-reference with HSP60 interactome for MitoCarta v2.0→v3.0 localization changes

**Result:** 1,136 mitochondrial proteins (Dataset 3), 525 matrix proteins (Dataset 7). Documented 70 localization changes between MitoCarta versions.

### Module B: Dataset Construction

**Purpose:** Download reference proteomes from UniProt and extract FASTA sequences for all substrate proteins.

#### Step B1: UniProt Proteome Download

Downloaded from UniProt FTP:
- *E. coli* K-12 (UP000000625): 4,403 reviewed proteins → `data/raw/uniprot/ecoli_k12_proteome.fasta`
- Human (UP000005640): 20,416 canonical proteins → `data/raw/uniprot/human_proteome.fasta`

#### Step B2: FASTA Extraction

**Script:** `scripts/module_c_extract_fasta.py`

**What it does:** Extracts sequences for all GroEL and HSP60 substrates from their respective proteome FASTA files, with fallback to UniProt REST API for missing entries.

**Algorithm:**
1. Index proteome FASTA by accession (parse `sp|ACC|ENTRY` headers)
2. For each substrate accession, look up in index
3. If missing: fetch from UniProt REST API with retry logic (3 attempts, exponential backoff)
4. Write output FASTA in the same order as the input TSV

**Result:** `groel_substrates.fasta` (252 sequences), `hsp60_tier1_substrates.fasta` (266 sequences).

### Module C: Orthology and Homology

**Purpose:** Identify evolutionary relationships between GroEL and HSP60 substrates to enable cross-species comparisons.

#### Step C1: Reciprocal Best Hit (RBH)

**Tool:** MMseqs2 `easy-rbh`

**What it does:** Finds pairs of proteins where each is the best BLAST hit of the other — the simplest orthology criterion.

**Command:** `mmseqs easy-rbh groel_substrates.fasta hsp60_tier1_substrates.fasta rbh_output tmp --format-output "query,target,evalue,bits,alnlen,qcov,tcov,pident,qlen,tlen"`

**Result:** 40 RBH pairs. Median sequence identity: 35.8%. E-value range: 2.3×10⁻¹⁵⁶ to 7.8×10⁻³.

**Key observation:** Only 8/84 (9.5%) Class III obligate GroEL substrates have HSP60 orthologs by RBH, compared to 18.4% of Class I and 19.0% of Class II. This suggests that the most GroEL-dependent proteins are the most divergent from their human counterparts.

#### Step C2: RBH Annotation

**Script:** `scripts/module_c_analyze_rbh.py`

Merges RBH results with GroEL class annotations and HSP60 metadata. Generates a comprehensive report with e-value distributions, identity distributions, coverage statistics, and top 20 pairs by bit score.

#### Step C3: OrthoFinder-Style Orthology

**Script:** `workflow/scripts/run_orthology.py`

**What it does:** Performs bidirectional all-vs-all MMseqs2 search on the full proteomes (4,403 × 20,416) and clusters reciprocal hits into orthogroups using connected components.

**Algorithm:**
1. Create MMseqs2 databases for both proteomes
2. Run bidirectional search (E. coli → Human, Human → E. coli) with e-value ≤ 1e-5, coverage ≥ 50%, pident ≥ 25%
3. Identify reciprocal pairs (found in both directions)
4. Build a graph where nodes are proteins and edges are reciprocal hits
5. Find connected components (union-find algorithm) — each component is an orthogroup
6. A single orthogroup can contain multiple E. coli and multiple human proteins (capturing paralogs)
7. Intersect with substrate lists to find shared substrate orthogroups

**Result:** 422 orthogroups total. 34 contain both GroEL and HSP60 substrates (= "shared"). These 34 orthogroups yield 62 substrate-substrate pairs.

#### Step C4: Dataset 6 Construction

**Script:** `workflow/scripts/build_dataset6_homologs.py`

**What it does:** Merges RBH pairs (40) with orthogroup pairs (62) into a unified Dataset 6.

**Algorithm:**
1. Expand orthogroup many-to-many relationships into all pairwise combinations
2. Match against RBH pairs to identify overlap
3. Assign evidence: "both" (33 pairs), "rbh_only" (7), "orthogroup_only" (29)
4. Sort by evidence priority, then by e-value

**Result:** 69 unique cross-species pairs. 79.7% share the same top CATH superfamily.

### Module D: Structure Acquisition and Quality Control

**Purpose:** Download AlphaFold predicted structures for all proteins and assess their quality.

#### Step D1-D2: AlphaFold Download

**Script:** `workflow/scripts/download_alphafold_pilot.py`

**What it does:** Downloads CIF-format AlphaFold structures from the EBI database for all pilot proteins.

**Algorithm:**
1. Collect unique accessions across all 4 datasets (GroEL, HSP60, mito, matrix)
2. Download from `https://alphafold.ebi.ac.uk/files/AF-{ACC}-F1-model_v6.cif` (fallback to v4)
3. Polite batching: 0.5s delay between batches of 50 requests
4. Parse each CIF file to extract: residue count, CA B-factors (= pLDDT scores), min/mean/frac_gt70 pLDDT
5. Build structure index TSV

**Result:** 1,382/1,390 structures downloaded (99.4%). 8 proteins lack AlphaFold structures: P07203, P30042, P36969, Q16881, Q5THJ4, Q86UA3, Q9BVL4, Q9NNW7. Mean pLDDT: 85.8. Total size: 466 MB.

#### Step D3: DSSP Secondary Structure Assignment

**Script:** `workflow/scripts/run_dssp.py`

**What it does:** Runs the mkdssp program on all AlphaFold structures to assign per-residue secondary structure (helix, strand, coil).

**DSSP codes:**
- **H** (α-helix), **G** (3₁₀-helix), **I** (π-helix) → grouped as "helix"
- **E** (β-strand), **B** (β-bridge) → grouped as "strand"
- **T** (turn), **S** (bend), **-** (coil) → grouped as "coil"

**Result:** 1,382/1,382 processed. Mean composition: 43.5% helix, 14.2% strand, 42.2% coil. Per-residue data stored in `dssp_per_residue.tsv` (used by Module F for regional analysis).

#### Step D4: Structure Quality Validation

**Script:** `workflow/scripts/validate_structure_quality.py`

**What it does:** Categorizes AlphaFold models into quality tiers and flags potentially unreliable structures.

**Quality tiers (by mean pLDDT):**
- **Very high** (≥ 90): Highly reliable backbone AND side chains
- **High** (80-90): Reliable backbone; side chains may have errors
- **Moderate** (70-80): Backbone generally correct; use with caution
- **Low** (50-70): Only fold-level topology is reliable
- **Very low** (< 50): Likely disordered or unreliable

**Flags:**
1. `flag_very_low_plddt`: mean pLDDT < 50 (structure essentially meaningless)
2. `flag_majority_unreliable`: mean pLDDT < 70 AND frac_plddt_gt70 < 0.5
3. `flag_few_usable`: frac_plddt_gt70 < 0.3

**Result:** 77.4% high or very high quality. 63 proteins flagged (4.6%). Core substrate datasets (GroEL, HSP60) have excellent quality: mean pLDDT 92-93.

### Module E: Structural Domain Assignment

**Purpose:** Assign structural domain boundaries to all proteins, combining curated CATH data with ML predictions, and perform structural clustering.

#### Step E1: CATH Domain Assignment via Gene3D

**Script:** `workflow/scripts/get_cath_domains.py`

**What it does:** Queries the InterPro REST API for Gene3D (CATH) domain annotations for all pilot proteins.

**CATH classification hierarchy:**
- **Class** (C): All-alpha, all-beta, alpha-beta, few-secondary-structures
- **Architecture** (A): Overall shape (e.g., barrel, sandwich, roll)
- **Topology** (T): Connectivity of secondary structure elements (e.g., TIM barrel = 3.20.20)
- **Homologous superfamily** (H): Evolutionary relationship (e.g., 3.20.20.70 = Aldolase class I)

**Algorithm:**
1. Query InterPro API per accession with rate limiting (1 req/sec)
2. Parse Gene3D entries from InterPro JSON
3. Handle discontinuous domains (fragments) — some domains consist of multiple non-contiguous segments
4. Checkpoint every 50 proteins for resume capability
5. Build protein summary: domain count, architecture string (CATH codes in positional order)

**Result:** 1,151/1,390 proteins (82.8%) have CATH domain assignments. 2,141 total domains. Mean pLDDT of assigned domains: 92.1 (higher than overall mean, because well-structured regions are easier to classify).

#### Step E2: Chainsaw ML Predictions for Remaining Proteins

**Script:** `workflow/scripts/run_chainsaw_e2.py`

**What it does:** Runs Chainsaw (Wells et al., 2024), a deep learning model trained on CATH domain structures, to predict domain boundaries for the 239 proteins that lack CATH annotations.

**Why Chainsaw?** Chainsaw predicts structural domain boundaries directly from AlphaFold 3D coordinates using STRIDE secondary structure assignments. Unlike sequence-based methods (HMMs), it captures structural features that define where one independently folding unit ends and another begins. It was specifically trained on CATH domains, ensuring consistency with our Gene3D annotations.

**Algorithm:**
1. Identify 239 proteins without CATH, filter to those with AlphaFold structures (236)
2. Symlink CIF files to a temporary directory (Chainsaw expects a directory of structures)
3. Run Chainsaw batch pipeline (~40 seconds per protein)
4. Parse chopping strings (e.g., "1-100_150-300" = two-segment domain) into domain records
5. Merge with CATH: proteins with CATH keep their Gene3D annotations; Chainsaw fills gaps
6. Write unified domain assignment table

**Result:** 236 additional proteins assigned. Unified coverage: 1,387/1,390 (99.8%). Only 3 proteins remain unassigned (no AlphaFold structure or Chainsaw failure).

#### Step E3: Per-Domain Structural Metrics

**Script:** `workflow/scripts/compute_domain_structural_metrics.py`

**What it does:** Computes secondary structure fractions and pLDDT confidence for each individual domain (not whole-protein).

**Algorithm:**
1. For each domain (defined by start/end residues), extract the corresponding DSSP secondary structure codes
2. From the CIF file, extract CA B-factors (pLDDT) for the domain residues
3. Compute: frac_helix, frac_strand, frac_coil, mean_pLDDT, min_pLDDT, frac_plddt_gt70, frac_plddt_gt90

**Result:** 2,141 domain-level records. Mean domain pLDDT: 92.1 (domains are better-predicted than full proteins because they exclude disordered regions).

#### Step E4: Foldseek Structural Clustering

**Tool:** Foldseek (van Kempen et al., 2024)
**Script:** `workflow/scripts/analyze_foldseek.py`

**What it does:** Performs all-vs-all structural similarity search using Foldseek's 3Di+AA alphabet, then clusters similar structures.

**Why Foldseek?** Traditional sequence-based clustering (e.g., CD-HIT) groups proteins by sequence identity. But proteins can have very different sequences yet identical 3D structures (< 20% sequence identity). Foldseek uses a structural alphabet (3Di) that encodes the local 3D arrangement of residues, allowing structural comparison at BLAST-like speed.

**Pipeline:**
1. `foldseek createdb` — Build structure database from CIF files
2. `foldseek search` — All-vs-all structural alignment (sensitivity 7.5, E-value 0.001)
3. `foldseek cluster` — Set-cover clustering on search results
4. `foldseek createtsv` — Export cluster membership to TSV

**Result:** 1,155 structural clusters from 1,382 proteins. 24 clusters are "shared" (contain at least one GroEL substrate AND one HSP60 substrate), meaning the same structural fold is a chaperonin substrate in both organisms.

#### Step E5: Domain Distribution Summary

**Script:** `workflow/scripts/domain_distribution_summary.py`

**What it does:** Summarizes domain architecture across all four datasets (GroEL, HSP60, matrix, mito).

**Key findings from Phase 1 pilot:**

| Metric | GroEL | HSP60 | Matrix | Mito |
|--------|-------|-------|--------|------|
| N proteins | 252 | 266 | 524 | 1,132 |
| % with CATH | 98% | 90.6% | 89.9% | 79.3% |
| Mean domains/protein | 2.02 | 1.95 | 2.01 | 1.81 |
| % alpha-beta | 66.3% | 67.5% | 71.3% | 60.1% |
| % mainly-alpha | 19.8% | 21.9% | 17.2% | 28.2% |
| % mainly-beta | 11.6% | 8.5% | 9.3% | 8.8% |
| Mean domain pLDDT | 93.4 | 92.7 | 92.8 | 91.7 |

**Interpretation:** Alpha-beta proteins dominate all datasets (60-71%), consistent with the general proteome. Chaperonin substrates have slightly higher domain pLDDT than background, suggesting they are well-structured proteins (not disordered). The differences between substrates and background are in specific superfamilies (tested in Module H), not broad fold classes.

### Module F: N-vs-C Terminus Stability Analysis

**Purpose:** Test whether N-terminal domains differ structurally from C-terminal regions, and whether this asymmetry is greater in chaperonin substrates.

#### Script: `workflow/scripts/module_f_n_vs_c_analysis.py`

This is the core analysis module addressing Goal 2. It implements a five-step pipeline:

#### Step F1: Region Definition

For each multi-domain protein:
1. Identify the **first structural domain** (N-domain) from CATH annotations, sorted by domain_start position
2. Define three regions:
   - **Pre-domain tail**: residues 1 to (first_domain_start - 1) — may include transit peptide, signal sequences, unstructured N-terminal extensions
   - **N-domain**: residues first_domain_start to first_domain_end — the first independently folding structural unit
   - **C-region**: residues (first_domain_end + 1) to protein_end — all subsequent domains and linkers

For single-domain proteins: the entire protein is treated as the N-domain; there is no C-region for paired comparison.

**Result:** 5,322 region boundaries defined.

#### Step F2: Sequence Metrics Per Region

For each region, compute from the amino acid sequence:
- **Net charge at pH 7**: Count of (K + R) - (D + E) per residue
- **Fraction hydrophobic**: (A, V, I, L, M, F, W, P) / total
- **Fraction charged**: (K, R, D, E) / total
- **Fraction polar**: (S, T, N, Q, Y, H, C) / total
- **Fraction aromatic**: (F, W, Y) / total

#### Step F3: Structure Metrics Per Region

For each region, extract from DSSP and AlphaFold CIF:
- **Secondary structure fractions** (helix, strand, coil) from DSSP
- **pLDDT statistics** (mean, min, frac_gt70, frac_gt90) from CIF B-factors

#### Step F4: Contact Order Calculation

**This is the most important metric.** Contact order is computed following the Plaxco-Simons definition:

```
Absolute Contact Order (ACO) = (1/N_contacts) × Σ |i - j|
```

where the sum is over all pairs (i, j) of CA atoms within 8 Å that are separated by at least 6 residues in sequence. The **relative contact order** (RCO) normalizes by the number of residues:

```
RCO = ACO / N_residues
```

High RCO means the fold requires many long-range contacts — it is topologically complex and folds slowly. The key insight is that RCO correlates strongly with experimental folding rate (Plaxco et al., 1998; r ≈ -0.75 for two-state folders).

**Implementation:** Parse CA coordinates from CIF files. Build distance matrix. Identify all CA-CA pairs within 8 Å with |i-j| ≥ 6. Compute ACO and RCO.

#### Step F5: Paired N-vs-C Comparisons

For each multi-domain protein with both N-domain and C-region:
1. Compute the difference: (N-domain metric) - (C-region metric) for each metric
2. Aggregate across all proteins in each dataset
3. Wilcoxon signed-rank test (paired, non-parametric) for each metric within each dataset
4. Rank-biserial effect size r = |W+ - W-| / (W+ + W-)

**Result:** 2,648 paired comparisons. The key finding: **N-domains have significantly HIGHER relative contact order than C-regions across ALL datasets** (not just substrates).

#### Extension to Chainsaw Proteins

**Script:** `workflow/scripts/module_f_extension_chainsaw.py`

Repeats the same F1-F5 analysis for the 236 proteins with Chainsaw domain boundaries, then merges with the original CATH results. This increases coverage from 1,151 to 1,387 proteins.

### Module G: Mitochondrial Targeting Signal Analysis

**Purpose:** Characterize how mitochondrial targeting signals relate spatially to structural domains.

**Script:** `workflow/scripts/module_g_mts_analysis.py`

#### Step G1-G2: UniProt Feature Extraction

Query UniProt REST API for transit peptide and signal peptide annotations:
- **Transit peptide** (TRANSIT 1..n): N-terminal extension cleaved upon mitochondrial import
- **Signal peptide** (SIGNAL 1..n): N-terminal extension for ER targeting (mutually exclusive with transit peptide in most cases)
- **Subcellular location**: Text description of known localization

Cache results in `uniprot_transit_signal_cache.tsv` to avoid re-querying.

#### Step G3: MitoCarta Integration

Build a lookup table from MitoCarta 3.0 with compartment flags (is_matrix, is_im, is_ims, is_om).

#### Step G4: Targeting Classification

Classify each protein into one of:
- **High-confidence matrix**: MitoCarta matrix + transit peptide annotation
- **Non-canonical matrix import**: MitoCarta matrix but NO transit/signal peptide (imported by alternative mechanisms)
- **Inner membrane / Outer membrane / IMS**: Other mitochondrial compartments
- **Secretory pathway**: Has signal peptide (not mitochondrial)
- **Non-mitochondrial**: No MitoCarta annotation, no transit peptide

#### Step G5: MTS-Domain Relationship

For proteins with both a transit peptide AND a CATH domain, compute:
- **gap_length** = first_domain_start - transit_peptide_end
- **mts_overlaps_domain**: transit peptide end > first domain start (overlap)
- **mts_is_pre_domain**: transit peptide ends before first domain starts (separate)

**Result:** 84.4% of transit peptides are pre-domain extensions (median gap: 18 residues). This means the MTS is cleaved off *before* the first structural domain, leaving the entire first domain to emerge unfolded into the matrix — exactly where HSP60 can engage it.

### Module H: Comparative Statistics

**Purpose:** Formal hypothesis testing across all three goals with rigorous multiple testing correction.

**Script:** `workflow/scripts/module_h_comparative_stats.py`

#### Goal 1 Tests: Domain Architecture Enrichment

For each CATH superfamily present in at least 3 substrate proteins:
1. Build a 2×2 contingency table:
   |  | Has superfamily | Lacks superfamily |
   |--|----------------|-------------------|
   | Substrates | a | b |
   | Background | c | d |
2. Fisher's exact test → odds ratio, 95% CI, p-value
3. Size-matched background: stratify by 10 kDa bins

**Phase 1 results:** GroEL enriched for 1.10.10.10 (Winged helix, OR=42.8, p=2.35×10⁻⁶), 3.20.20.70 (Aldolase/TIM barrel, OR=9.16, p=2.35×10⁻⁶). HSP60 shows depletion of 1.50.40.10 (mitochondrial carrier, OR=0.16).

#### Goal 2 Tests: N-vs-C Stability Asymmetry

Three levels of testing:
1. **Within-protein paired test** (H2.1): Wilcoxon signed-rank on N-domain RCO vs C-region RCO
2. **Substrate vs background** (H2.2): Mann-Whitney U comparing RCO difference (N-C) between substrates and background
3. **GroEL class comparison** (H2.3): Kruskal-Wallis H-test comparing RCO difference across Class I/II/III

#### Goal 3 Tests: MTS Targeting

1. **Matrix enrichment** (H3.1): Fisher's exact test, HSP60 substrates vs mito background for matrix localization
2. **MTS pre-domain dominance** (H3.3): Binomial test, proportion of MTS that are pre-domain extensions vs 50% null

#### Multiple Testing Correction

**Hierarchical Benjamini-Hochberg:**
1. Within each family (H1, H2, H3), apply BH correction to all tests
2. Compute a family-level summary p-value using the Simes method
3. Apply BH correction across the three family summary p-values
4. A test is "significant overall" only if both its within-family BH p-value < 0.05 AND its family is significant across families

**Phase 1 result:** 281 tests total, 22 significant after hierarchical correction.
**Phase 2 result:** 56 tests total, 25 significant after hierarchical correction.

### Module I: Publication Figures

**Purpose:** Generate 6 publication-quality figures summarizing all findings.

**Script:** `workflow/scripts/generate_figures.py` (Phase 1), `workflow/phase2/scripts/module_i_polished.py` (Phase 2)

**Design principles:**
- Colorblind-friendly palette (Okabe-Ito / Wong 2011)
- 300 DPI for publication quality
- Both PDF (vector) and PNG (raster) formats
- Real p-values and sample sizes annotated on all plots
- Consistent color coding: GroEL = blue, HSP60 = orange, Matrix = green, Mito = red

**The 6 figures:**

| Figure | Title | Content |
|--------|-------|---------|
| Fig 1 | Domain Architecture | CATH class distribution (stacked bars), top 10 superfamilies (horizontal bars), domain count distribution, TIM barrel enrichment visualization |
| Fig 2 | N-vs-C Stability | Split violin plots of relative contact order (N-domain vs C-region) by dataset; pLDDT violins; paired difference heatmap |
| Fig 3 | GroEL Class Comparison | Contact order and pLDDT difference by GroEL class (I/II/III); shows NO class effect (p=0.77) |
| Fig 4 | MTS Targeting | Sub-mitochondrial localization breakdown; MTS gap-to-first-domain histogram; transit peptide cleavage vs domain scatter |
| Fig 5 | Orthology & Conservation | Orthogroup categories (shared/GroEL-only/HSP60-only); N-domain RCO conservation scatter (r=0.84) |
| Fig 6 | Summary | Key findings overview; enriched fold heatmap; N-C asymmetry summary bars; MTS architecture pie |

---

## 8. Phase 2: Full-Scale HPC Pipeline

### Why Phase 2?

Phase 1 analyzed only the substrate proteins and their immediate mitochondrial backgrounds (1,390 proteins). Phase 2 scales to the **full proteomes** (4,403 E. coli + 20,416 human = ~25,000 proteins) to:
1. Compute structural domains for ALL proteins, not just substrates
2. Perform Foldseek structural clustering on the full set
3. Calculate FoldX thermodynamic stability for ALL proteins
4. Run statistical tests with properly powered backgrounds
5. Generate publication-quality figures with full-scale data

### HPC Infrastructure

**Cluster:** CSIR-IGIB HPC "Tejas" (tejas.igib.res.in)

| Resource | Value |
|----------|-------|
| Partition | `compute` |
| CPUs per node | 40 |
| RAM per node | 380 GB |
| QOS | `common` (max 300 GB, 5-day walltime, 5 concurrent jobs) |
| Filesystem | Lustre (`/lustre/vishal.bharti/`) |
| Scheduler | SLURM |
| Conda | Available (with `LD_LIBRARY_PATH` fix required) |

**Critical HPC lesson — LD_LIBRARY_PATH:** The HPC uses system gcc-8.4, whose `libstdc++.so` is older than what conda's numpy/scipy require. Every SLURM script must include:
```bash
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```
Without this, Python imports fail silently with obscure `ImportError` messages.

**Critical HPC lesson — Python buffering:** SLURM captures stdout/stderr at job completion, not in real-time. All scripts use `python3 -u` (unbuffered) to enable real-time log monitoring via `tail -f`.

### Pipeline Architecture

The Phase 2 pipeline has two execution modes:

1. **Snakemake** (automated): `Snakefile` + `config.yaml` define all rules and dependencies. Ideal for dry runs and reproducibility documentation.

2. **Manual SLURM submission** (used in practice): 19 SLURM scripts in `slurm_jobs/` with explicit dependency chains. This was preferred because:
   - The FoldX array job (501 tasks) required manual submission after generation
   - Debugging individual steps was easier with separate log files
   - The QOS limit (5 concurrent jobs) was easier to manage manually

**Dependency graph:**

```
00_setup.sh (validate environment)
    ↓
01_download_alphafold.sh (download 25,007 structures, ~22 GB)
    ↓
    ├── 02_foldseek_createdb.sh (build structure databases)
    │       ↓
    │   03_foldseek_search.sh (all-vs-all structural search, 64 GB RAM)
    │       ↓
    │   04_foldseek_cluster.sh (cluster similar structures)
    │
    ├── 05_chainsaw.sh (ML domain prediction, 72 hours)
    │       ↓
    │   08_module_e_domains.sh (unify CATH + Chainsaw)
    │
    └── 06_foldx_generate.sh (generate array job script)
            ↓
        [MANUAL] sbatch foldx_array.slurm (501 array tasks)
            ↓
        07_foldx_collect.sh (merge per-protein results)
            ↓
        09_module_f_stability.sh (N-vs-C with FoldX ΔG)
            ↓
        10_module_h_stats.sh (statistical testing)
            ↓
        11_module_i.sh (publication figures)
```

### Phase 2 Directory Structure (Two `phase2` Folders)

As noted in Section 6, the project has two `phase2` directories:

#### `workflow/phase2/` — The Pipeline Code

```
workflow/phase2/
├── Snakefile                      # Snakemake workflow definition
├── config.yaml                    # All paths, parameters, SLURM settings
├── README.md                      # Deployment guide
│
├── download_alphafold_full.py     # Bulk download from EBI FTP + individual fallback
├── run_foldseek_full.py           # createdb → search → cluster → export
├── run_foldx.py                   # CIF→PDB conversion + RepairPDB + Stability
│
├── scripts/
│   ├── module_f_full.py           # Three-region decomposition + contact order
│   ├── module_h_full.py           # Hierarchical BH statistical testing
│   ├── module_i_full.py           # Publication figures (6 panels)
│   └── module_i_polished.py       # Polished figures with real p-values
│
├── slurm_jobs/                    # 19 SLURM submission scripts
│   ├── 00_setup.sh                # Environment validation
│   ├── 01_download_alphafold.sh   # 4 GB RAM, 1 CPU, 4 hours
│   ├── 02_foldseek_createdb.sh    # 16 GB RAM, 4 CPUs, 2 hours
│   ├── 03_foldseek_search.sh      # 64 GB RAM, 16 CPUs, 24 hours
│   ├── 04_foldseek_cluster.sh     # 64 GB RAM, 16 CPUs, 8 hours
│   ├── 05_chainsaw.sh             # 16 GB RAM, 4 CPUs, 72 hours
│   ├── 05a_decompress_human.sh    # Decompress human AlphaFold tar
│   ├── 06_foldx_generate.sh       # 2 GB RAM, generates array script
│   ├── 07_foldx_collect.sh        # 4 GB RAM, merges per-protein JSONs
│   ├── 08_module_e_domains.sh     # 16 GB RAM, unified domain table
│   ├── 09_module_f.sh             # 16 GB RAM, N-vs-C stability
│   ├── 09_module_f_stability.sh   # Alternative stability wrapper
│   ├── 10_module_h.sh             # 16 GB RAM, statistics
│   ├── 10_module_h_stats.sh       # Alternative stats wrapper
│   ├── 11_module_i.sh             # 8 GB RAM, figures
│   ├── 11_module_i_figures.sh     # Alternative figures wrapper
│   ├── submit_pipeline.sh         # Master submission with dependencies
│   ├── submit_analysis.sh         # Analysis chain (F→H→I)
│   └── test_stride_fix.sh         # STRIDE binary validation
│
├── envs/                          # Conda environment definitions
└── rules/                         # Snakemake rule files
```

#### `results/phase2/` — The Output Data

```
results/phase2/
├── structures/
│   ├── structure_index_full.tsv           # 25,007 proteins with metadata
│   ├── download_checkpoint.json           # Download progress tracking
│   ├── missing_accessions_ecoli.txt       # 32 E. coli proteins without AF
│   └── missing_accessions_human.txt       # Proteins without AF
│
├── domains/
│   ├── chainsaw_full_predictions.tsv      # Chainsaw output (23,868 proteins)
│   ├── chainsaw_full_predictions_annotated.tsv  # With dataset labels
│   ├── unified_domain_assignments_full.tsv      # CATH + Chainsaw merged (25,258 records)
│   └── domain_distribution_full.tsv       # Domain counts by dataset
│
├── stability/
│   ├── region_boundaries_full.tsv         # Three-region boundaries (5,322 records)
│   ├── n_vs_c_paired_full.tsv             # Paired N-vs-C comparisons (2,648 pairs)
│   ├── contact_order_full.tsv             # RCO per region (11,824 records)
│   ├── sequence_metrics_full.tsv          # Charge, hydrophobicity per region
│   └── structure_metrics_full.tsv         # pLDDT, SS fractions per region
│
├── stats/
│   ├── corrected_pvalues_full.tsv         # 56 tests, 25 significant
│   ├── stability_comparisons_full.tsv     # All stability test details
│   └── statistics_summary_full.txt        # Human-readable report
│
├── figures/
│   ├── fig1_domain_architecture.{pdf,png}
│   ├── fig2_n_vs_c_stability.{pdf,png}
│   ├── fig3_groel_class_comparison.{pdf,png}
│   ├── fig4_mts_targeting.{pdf,png}
│   ├── fig5_orthology.{pdf,png}
│   └── fig6_summary.{pdf,png}
│
├── foldseek/
│   └── analysis/
│       ├── foldseek_clusters_full.tsv     # 16,193 clusters, 27,063 proteins
│       ├── combined_cluster_membership.tsv
│       └── foldseek_full_summary.txt
│
└── foldx/                                 # IN PROGRESS on HPC
    └── per_protein/                       # {ACC}.json files (10,775 done as of 2026-03-24)
```

### Step-by-Step HPC Execution

#### Step 0: Setup (`00_setup.sh`)
Validates the HPC environment: checks conda, Python packages, disk space, Phase 1 data availability. Creates output directories on Lustre.

#### Step 1: AlphaFold Download (`01_download_alphafold.sh` → `download_alphafold_full.py`)

Downloads all AlphaFold structures for both proteomes:
- **E. coli**: Bulk tar archive from EBI FTP (~2 GB) → extracts 4,371 CIF files
- **Human**: Bulk tar archive from EBI FTP (~20 GB) → extracts 20,636 CIF files
- **Fallback**: For any missing accessions, downloads individually from the API
- **Checkpoint**: JSON-based progress tracking allows resume after interruption

**Result:** 25,007 structures downloaded. Coverage: 99.3% of E. coli, ~100% of human.

#### Step 2-4: Foldseek Pipeline (`02-04_foldseek_*.sh` → `run_foldseek_full.py`)

All-vs-all structural similarity search and clustering:

**Step 2 — createdb** (16 GB, 4 CPUs, 2 hours):
Build Foldseek structure databases from CIF files (both organisms separately, then combined).

**Step 3 — search** (64 GB, 16 CPUs, 24 hours):
The most resource-intensive step. All-vs-all structural alignment using Foldseek's 3Di+AA alphabet.
- Sensitivity: 7.5 (high)
- E-value: 0.001
- Min coverage: 50%
- Split memory: 64 GB (prevents out-of-memory on large databases)

**Step 4 — cluster** (64 GB, 16 CPUs, 8 hours):
Set-cover clustering on search results, export to TSV, analyze cluster membership.

**Result:** 16,193 structural clusters from 27,063 proteins. 75.6% are singletons (unique structural folds). 25 clusters are "shared" (contain both GroEL and HSP60 substrates).

#### Step 5: Chainsaw Domain Prediction (`05_chainsaw.sh`)

Runs Chainsaw ML model on all 25,007 proteins in batches of 500:
- Requires STRIDE (heiniglab/stride, NOT bioconda stride) for secondary structure input
- Creates temporary symlink directories on Lustre for each batch
- Merges batch outputs into single TSV
- Wall time: 72 hours (4 CPUs, 16 GB)

**Result:** 25,007 proteins processed. 0 STRIDE failures (after the correct binary was installed). 93.6% assigned at least one domain.

#### Step 6-7: FoldX Thermodynamic Stability

**Step 6 — Generate array job** (`06_foldx_generate.sh` → `run_foldx.py --generate-slurm`):
Splits 25,007 proteins into 501 chunks of 50, generates a SLURM array job script.

**Step 6b — Submit array job** (MANUAL):
```bash
sbatch foldx_array.slurm  # Job 94152, 501 array tasks
```

**For each protein in each chunk:**
1. Convert CIF → PDB (using gemmi or BioPython)
2. Run FoldX RepairPDB (optimize rotamers, fix clashes)
3. Run FoldX Stability (compute ΔG at 298.15 K, pH 7.0, ionic strength 0.05)
4. Save result as JSON: `per_protein/{ACC}.json`

**FoldX 5.1 details:**
- Binary: `/lustre/vishal.bharti/software/foldx5/foldx` (symlink to `foldx_20270131`)
- License: Academic, expires 2026-12-31
- No rotabase.txt needed (compiled into binary)
- Output: `total_energy` (ΔG in kcal/mol). Component energies are null in FoldX 5.1 (simplified output format).
- ~40 seconds per protein average
- 300-second timeout per protein

**Step 7 — Collect results** (`07_foldx_collect.sh`):
Merges all per-protein JSON files into a single `foldx_stability_all.tsv`.

**Current status (2026-03-24):** 211/501 chunks completed (42%), 10,775/25,007 proteins processed. Two timeout chunks (125, 139) resubmitted with extended 6-hour wall time. Estimated completion: April 1-2, 2026.

#### Step 8: Module E — Unified Domain Architecture (`08_module_e_domains.sh`)

Integrates Phase 1 CATH assignments (1,390 proteins, curated) with Phase 2 Chainsaw predictions (23,868 proteins):
- CATH annotations take priority (higher confidence)
- Chainsaw fills gaps for proteins without CATH
- Computes domain count distribution per dataset (GroEL, HSP60, matrix, mito, proteome background)
- Loads Foldseek clusters for overlap analysis

**Result:** 25,258 proteins with domain assignments. Source column tracks whether each came from CATH or Chainsaw.

#### Step 9: Module F — N-vs-C Stability (`09_module_f_stability.sh`)

Full-scale version of Module F:
- Parses domain boundaries from both CATH (start/end columns) and Chainsaw (chopping notation)
- Defines three regions per protein
- Computes contact order using Plaxco definition (8 Å cutoff, min_seq_sep 6)
- Will integrate FoldX ΔG as primary stability metric once FoldX completes

**Result:** 2,648 paired N-vs-C comparisons (multi-domain proteins only). 5,322 region boundaries. 11,824 contact order records.

#### Step 10: Module H — Statistics (`10_module_h_stats.sh`)

Full-scale hypothesis testing with hierarchical BH correction:
- **56 tests** across 3 families (reduced from 281 in Phase 1 pilot by consolidating redundant tests)
- Compartment-matched AND size-matched backgrounds
- Effect sizes reported for every test

**Result:** 25/56 tests significant after hierarchical correction.

#### Step 11: Module I — Figures (`11_module_i.sh`)

Publication-quality figures with full-scale data:
- Real p-values (not placeholder)
- Real sample sizes
- Colorblind-friendly palette
- 300 DPI PDF + PNG

Locally polished using `module_i_polished.py` after transfer from HPC.

### FoldX Thermodynamic Stability

FoldX deserves special attention because it provides the thermodynamic stability dimension. **Completed April 1, 2026** (25,007 proteins, 0 failures). Key finding: GroEL substrates have significantly lower total energy (median -38.6 vs -15.2 for E. coli bg, p=2.9e-3, d=-0.07, compartment-matched vs background), while HSP60 substrates show no significant difference (p=0.77).

**What FoldX calculates:** The empirical free energy of folding (ΔG, kcal/mol) using a force field that includes:
- Van der Waals energy (attractive + repulsive)
- Solvation energy (hydrophobic burial + polar desolvation)
- Hydrogen bond energy (backbone + sidechain)
- Electrostatic energy (charge-charge interactions, ionic strength dependent)
- Entropy cost (sidechain + backbone conformational entropy lost upon folding)

A **more negative ΔG** = more thermodynamically stable fold.

**Why it matters:** Contact order tells us about folding *kinetics* (how fast/slow the protein folds). FoldX ΔG tells us about folding *thermodynamics* (how stable the folded state is). A protein with high contact order and strongly negative ΔG is one that folds slowly but ends up in a very stable state — potentially an ideal chaperonin substrate (needs help getting there, but once folded, stays folded).

**What changes after FoldX completes:**
1. Module F will be re-run with ΔG as an additional column in `n_vs_c_paired_full.tsv`
2. Module H will test: "Do N-domains have different ΔG than C-regions?" and "Is the ΔG difference substrate-specific?"
3. Module I will generate updated figures including FoldX violin plots
4. The complete analysis will then be manuscript-ready

---

## 9. Results and Interpretations

### Goal 1: Domain Architecture Results

#### Phase 1 Pilot (1,390 proteins)

**GroEL substrate enrichments (Fisher's exact test, size-matched E. coli cytoplasm background):**

| Superfamily | Name | Substrates | Background | OR | p (BH) | Interpretation |
|-------------|------|-----------|------------|-----|---------|---------------|
| 1.10.10.10 | Winged helix | 15 | ~1 expected | 42.8 | 2.35×10⁻⁶ | DNA-binding proteins with complex helix-turn-helix topology |
| 3.20.20.70 | Aldolase I / TIM barrel | 22 | ~5 expected | 9.16 | 2.35×10⁻⁶ | 8-stranded beta/alpha barrel — the iconic chaperonin substrate fold |

**HSP60 substrate enrichments (Fisher's exact, matrix background):**

| Superfamily | Name | Substrates | Background | OR | p (BH) | Interpretation |
|-------------|------|-----------|------------|-----|---------|---------------|
| 3.30.830.10 | - | High | Low | 5.4 | 1.98×10⁻³ | Matrix enzyme fold |
| 3.90.226.10 | - | High | Low | 4.8 | 2.78×10⁻³ | Mitochondrial fold |

**HSP60 depletions:**

| Superfamily | Name | OR | p (BH) | Interpretation |
|-------------|------|-----|---------|---------------|
| 1.50.40.10 | Mitochondrial carrier | 0.16 | 0.028 | Inner membrane carriers avoid HSP60 — they fold in the membrane, not the matrix |

#### Phase 2 Full-Scale (25,007 proteins)

The Phase 2 results confirmed and strengthened the pilot findings:

| Test | GroEL Phase 1 OR | GroEL Phase 2 OR | Interpretation |
|------|:---:|:---:|---|
| TIM barrel enrichment | 9.16 | 8.4 | Robust — confirmed at full scale |
| 1.10.10.10 enrichment | 42.8 | 50.9 | Even stronger with larger background |

**Key interpretation:** GroEL substrates are NOT just "difficult proteins." They are specifically enriched for certain topologically complex folds — particularly TIM barrels, which have a complex 8-fold symmetry that requires many sequential strand-helix-strand elements to come together correctly. The barrel cannot form until all 8 strands are present, making co-translational folding impossible — the protein MUST be fully synthesized before it can fold, creating a window of vulnerability where GroEL is essential.

**What about multi-domain enrichment?** Neither substrate set is enriched for having MORE domains (GroEL OR=1.13, p=0.35; HSP60 OR=0.85, p=0.16). Chaperonin substrate identity is determined by fold topology, not domain count.

### Goal 2: N-vs-C Structural Asymmetry Results

This is the project's most important and scientifically surprising finding.

#### The Universal N > C Contact Order Asymmetry

**Hypothesis H2.1 — CONFIRMED:** N-terminal domains have significantly higher relative contact order than C-terminal regions:

| Dataset | n (pairs) | Median N RCO | Median C RCO | Median diff | r (effect) | p-value |
|---------|:---------:|:---:|:---:|:---:|:---:|:---:|
| GroEL substrates | 124 | — | — | 0.043 | 0.41 | 8.9×10⁻⁵ |
| HSP60 substrates | 131 | — | — | 0.059 | 0.46 | 5.3×10⁻⁶ |
| Matrix background | 251 | — | — | 0.069 | 0.43 | 2.4×10⁻⁹ |
| **Mito background** | **425** | — | — | **0.064** | **0.48** | **7.1×10⁻¹⁸** |

The effect is medium-to-large (r = 0.41-0.48) and highly significant in every group.

#### The Crucial Negative Result: NOT Substrate-Specific

**Hypothesis H2.2 — REJECTED:** The N-vs-C asymmetry is NOT greater in chaperonin substrates than in background proteomes.

Mann-Whitney U test comparing the N-C contact order difference between substrates and background:
- GroEL substrates vs E. coli background: p = 0.058, r < 0.11 (NOT significant)
- HSP60 substrates vs mito background: p = 0.536, r < 0.11 (NOT significant)

**This means:** Background proteins show the SAME N-vs-C asymmetry as chaperonin substrates. The effect is universal to all multi-domain proteins.

#### No GroEL Class Effect

**Hypothesis H2.3 — REJECTED:** Class III obligate substrates do NOT show greater asymmetry than Class I spontaneous folders.

Kruskal-Wallis H-test: p = 0.77, η² = -0.014 (essentially zero effect size).

**Interpretation:** If N-vs-C asymmetry drove chaperonin dependence, Class III proteins (which absolutely require GroEL) should show the greatest asymmetry. They don't. The asymmetry is a general feature of protein architecture that has nothing to do with chaperonin biology specifically.

#### What Does This Mean Biologically?

The N > C contact order asymmetry reflects **co-translational folding constraints:**
1. Proteins are synthesized N-to-C on the ribosome
2. The N-terminal region emerges first and begins to fold while the rest is still being synthesized
3. Evolutionary pressure has placed the most topologically complex folds (high contact order) at the N-terminus, where they have the most time to fold during translation
4. C-terminal regions, emerging last, tend to adopt simpler folds or extend existing domains

This is a fundamental property of how proteins are organized, conserved across 2 billion years of evolution, and NOT specific to chaperonin substrates. It is a "gravitational constant" of protein architecture — always present, not caused by chaperonins.

### Goal 3: Mitochondrial Targeting Results

#### HSP60 Substrates Are Matrix-Enriched

**Hypothesis H3.1 — CONFIRMED:**

| Comparison | HSP60 substrates | Background | OR | 95% CI | p-value |
|-----------|:---:|:---:|:---:|:---:|:---:|
| Matrix vs non-matrix (within mito) | 46.6% matrix | 28.5% matrix | 3.29 | [2.46, 4.40] | 1.6×10⁻¹⁶ |

**Interpretation:** HSP60 substrates are 3.3× more likely to be matrix proteins than expected from the general mitochondrial proteome. This makes biological sense: HSP60 is a matrix chaperonin, so it primarily assists matrix-targeted proteins. However, 21.1% of HSP60 substrates enter the matrix by "non-canonical" pathways (no detectable transit peptide), suggesting alternative import mechanisms exist.

#### MTS Architecture: The "Landing Pad" Model

**Hypothesis H3.3 — CONFIRMED:**

Of all proteins with a transit peptide AND a CATH domain:
- **84.4%** have the MTS as a pre-domain extension (transit peptide cleaved BEFORE the first domain starts)
- **15.6%** have partial overlap between the MTS and the first domain
- Binomial test: p = 3.4×10⁻⁵¹ (overwhelmingly pre-domain)

**Gap statistics:** Median gap between transit peptide cleavage site and first domain start: **18 residues** (mean: 37.4, range: 0-579).

**The "landing pad" model:** After mitochondrial import and MTS cleavage, the protein emerges into the matrix with a short unstructured linker followed by its first structural domain in an unfolded state. This ~18-residue gap creates a "landing pad" that HSP60 can grip. The first domain, now free of the transit peptide and unfolded in the matrix environment, is optimally positioned for HSP60 to capture and assist its folding.

### Cross-Species Conservation Results

The 69 homolog pairs in Dataset 6 reveal remarkable conservation:

**N-domain contact order conservation:**
- Spearman correlation between GroEL substrate N-domain RCO and its HSP60 ortholog's N-domain RCO: **r = 0.84, p = 5.3×10⁻¹³** (n = 45 pairs with both values)

This means: across 2 billion years of evolution, the topological complexity of the N-terminal domain has been strongly conserved. A protein that has a complex N-terminal fold in *E. coli* has a correspondingly complex N-terminal fold in humans. The selective pressure maintaining N-terminal folding complexity is ancient and strong.

**Fold type conservation:**
- 55/69 homolog pairs (79.7%) share the same top CATH superfamily
- This confirms that chaperonin substrates are constrained to the same structural fold across species

### Statistical Summary

#### Phase 1 Pilot

| Family | Tests | Significant (BH) | Key Findings |
|--------|:-----:|:-----------------:|-------------|
| Domain architecture | ~200 | 5 | GroEL: TIM barrel (OR=9.2), Winged helix (OR=42.8) |
| Stability asymmetry | ~70 | 11 | N>C contact order universal; NOT substrate-specific |
| MTS targeting | ~11 | 6 | Matrix enrichment (OR=3.3); 84.4% pre-domain MTS |
| **Total** | **281** | **22** | |

#### Phase 2 Full-Scale

| Family | Tests | Significant (BH) | Key Findings |
|--------|:-----:|:-----------------:|-------------|
| Domain architecture | 24 | 9 | TIM barrel (OR=8.4), 1.10.10.10 (OR=50.9) |
| Stability asymmetry | 30 | 14 | N>C universal (r=0.41-0.48); substrate NS |
| MTS targeting | 2 | 2 | Matrix OR=3.29; pre-domain 84.4% |
| **Total** | **56** | **25** | |

**Why fewer tests in Phase 2?** The Phase 1 pilot tested many exploratory superfamilies. Phase 2 consolidated to pre-registered hypotheses only, reducing the multiple testing burden from 281 to 56 tests while increasing statistical power through larger sample sizes.

---

## 10. Biological Synthesis and Interpretation

The three goals converge into a coherent biological narrative:

### What Makes a Protein a Chaperonin Substrate?

**It's the fold, not the general difficulty.** Chaperonin substrates are enriched for specific topologies — particularly TIM barrels (GroEL OR=8.4) and winged helix domains (OR=50.9). These are folds with inherently complex topologies that cannot form co-translationally because they require all secondary structure elements to be present before the final fold can snap into place.

**It's NOT about N-vs-C asymmetry.** The N > C contact order asymmetry is real and strong (p = 10⁻¹⁸), but it is universal — present in ALL multi-domain proteins, not just chaperonin substrates. This asymmetry reflects the physics of vectorial translation (N-terminal synthesized first), not chaperonin biology.

**It's NOT about obligate dependence.** GroEL Class III (obligate) substrates show the same N-vs-C asymmetry as Class I (spontaneous) substrates (p=0.77 for class effect). Whatever makes a protein obligately dependent on GroEL is not captured by contact order differences between N and C.

### How Does the Mitochondrial System Work?

**HSP60 substrates are matrix residents.** 3.3× enriched for matrix localization vs the general mitochondrial proteome. This is expected — HSP60 is in the matrix.

**The MTS creates a "landing pad."** 84.4% of transit peptides are pre-domain extensions with a ~18-residue gap to the first domain. After import and MTS cleavage, the protein's first domain emerges unfolded into the matrix — exactly where HSP60 can engage it. This architecture is not coincidental; it facilitates chaperonin-assisted folding upon import.

### What is Conserved Across 2 Billion Years?

**The fold itself.** 79.7% of homolog pairs share the same CATH superfamily. TIM barrels that need GroEL in *E. coli* have TIM barrel orthologs that need HSP60 in humans.

**The N-domain complexity.** N-domain contact order is correlated across homologs (r=0.84). Whatever makes the N-terminal fold complex in *E. coli* is equally complex in the human ortholog.

### The Title Explained

"The End is the Beginning": The C-terminus (the "end" of translation) represents the "beginning" of the chaperonin's work — but NOT because C-terminal regions are more complex. Rather, the entire protein, having been synthesized N-to-C, emerges from the ribosome (or mitochondrial import channel) with its most complex domain already at the N-terminus. The chaperonin doesn't preferentially assist the C-terminal fold; it captures the entire protein and helps the globally complex topology snap into place. The asymmetry we observe is a footprint of translational history, conserved across all of life.

---

## 11. Limitations and Caveats

### Methodological Limitations

1. **pLDDT is confidence, not stability.** We use contact order for folding kinetics and FoldX for thermodynamics. FoldX (now complete) shows GroEL substrates are thermodynamically MORE stable (median -38.6 kcal/mol) — chaperonin assistance is needed for kinetic reasons. Caveat: FoldX was parameterized on experimental structures, not AlphaFold models.

2. **Co-IP captures interactions, not function.** HSP60 substrates are proteins that physically associate with HSP60 in a co-IP experiment. Some may be transient interactors, not true folding substrates. SILAC filtering mitigates this but cannot fully eliminate it.

3. **AlphaFold structures are predictions.** While generally excellent (mean pLDDT 85-93), some regions may be inaccurately predicted, especially for disordered regions or multi-state proteins.

4. **Contact order is an imperfect proxy.** RCO correlates with folding rate for two-state folders (r ≈ -0.75) but is less predictive for multi-state folders, knotted proteins, or proteins with complex co-translational folding pathways.

5. **MitoCarta annotations evolve.** 70 proteins changed compartment between v2.0 and v3.0. Future versions may further reclassify proteins.

### Statistical Limitations

6. **Multiple testing.** Despite hierarchical BH correction, 56 tests still carry risk of false positives. Effect sizes and biological plausibility should be considered alongside p-values.

7. **Power limitations for Dataset 6.** Only 69 homolog pairs (45 with paired contact order). Cross-species comparisons have limited statistical power for subtle effects.

8. **Size matching is approximate.** 10 kDa bins may not fully control for protein size effects, especially at the extremes of the size distribution.

### Biological Limitations

9. **No Class III equivalents for HSP60.** GroEL substrates have class annotations (I/II/III dependence). HSP60 substrates lack equivalent dependence classifications, preventing a direct comparison of "obligateness."

10. **MTS prediction gaps.** ~30-40% of matrix proteins lack a detectable transit peptide. These may use alternative import pathways (internal targeting signals, carrier-mediated import) that we cannot analyze without specialized predictors (TargetP 2.0, which requires a DTU license).

---

## 12. Software and Reproducibility

### Local Environment

| Component | Version |
|-----------|---------|
| Hardware | Apple M1, 8 GB RAM |
| OS | macOS Darwin 25.2.0 |
| Python | 3.9.16 (Anaconda) |
| pandas | 2.2.2 |
| scipy | 1.9.2 |
| matplotlib | (installed) |
| seaborn | (installed) |
| biopython | 1.78 |
| openpyxl | (installed) |
| Snakemake | 7.32.4 |
| mkdssp | 2.2.1 |

### Conda Environment `proteomics`

| Tool | Version | Architecture | Notes |
|------|---------|-------------|-------|
| MMseqs2 | 18.8cc5c | x86_64 (Rosetta 2) | Sequence search + orthology |
| Foldseek | 10.941cd33 | x86_64 (Rosetta 2) | Structural search + clustering |

### HPC Software

| Tool | Version | Path | Notes |
|------|---------|------|-------|
| FoldX | 5.1 | `/lustre/vishal.bharti/software/foldx5/foldx` | Academic license, expires 2026-12-31 |
| Chainsaw | latest (2026) | `/lustre/vishal.bharti/software/chainsaw` | ML domain prediction |
| STRIDE | heiniglab/stride | Compiled from source | NOT bioconda stride |

### Reproducibility

All analysis is fully reproducible:
1. **Phase 1**: Run scripts in `workflow/scripts/` in module order (A → I)
2. **Phase 2**: Deploy `workflow/phase2/` on HPC following `docs/HPC_DEPLOYMENT_GUIDE.md`
3. **Configuration**: All parameters in `workflow/phase2/config.yaml`
4. **Raw data**: All input data in `data/raw/` (downloadable from public databases)
5. **Processed data**: All intermediate files in `data/processed/` and `results/`

---

## 13. Session History

The project was developed across 6 interactive sessions:

| Session | Date | Key Accomplishments |
|:-------:|------|-------------------|
| 1 | 2026-03-10 | Project planning, dataset curation, GroEL standardization, HSP60 filtering, MitoCarta parsing |
| 2 | 2026-03-11 | FASTA extraction, RBH analysis, AlphaFold download, DSSP, CATH domains, Chainsaw, orthology, Dataset 6 |
| 3 | 2026-03-12 | Module F (N-vs-C), Module G (MTS), Module H (statistics), Module I (figures), Phase 1 complete |
| 4 | 2026-03-14 | Phase 2 pipeline development, HPC deployment, bug fixes (column names, STRIDE binary), Chainsaw full-scale |
| 5 | 2026-03-17 | Foldseek full-scale, Module E verified, F→H→I chain submitted, results transferred to Mac |
| 6 | 2026-03-22 | FoldX 5.1 installation, test run, full array job submitted (501 tasks), timeout fixes, collaborator deliverables |
| 7 | 2026-03-24 | PPT (35 slides) + Q&A guide created. FoldX ~42%. |
| 8 | 2026-03-25 | GitHub repo + reproducibility infra. FoldX ~48%. |
| 9 | 2026-04-01 | FoldX 100% COMPLETE (25,007, 0 failures). F→H→I chain re-run. 60 tests, 28 sig. Column bugs fixed. Manuscript next. |

---

## 14. Appendix: Complete File Inventory

### Input Data Files

| File | Path | Size | Description |
|------|------|------|-------------|
| E. coli proteome FASTA | `data/raw/uniprot/ecoli_k12_proteome.fasta` | ~2.5 MB | 4,403 proteins |
| E. coli proteome TSV | `data/raw/uniprot/ecoli_k12_proteome.tsv` | ~1 MB | Metadata |
| Human proteome FASTA | `data/raw/uniprot/human_proteome.fasta` | ~12 MB | 20,416 proteins |
| Human proteome TSV | `data/raw/uniprot/human_proteome.tsv` | ~5 MB | Metadata |
| MitoCarta 3.0 | `data/raw/mitocarta/Human.MitoCarta3.0.xls` | ~5 MB | Mitochondrial proteome |
| AlphaFold structures (pilot) | `data/raw/alphafold/pilot/AF-*.cif` | ~466 MB | 1,382 CIF files |
| GroEL substrates (standardized) | `data/processed/groel_substrates_standardized.tsv` | ~60 KB | 252 proteins, 16 columns |
| GroEL FASTA | `data/processed/groel_substrates.fasta` | ~100 KB | 252 sequences |
| HSP60 Tier-1 substrates | `data/processed/hsp60_tier1_substrates.tsv` | ~80 KB | 266 proteins, 28 columns |
| HSP60 Tier-1 FASTA | `data/processed/hsp60_tier1_substrates.fasta` | ~120 KB | 266 sequences |
| HSP60 full standardized | `data/processed/hsp60_interactome_standardized.tsv` | ~100 KB | All tiers |
| Mitochondrial proteome | `data/processed/human_mito_proteome.tsv` | ~120 KB | 1,136 proteins |
| Matrix proteome | `data/processed/human_matrix_proteome.tsv` | ~55 KB | 525 proteins |
| Homolog pairs | `data/processed/groel_hsp60_homologs.tsv` | ~10 KB | 69 pairs |

### Phase 1 Result Files

| File | Path | Records | Description |
|------|------|:-------:|-------------|
| RBH annotated | `results/homology/rbh_groel_hsp60_annotated.tsv` | 40 | Reciprocal best hits |
| Orthogroups | `results/homology/orthogroups_ecoli_human.tsv` | 422 | Full orthogroups |
| Substrate orthogroups | `results/homology/substrate_orthogroups.tsv` | 34 | Shared substrate orthogroups |
| Structure index | `results/structures/structure_index.tsv` | 1,382 | AlphaFold metadata |
| DSSP summary | `results/structures/dssp_summary.tsv` | 1,382 | SS composition |
| DSSP per-residue | `results/structures/dssp_per_residue.tsv` | ~500K | Per-residue SS |
| Quality validation | `results/structures/structure_quality_validation.tsv` | 1,382 | Quality tiers |
| Flagged structures | `results/structures/flagged_low_quality_structures.tsv` | 63 | Low-quality flags |
| CATH domains | `results/domains/cath_domain_assignments.tsv` | 2,141 | Gene3D domains |
| CATH protein summary | `results/domains/cath_protein_summary.tsv` | 1,390 | Per-protein summary |
| Chainsaw predictions | `results/domains/chainsaw_domain_predictions.tsv` | ~500 | ML predictions |
| Unified domains | `results/domains/ml_domain_assignments.tsv` | 1,387 | CATH + Chainsaw |
| Domain metrics | `results/domains/domain_structural_metrics.tsv` | 2,141 | Per-domain SS + pLDDT |
| Foldseek clusters | `results/domains/foldseek_clusters.tsv` | 1,382 | Structural clusters |
| Domain distribution | `results/domains/domain_distribution_summary.tsv` | ~200 | Summary statistics |
| Region boundaries | `results/termini/region_boundaries.tsv` | ~3,000 | Three-region definition |
| N-vs-C paired | `results/termini/n_vs_c_paired.tsv` | ~1,200 | Paired comparisons |
| Contact order | `results/termini/contact_order.tsv` | ~4,000 | RCO per region |
| Combined targeting | `results/mts/combined_targeting.tsv` | ~1,100 | Targeting classification |
| MTS-domain relationship | `results/mts/mts_domain_relationship.tsv` | ~400 | MTS gap analysis |
| Corrected p-values | `results/stats/corrected_pvalues.tsv` | 281 | All Phase 1 tests |
| Domain enrichment | `results/stats/domain_enrichment.tsv` | ~200 | Superfamily enrichment |
| Stability comparisons | `results/stats/stability_comparisons.tsv` | ~70 | N-vs-C test details |
| Targeting stats | `results/stats/targeting_stats.tsv` | ~10 | MTS statistics |
| Figures | `results/figures/fig[1-6]_*.{pdf,png}` | 12 | Phase 1 figures |

### Phase 2 Result Files

| File | Path | Records | Description |
|------|------|:-------:|-------------|
| Structure index (full) | `results/phase2/structures/structure_index_full.tsv` | 25,007 | All proteins |
| Chainsaw (full) | `results/phase2/domains/chainsaw_full_predictions.tsv` | 23,868 | ML predictions |
| Unified domains (full) | `results/phase2/domains/unified_domain_assignments_full.tsv` | 25,258 | All assignments |
| Domain distribution (full) | `results/phase2/domains/domain_distribution_full.tsv` | ~500 | Distribution |
| Region boundaries (full) | `results/phase2/stability/region_boundaries_full.tsv` | 5,322 | Three regions |
| N-vs-C paired (full) | `results/phase2/stability/n_vs_c_paired_full.tsv` | 2,648 | Paired comparisons |
| Contact order (full) | `results/phase2/stability/contact_order_full.tsv` | 11,824 | Per-region RCO |
| Sequence metrics (full) | `results/phase2/stability/sequence_metrics_full.tsv` | 11,824 | Charge, hydrophobicity |
| Structure metrics (full) | `results/phase2/stability/structure_metrics_full.tsv` | 11,824 | pLDDT, SS fractions |
| Corrected p-values (full) | `results/phase2/stats/corrected_pvalues_full.tsv` | 56 | All Phase 2 tests |
| Statistics summary | `results/phase2/stats/statistics_summary_full.txt` | — | Human-readable |
| Foldseek clusters (full) | `results/phase2/foldseek/analysis/foldseek_clusters_full.tsv` | 16,193 | Structure clusters |
| Figures (full) | `results/phase2/figures/fig[1-6]_*.{pdf,png}` | 12 | Publication figures |
| FoldX results | `results/phase2/foldx/foldx_stability_all.tsv` | 25,007 | COMPLETE |

### Scripts Inventory

| Script | Module | Purpose |
|--------|--------|---------|
| `scripts/validate_uniprot_accessions.py` | A | GroEL accession demerging and standardization |
| `scripts/filter_hsp60_interactome.py` | A | HSP60 SILAC filtering and tier assignment |
| `scripts/module_c_extract_fasta.py` | B | FASTA sequence extraction from proteomes |
| `scripts/module_c_analyze_rbh.py` | C | RBH analysis and annotation |
| `workflow/scripts/parse_mitocarta.py` | A | MitoCarta 3.0 parsing |
| `workflow/scripts/download_alphafold_pilot.py` | D | AlphaFold CIF download (pilot) |
| `workflow/scripts/run_dssp.py` | D | DSSP secondary structure assignment |
| `workflow/scripts/validate_structure_quality.py` | D | Quality tier assignment |
| `workflow/scripts/get_cath_domains.py` | E | CATH/Gene3D domain query |
| `workflow/scripts/run_chainsaw_e2.py` | E | Chainsaw ML domain prediction |
| `workflow/scripts/compute_domain_structural_metrics.py` | E | Per-domain SS + pLDDT |
| `workflow/scripts/analyze_foldseek.py` | E | Foldseek cluster analysis |
| `workflow/scripts/domain_distribution_summary.py` | E | Domain architecture summary |
| `workflow/scripts/run_orthology.py` | C | OrthoFinder-style orthology |
| `workflow/scripts/build_dataset6_homologs.py` | C | Dataset 6 construction |
| `workflow/scripts/module_f_n_vs_c_analysis.py` | F | N-vs-C stability (CATH proteins) |
| `workflow/scripts/module_f_extension_chainsaw.py` | F | N-vs-C extension (Chainsaw proteins) |
| `workflow/scripts/module_g_mts_analysis.py` | G | MTS analysis |
| `workflow/scripts/module_h_comparative_stats.py` | H | Statistical testing |
| `workflow/scripts/generate_figures.py` | I | Publication figures |
| `workflow/phase2/download_alphafold_full.py` | D (full) | Bulk AlphaFold download |
| `workflow/phase2/run_foldseek_full.py` | E (full) | Full-scale Foldseek |
| `workflow/phase2/run_foldx.py` | F (full) | FoldX stability |
| `workflow/phase2/scripts/module_f_full.py` | F (full) | Full-scale N-vs-C |
| `workflow/phase2/scripts/module_h_full.py` | H (full) | Full-scale statistics |
| `workflow/phase2/scripts/module_i_full.py` | I (full) | Full-scale figures |
| `workflow/phase2/scripts/module_i_polished.py` | I (full) | Polished figures |

### SLURM Job Scripts

| Script | Resources | Purpose |
|--------|-----------|---------|
| `00_setup.sh` | 2 GB, 1 CPU | Environment validation |
| `01_download_alphafold.sh` | 4 GB, 1 CPU | Download 25,007 structures |
| `02_foldseek_createdb.sh` | 16 GB, 4 CPUs | Build structure databases |
| `03_foldseek_search.sh` | 64 GB, 16 CPUs | All-vs-all structural search |
| `04_foldseek_cluster.sh` | 64 GB, 16 CPUs | Cluster similar structures |
| `05_chainsaw.sh` | 16 GB, 4 CPUs | ML domain prediction (72h) |
| `06_foldx_generate.sh` | 2 GB, 1 CPU | Generate FoldX array script |
| `07_foldx_collect.sh` | 4 GB, 1 CPU | Merge FoldX results |
| `08_module_e_domains.sh` | 16 GB, 4 CPUs | Unified domain assignments |
| `09_module_f_stability.sh` | 16 GB, 4 CPUs | N-vs-C stability analysis |
| `10_module_h_stats.sh` | 16 GB, 4 CPUs | Statistical testing |
| `11_module_i.sh` | 8 GB, 2 CPUs | Publication figures |
| `submit_pipeline.sh` | — | Master submission script |
| `submit_analysis.sh` | — | Analysis chain (F→H→I) |

---

*Document generated: March 24, 2026*
*Project: Antah Asti Prarambh — Comparative Structural Proteomics of Chaperonin Substrates*
*Author: Vishal Bharti, CSIR-IGIB*
