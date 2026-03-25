# Antah Asti Prarambh: Pre-registered Primary Hypotheses

Document created: 2026-03-12
Module: H (Comparative Statistics)

## Goal 1: Domain Architecture

**H1.1** GroEL substrates (especially Class III) are enriched for specific CATH
superfamilies (TIM barrels, Rossmann folds) compared to size-matched cytoplasmic
E. coli background.

**H1.2** HSP60 substrates show enrichment for specific fold architectures compared
to size-matched mitochondrial matrix background.

**H1.3** Structural fold enrichment is conserved between GroEL and HSP60 systems
(shared orthogroup substrates show similar fold preferences).

## Goal 2: N-vs-C Stability Asymmetry

**H2.1** N-terminal domains of chaperonin substrates have different stability
proxies (contact order) compared to C-terminal regions (within-protein paired
test).

**H2.2** The N-vs-C asymmetry is greater in chaperonin substrates than in matched
background proteins.

**H2.3** Class III (obligate) GroEL substrates show greater N-vs-C asymmetry than
Class I (spontaneous folders).

## Goal 3: Matrix Targeting and MTS

**H3.1** HSP60 substrates are enriched for matrix localization compared to the
general mitochondrial proteome.

**H3.2** MTS-bearing HSP60 substrates have distinct first-domain properties
compared to non-MTS-bearing substrates.

**H3.3** The MTS is predominantly a pre-domain extension, not part of the first
structural domain.

## Statistical Framework

- **Multiple testing**: Hierarchical BH correction at two levels:
  - Level 1: Three goal families (domain, stability, targeting)
  - Level 2: Within each family, BH correction on individual tests
- **Primary test level**: alpha = 0.05 (two-sided)
- **Effect size reporting**: Odds ratios with 95% CI (Fisher tests),
  rank-biserial correlation (Wilcoxon), eta-squared (Kruskal-Wallis)
- **Controls**: Compartment-matched AND size-matched (10 kDa bins)
