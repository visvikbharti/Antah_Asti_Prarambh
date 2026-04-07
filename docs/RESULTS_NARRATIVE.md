# Antah Asti Prarambh: Results Narrative

*A comparative structural analysis of Group I chaperonin substrates across the prokaryote-eukaryote divide*

---

## 1. Dataset Assembly and Quality Control

The study assembled six primary datasets spanning two Group I chaperonin systems and their respective compartment-matched backgrounds. For the bacterial system, 252 *E. coli* GroEL substrates were drawn from the Kerner et al. (2005) co-immunoprecipitation interactome, classified into obligate dependence tiers: 38 Class I (spontaneous folders), 126 Class II (partial dependence), 84 Class III (obligate substrates), and 4 ambiguous Class I/II proteins. All 252 UniProt accessions were validated against the current UniProt database via REST API; 149 accessions required remapping due to demerged strain-specific entries, and all were successfully resolved to K-12 reference proteome entries with zero obsolete accessions. The compartment-matched background comprised the full *E. coli* cytoplasmic proteome (663 proteins after size-matching in 10 kDa bins).

For the eukaryotic system, the HSP60 interactome was derived from Bruderer et al. SILAC co-immunoprecipitation data. The raw dataset contained 325 proteins, which were filtered through a multi-criterion quality pipeline. Proteins were classified as bait (n = 2), co-chaperone (n = 4), likely contaminant (n = 4), or candidate substrate (n = 315). Quality tiering based on SILAC enrichment ratios yielded 266 Tier 1 substrates (median SILAC ratio = 22.2, range 5.1--91.0), 49 Tier 2 (median SILAC ratio = 8.7), and 10 excluded proteins. All 266 Tier 1 substrates were annotated in MitoCarta (100.0%), confirming their mitochondrial provenance.

A critical quality control step involved reconciling MitoCarta 2.0 annotations (used in the original HSP60 study) with MitoCarta 3.0 (Rath et al., 2021). Of the 325 HSP60 interactors, 278 (85.5%) were annotated in MitoCarta 2.0 but only 274 (84.3%) in MitoCarta 3.0; 6 proteins were lost (including GAPDH, P4HB, and two ribosomal proteins) while 2 were gained (CDK5RAP1, NSUN2). More consequentially, 70 proteins underwent sub-mitochondrial localization reclassification between versions. Of these, 51 were reclassified from matrix to inner membrane (MIM) -- primarily respiratory chain subunits (SDHA, NDUFS1-4, NDUFV1-2, UQCRC1-2) and ATP synthase components (ATP5A1, ATP5B, ATP5D, ATP5F1, ATP5H, ATP5L, ATP5O). Conversely, 18 proteins were newly assigned to matrix in MitoCarta 3.0 (including PC, PCCB, PCK2, SUCLG1-2, PDHB). This reclassification reduced the count of HSP60 interactors annotated as matrix-localized from 227 (MitoCarta 2.0) to 194 (MitoCarta 3.0), a substantial shift that underscores the importance of using updated annotations.

Two background datasets were assembled for the human system: the full mitochondrial proteome from MitoCarta 3.0 (1,136 proteins, of which 525 are matrix-localized) and the matrix-only subset (524 proteins after filtering). These served as compartment-matched controls for the HSP60 enrichment analyses.

## 2. Cross-Species Orthology

### 2.1 Reciprocal Best Hit Analysis

To identify conserved chaperonin substrates across the prokaryote-eukaryote divide, we performed reciprocal best hit (RBH) analysis using MMseqs2 (v18.8, default parameters, e-value < 0.001). Of the 252 GroEL substrates and 266 HSP60 Tier 1 substrates, 40 RBH pairs were identified, representing 15.9% of GroEL and 15.0% of HSP60 substrates (Fig 5).

The 40 RBH pairs showed strong sequence conservation: median percent identity 35.8% (range 21.6--60.1%), median e-value 3.39 x 10^-58, and median query coverage 0.924. The e-value distribution was heavily skewed toward high-confidence matches, with 15 pairs (37.5%) at e-values below 10^-100 and only 4 pairs (10.0%) in the moderate range (10^-20 to 10^-10).

The distribution of RBH matches across GroEL dependence classes was informative: Class I contributed 7/38 substrates (18.4%), Class II contributed 24/126 (19.0%), and Class III contributed 8/84 (9.5%). The lower fraction of Class III matches may reflect the specialized nature of obligate substrates, many of which are taxonomically restricted to bacteria. The 8 Class III orthologs included metabolically essential enzymes such as SDHA (succinate dehydrogenase, 55.4% identity), ALDH2 (aldehyde dehydrogenase, 42.6% identity), and LAP3 (leucyl aminopeptidase, 35.1% identity). The top-scoring RBH pair overall was sucA/OGDH (2-oxoglutarate dehydrogenase, 652 bits, e-value 6.17 x 10^-204).

### 2.2 Orthogroup Analysis

To capture many-to-many orthology relationships arising from gene duplication, we performed MMseqs2 all-vs-all reciprocal searches followed by connected-component clustering (e-value < 10^-5, minimum 25% identity, minimum 50% coverage). This identified 422 orthogroups with a median size of 3 proteins (range 2--101). Of these, 199 were 1-to-1 orthogroups, 114 were 1-to-many (*E. coli* to *H. sapiens*), 40 were many-to-1, and 69 were many-to-many.

Among the 422 orthogroups, 143 contained at least one chaperonin substrate. Of these, 34 were shared (containing both GroEL and HSP60 substrates), 51 contained GroEL substrates only, and 58 contained HSP60 substrates only. The 34 shared orthogroups encompassed 62 substrate pairs, a substantial increase over the 40 RBH pairs. The 29 additional pairs arose from paralogous expansions; for example, GroEL substrate aldB (Class III) mapped to an orthogroup containing five human aldehyde dehydrogenases (ALDH2, ALDH1B1, ALDH9A1, ALDH7A1, ALDH5A1), of which four are HSP60 substrates. Similarly, GroEL substrate fadA mapped to an orthogroup containing three human thiolases (ACAT1, ACAA2, HADHB).

### 2.3 The Merged 69-Pair Dataset (Dataset 6)

Merging RBH and orthogroup evidence produced Dataset 6, comprising 69 unique cross-species substrate pairs. Of the original 40 RBH pairs, 33 were confirmed by orthogroup membership; the remaining 7 RBH pairs were not captured by orthogroup clustering (likely due to the more stringent coverage thresholds). The orthogroup approach contributed 29 additional pairs, yielding 69 total. This dataset served as the basis for cross-species conservation analyses in Goals 1 and 2.

Among these 69 homolog pairs, 55 (79.7%) shared the same top CATH superfamily assignment (Jaccard index for shared superfamilies = 0.202; hypergeometric p = 1.0), indicating that fold conservation tracks closely with sequence orthology for chaperonin substrates.

### 2.4 Biological Interpretation

The conserved chaperonin clients are dominated by core metabolic enzymes: TCA cycle components (OGDH/sucA, SDHA/sdhA, DLST/sucB), amino acid biosynthetic enzymes (SHMT2/glyA, OAT/astC), and mitochondrial translation machinery (GFM1/fusA, TARS2/thrS, NARS2/asnS). The 51 GroEL-only orthogroups likely represent cases of lost chaperonin dependence in eukaryotes or incomplete coverage of the HSP60 interactome, while the 58 HSP60-only orthogroups may reflect gained dependence in the mitochondrial context or undetected GroEL substrates.

## 3. Goal 1: Structural Domain Architecture

### 3.1 Domain Coverage and Architecture

CATH domain assignments were obtained for 18,855 proteins (75.3% of the 25,007 total) using the InterPro API, yielding 51,667 structural domains. Chainsaw ML predictions provided domain boundaries for the remaining 6,164 proteins (24.7%), giving 100% coverage (25,019 unified records). Assignment rates varied by dataset: GroEL substrates had the highest coverage at 247/252 (98.0%), followed by HSP60 at 241/266 (90.6%), matrix background at 471/524 (89.9%), and full mitochondrial proteome at 898/1,132 (79.3%).

Mean domain counts were remarkably consistent across datasets: 2.02 domains per protein for GroEL, 1.95 for HSP60, 2.01 for matrix background, and 1.81 for the full mitochondrial proteome. The domain count distributions are presented in Fig 1.

### 3.2 CATH Class Distribution

All datasets were dominated by alpha-beta domains, consistent with the known preference of Group I chaperonins for this structural class (Fig 1). GroEL substrates comprised 66.3% alpha-beta, 19.8% mainly-alpha, 11.6% mainly-beta, 1.2% few-secondary-structure, and 1.0% class 6 domains. HSP60 substrates showed a similar profile: 67.5% alpha-beta, 21.9% mainly-alpha, 8.5% mainly-beta.

Formal chi-squared tests confirmed that the CATH class distribution of chaperonin substrates differed significantly from their respective backgrounds. For GroEL versus the *E. coli* cytoplasmic background: chi-squared p = 5.23 x 10^-21, Cramer's V = 0.120. For HSP60 versus the mitochondrial background: chi-squared p = 2.39 x 10^-24, Cramer's V = 0.101. Both are highly significant, with modest effect sizes (Cramer's V approximately 0.1) indicating that alpha-beta dominance is a general cellular trend that chaperonin substrates accentuate modestly.

A notable difference between datasets was the representation of mainly-alpha domains, which constituted 28.2% of the full mitochondrial proteome but only 19.8% of GroEL substrates -- consistent with the observation that helical membrane proteins are abundant in the mitochondrial proteome but are generally not chaperonin clients.

### 3.3 Superfamily Enrichment

#### GroEL Substrates

Fisher's exact tests with Benjamini-Hochberg correction across 123 superfamilies identified 5 significantly enriched or depleted superfamilies in GroEL substrates relative to the size-matched *E. coli* background (Fig 1):

1. **Winged helix-like DNA-binding domain (1.10.10.10)**: OR = 42.80 [95% CI: 5.62, 325.83], p_BH = 2.35 x 10^-6. This was the most dramatically enriched superfamily, with 15 GroEL substrates versus only 1 background protein carrying this fold. These winged-helix domains are commonly found as small regulatory modules appended to larger catalytic domains.

2. **Aldolase class I / TIM barrel (3.20.20.70)**: OR = 22.6, p_BH = 2.35 x 10^-6. The canonical TIM barrel is the archetype of a GroEL-dependent fold. This confirms decades of biochemical observation that TIM barrels are preferential chaperonin substrates.

3. **Mitochondrial carrier domain (1.50.40.10)**: OR = 0.00 [0.00, 0.42], p_BH = 1.20 x 10^-5. This superfamily was completely absent from GroEL substrates (0/247) but present in 48 background proteins, reflecting the fact that mitochondrial carriers are eukaryotic innovations not found in the *E. coli* cytoplasm.

4. **Superfamily 3.40.50.2300**: OR = infinity [2.00, 636.43], p_BH = 1.18 x 10^-2. Six GroEL substrates versus zero background proteins.

5. **Metal-dependent hydrolases (3.20.20.140)**: OR = infinity [1.66, 546.35], p_BH = 3.52 x 10^-2. Five GroEL substrates versus zero background proteins.

#### HSP60 Substrates

Among 119 superfamilies tested for HSP60, only one reached significance:

1. **Mitochondrial carrier domain (1.50.40.10)**: OR = 0.16 [0.05, 0.52], p_BH = 2.79 x 10^-2. HSP60 substrates contained only 3 proteins with this domain versus 48 in the mitochondrial background. This depletion is biologically coherent: mitochondrial carriers are integral inner membrane proteins that follow the carrier import pathway (TIM22 complex) rather than matrix import via TIM23/HSP60.

The asymmetry between GroEL (5 significant superfamilies) and HSP60 (1 significant) likely reflects both differences in statistical power (larger background set for GroEL) and the greater structural heterogeneity of the mitochondrial matrix proteome.

#### Cross-Species Conservation of Fold Preferences

Among the 69 homolog pairs in Dataset 6, 55 (79.7%) shared the same top CATH superfamily assignment. The two systems shared 85 CATH superfamilies (GroEL-only: 179; HSP60-only: 156; Jaccard index = 0.202). A hypergeometric test for excess sharing was not significant (p = 1.0), indicating that the overlap in superfamily repertoire is consistent with the general overlap between bacterial and mitochondrial proteomes rather than specific convergence on chaperonin-preferred folds.

### 3.4 Structural Clustering (Foldseek)

Foldseek structural clustering (cascaded set-cover algorithm; minimum 30% sequence identity, 50% coverage, e-value < 0.01) across all 25,007 proteins yielded 16,242 clusters from 27,063 entries (Fig 1). The distribution was heavily skewed toward singletons.

Per-dataset cluster counts were: GroEL 239 clusters, HSP60 240 clusters.

Twenty-five clusters contained members from both the GroEL and HSP60 substrate sets. The largest shared clusters included Cluster 2 (aldehyde dehydrogenases: GroEL aldB with HSP60 ALDH2, ALDH1B1, ALDH9A1, ALDH5A1; 7 members total), Cluster 3 (short-chain dehydrogenases: GroEL P0AEK2 with HSP60 Q7Z4W1, Q8N4T8, Q92506, Q99714; 7 members), and Cluster 8 (thiolases: GroEL fadA with HSP60 ACAT1, ACAA2, HADHB; 5 members). These shared clusters recapitulate the orthology-based findings, confirming that structural similarity underlies the conservation of chaperonin dependence.

## 4. Goal 2: N-terminal Domain vs C-terminal Region

### 4.1 Region Definition

For multi-domain proteins, the first CATH domain (by sequence position) was designated the "N-terminal domain" and all sequence C-terminal to the last residue of the first domain constituted the "C-terminal region." For single-domain proteins, the sole domain served as the N-terminal domain and any C-terminal tail beyond the domain boundary served as the C-terminal region. The paired N-vs-C analysis was restricted to proteins with both regions of sufficient length. The resulting paired datasets comprised: GroEL n = 121 (for contact order) to n = 144 (for other metrics), HSP60 n = 121 to n = 135, matrix background n = 144 to n = 168, and full mitochondrial background n = 110 to n = 118.

### 4.2 Sequence-Derived Properties

Columns in the paired dataset included length, net charge, and fraction hydrophobic residues for each region. These were computed for each protein but were not the primary focus of statistical testing. Compositional differences between N-terminal domains and C-terminal regions were generally modest, with no consistent substrate-specific patterns emerging beyond what is expected from the known N-terminal enrichment of structured content.

### 4.3 Structure-Derived Properties

#### pLDDT Confidence

AlphaFold predicted Local Distance Difference Test (pLDDT) scores were computed for each region as a proxy for structural confidence and order. For GroEL substrates, N-terminal domains had significantly higher mean pLDDT than C-terminal regions (median N = 95.0, median C = 93.3; Wilcoxon signed-rank W = 3,750, p_BH = 6.13 x 10^-3, r = 0.282, direction N > C). The fraction of residues with pLDDT > 70 was also significantly higher in N-terminal domains (W = 1,623, p_BH = 4.16 x 10^-4, r = 0.438).

For HSP60 substrates, mean pLDDT showed a trend but did not reach significance after correction (W = 3,752, p_BH = 0.101, r = 0.183). The fraction with pLDDT > 70 was marginal (W = 2,081.5, p_BH = 0.059, r = 0.238, direction C > N).

#### Secondary Structure

Neither helical content (frac_helix) nor strand content (frac_strand) showed significant N-vs-C differences in any dataset after correction. For GroEL: frac_helix W = 5,001, p_BH = 0.662, r = 0.042; frac_strand W = 4,471, p_BH = 0.239, r = 0.132. For HSP60: frac_helix W = 4,125, p_BH = 0.444, r = 0.088; frac_strand W = 3,857, p_BH = 0.239, r = 0.134.

### 4.4 Contact Order -- The Primary Stability Proxy

Relative contact order (RCO) served as the principal structural proxy for folding complexity. Higher contact order indicates more long-range contacts, which is associated with slower folding rates and greater dependence on chaperone assistance. Within-protein paired Wilcoxon signed-rank tests revealed a consistent and significant N > C asymmetry across all four datasets (Fig 2):

- **GroEL substrates** (n = 121): W = 2,262, p_BH = 8.20 x 10^-4, r = 0.387 (medium effect), direction N > C. Median N-domain RCO = 0.265; median C-region RCO = 0.234; median difference = 0.041.

- **HSP60 substrates** (n = 121): W = 1,776, p_BH = 7.34 x 10^-6, r = 0.519 (large effect), direction N > C. Median N-domain RCO = 0.266; median C-region RCO = 0.200; median difference = 0.065.

- **Matrix background** (n = 144): W = 3,195, p_BH = 3.59 x 10^-4, r = 0.388 (medium effect), direction N > C. Median N-domain RCO = 0.251; median C-region RCO = 0.193; median difference = 0.080.

- **Full mitochondrial background** (n = 110): W = 1,086, p_BH = 9.00 x 10^-8, r = 0.644 (large effect), direction N > C. Median N-domain RCO = 0.255; median C-region RCO = 0.201; median difference = 0.059.

Multivariate confirmation using Hotelling's T-squared tests (combining contact order, pLDDT, helical fraction, and strand fraction) was significant for all four datasets: GroEL T^2 = 18.08, F(4,117) = 4.41, p = 2.36 x 10^-3; HSP60 T^2 = 27.92, F(4,117) = 6.80, p = 5.84 x 10^-5; matrix background T^2 = 31.10, F(4,140) = 7.61, p = 1.41 x 10^-5; mitochondrial background T^2 = 60.19, F(4,106) = 14.63, p = 1.51 x 10^-9.

### 4.5 Critical Finding: Universal Asymmetry

The most important result of Goal 2 is a negative one for the original hypothesis. The N > C contact order asymmetry is NOT specific to chaperonin substrates; it is a universal feature of multi-domain protein architecture (Fig 2).

Mann-Whitney U tests comparing the magnitude of N-vs-C contact order asymmetry between substrate and background datasets revealed no significant differences in any comparison:

- GroEL vs. mitochondrial background: U = 5,921, p_BH = 0.570, r = 0.110
- HSP60 vs. matrix background: U = 8,544, p_BH = 0.949, r = 0.019
- HSP60 vs. mitochondrial background: U = 6,357, p_BH = 0.949, r = 0.045

Similarly, no substrate-specific differences were found for pLDDT asymmetry (GroEL vs. mito_bg: U = 7,696, p_BH = 0.570, r = 0.094; HSP60 vs. matrix_bg: U = 10,114, p_BH = 0.570, r = 0.108), helical asymmetry, or strand asymmetry (all p_BH > 0.57).

This result has significant implications. The observation that N-terminal domains universally have higher contact order than C-terminal regions suggests that this asymmetry is a fundamental feature of protein architecture -- likely reflecting the co-translational folding constraint whereby N-terminal domains must fold first and do so under the kinetic pressure of ribosomal emergence. Chaperonin substrates are embedded within this general architectural trend; they do not exhibit exaggerated asymmetry despite their known folding difficulties. The original hypothesis (H2.2) that chaperonin substrates would show greater N-vs-C asymmetry is therefore rejected.

### 4.6 GroEL Class Comparison

Kruskal-Wallis tests comparing the N-vs-C asymmetry across GroEL dependence classes (Class I, n = 19--25; Class II, n = 64--70; Class III, n = 35--46) revealed no significant differences for any metric (Fig 3):

- Contact order asymmetry: H = 0.32, p_BH = 0.964, eta-squared = -0.014
- pLDDT asymmetry: H = 0.95, p_BH = 0.964, eta-squared = -0.008
- Helical fraction asymmetry: H = 0.07, p_BH = 0.964, eta-squared = -0.014
- Strand fraction asymmetry: H = 0.18, p_BH = 0.964, eta-squared = -0.013

The negative eta-squared values indicate effect sizes indistinguishable from zero. This rejects hypothesis H2.3: obligate GroEL substrates (Class III) do not show greater N-vs-C structural asymmetry than spontaneous folders (Class I). Whatever determines the degree of chaperonin dependence, it is not captured by the N-vs-C contact order differential.

### 4.7 Cross-Species Conservation

Among the 69 homolog pairs in Dataset 6, we assessed whether the relative contact order of N-terminal domains was conserved between the GroEL and HSP60 orthologs. The Spearman correlation was r = 0.82 (Fig 5), indicating strong cross-species conservation of first-domain folding complexity. Homolog pairs that are jointly chaperonin substrates in both *E. coli* and humans tend to have similarly complex N-terminal folds, consistent with the conservation of the structural features that necessitate chaperonin assistance.

## 5. Goal 3: Mitochondrial Matrix Targeting

### 5.1 Targeting Classification

Module G analyzed the targeting signals and sub-mitochondrial localization of all 1,139 proteins (266 HSP60 Tier 1 substrates, 1,133 MitoCarta entries, with 260 overlapping). The overall targeting breakdown was: inner membrane 391 (34.3%), high-confidence matrix 325 (28.5%), non-canonical matrix import 192 (16.9%), outer membrane 110 (9.7%), intermembrane space 52 (4.6%), and other/unspecified 50 (4.4%) (Fig 4).

Among the 266 HSP60 Tier 1 substrates specifically, the distribution was markedly shifted toward matrix localization: high-confidence matrix 124 (46.6%), inner membrane 71 (26.7%), non-canonical matrix import 56 (21.1%), with smaller fractions in IMS (3, 1.1%), outer membrane (2, 0.8%), and other categories. A total of 181 HSP60 substrates (68.0%) were classified as matrix-targeted (combining high-confidence matrix, non-canonical matrix, and probable matrix categories). Transit peptides were detected in 168 of 266 HSP60 substrates (63.2%), compared to 494/1,139 (43.4%) overall.

### 5.2 Matrix Enrichment

HSP60 substrates were significantly enriched for matrix localization compared to the general mitochondrial proteome (Fisher's exact test; Fig 4): 181 HSP60 substrates were matrix-localized versus 85 non-matrix, compared to 343 non-HSP60 matrix versus 530 non-HSP60 non-matrix proteins. This yielded OR = 3.29 [95% CI: 2.46, 4.40], p = 1.60 x 10^-16. This strong enrichment confirms hypothesis H3.1 and is biologically expected: HSP60 is a matrix-resident chaperonin, and its substrates should disproportionately localize to the matrix.

MTS prevalence was also modestly higher among HSP60 matrix substrates compared to non-HSP60 matrix proteins: 124/181 HSP60 matrix substrates had transit peptides (68.5%) versus 201/343 non-HSP60 matrix proteins (58.6%), yielding OR = 1.54 [1.05, 2.25], p = 0.029. This suggests that HSP60 substrates are somewhat more dependent on canonical MTS-mediated import.

### 5.3 MTS Architecture

Among the 436 proteins with both a detectable transit peptide and at least one CATH domain assignment, the MTS was overwhelmingly a pre-domain extension: 368/436 (84.4%) had MTS endpoints upstream of the first domain boundary, while only 68/436 (15.6%) showed MTS-domain overlap (binomial test, H0: p = 0.5, p = 3.42 x 10^-51) (Fig 4). This confirms hypothesis H3.3 with very high confidence.

For the 368 non-overlapping cases, the gap between MTS cleavage site and first domain start had a median of 18.0 residues (mean 37.4, range 0--579 residues). This gap represents the "pre-domain tail" -- a linker region that is neither part of the targeting signal nor part of the first structural domain. For the 68 overlapping cases, the mean overlap extent was 10.3 residues (maximum 48 residues).

### 5.4 MTS-Domain Relationship

The predominance of pre-domain MTS architecture has important implications for the folding pathway of mitochondrial matrix proteins. After import through the TIM23 translocase and MTS cleavage by MPP, the mature protein begins with an unstructured pre-domain tail followed by the first structural domain. This architecture means that the N-terminal domain emerges into the matrix in a largely unfolded state -- exactly the context in which HSP60 would be required to assist folding. The 18-residue median gap between MTS cleavage and first domain boundary provides a minimal linker that may serve as a flexible tether during chaperonin engagement. The 15.6% of cases where MTS overlaps the first domain suggests that in a minority of substrates, the MTS itself may partially constitute a structural element, though more commonly the targeting signal is fully excised before the first domain begins.

## 6. Integrated Findings

The three goals of this study converge on a nuanced picture of chaperonin-substrate relationships across the prokaryote-eukaryote divide (Fig 6).

**Goal 1** establishes that chaperonin substrates are enriched for specific structural folds -- particularly TIM barrels (OR = 22.6) and winged-helix domains (OR = 42.80) in the GroEL system -- but that the overall fold landscape is dominated by the same alpha-beta architectures (66--68%) that characterize their respective compartments. The mitochondrial carrier domain is systematically depleted from both chaperonin substrate sets, reflecting the alternative import pathway used by these inner membrane proteins. Cross-species fold conservation is high among orthologous substrates (55/69 sharing the same top superfamily), and Foldseek identifies 25 structurally similar clusters spanning both systems.

**Goal 2** reveals a universal N > C contact order asymmetry that is present in chaperonin substrates (GroEL: r = 0.387, p = 8.20 x 10^-4; HSP60: r = 0.519, p = 7.34 x 10^-6) but equally present in background proteins (matrix background: r = 0.388, p = 3.59 x 10^-4; full mitochondrial: r = 0.644, p = 9.00 x 10^-8). This asymmetry is not exaggerated in chaperonin substrates (all cross-dataset comparisons p > 0.57) and does not vary with GroEL dependence class (all p > 0.96). The most parsimonious interpretation is that the N > C contact order gradient is a general feature of multi-domain protein evolution, shaped by co-translational folding constraints rather than chaperonin biology per se.

**Goal 3** confirms that HSP60 substrates are strongly enriched for matrix localization (OR = 3.29, p = 1.60 x 10^-16) and that the mitochondrial targeting sequence is architecturally separate from the first structural domain in 84.4% of cases (p = 3.42 x 10^-51). The 18-residue median gap between MTS cleavage and domain start defines a structural "landing pad" that may facilitate chaperonin engagement upon import.

Taken together, these results suggest that chaperonin substrate identity is determined by specific fold topologies (TIM barrels, winged-helix motifs) rather than by global features of N-vs-C structural polarity. The conservation of both substrate identity and first-domain contact order across 2 billion years of divergence (r = 0.82 for ortholog pairs) argues that the structural determinants of chaperonin dependence are deeply embedded in protein fold architecture, transcending the prokaryote-eukaryote boundary.

## 7. Limitations

Several limitations of this analysis should be noted:

1. **pLDDT as confidence, not stability**. AlphaFold's pLDDT score reflects prediction confidence (how well the model agrees with the sequence-structure relationship in the training set), not thermodynamic stability. Regions with low pLDDT may be intrinsically disordered, flexible, or simply poorly predicted. We use pLDDT as a proxy for structural order, not for folding free energy. True stability comparisons would require experimental data (e.g., hydrogen-deuterium exchange, thermal melts).

2. **Co-IP interactome vs. functional substrate list**. The HSP60 interactome from Bruderer et al. is derived from SILAC co-immunoprecipitation, which captures physical interactions but cannot distinguish between obligate folding substrates and transient or non-productive interactions. Even the stringent Tier 1 filtering (266 proteins) likely includes some non-substrate interactors. By contrast, the GroEL substrate classification (Kerner et al., 2005) is based on pulse-chase experiments with chaperonin depletion, providing more direct evidence of folding dependence. This asymmetry in data quality between the two systems limits the strength of cross-species comparisons.

3. **Statistical power considerations**. With 62 total tests and hierarchical BH correction, some true positives may be missed due to the multiple testing burden. The superfamily enrichment analysis is particularly affected: testing many superfamilies dilutes power for individual folds. Effect sizes for the N-vs-C asymmetry comparisons (H2.2) were small (r = 0.019--0.110), and post hoc power analysis suggests that substantially larger sample sizes would be needed to detect asymmetry differences of this magnitude, if they exist.

4. **Relative contact order as a stability proxy**. Contact order correlates with folding rate (Plaxco et al., 1998) but is an imperfect proxy for folding difficulty. Proteins with identical contact orders can have very different folding landscapes depending on specific interactions, frustration patterns, and kinetic traps. Experimental folding data for a subset of proteins would strengthen the biological interpretation.

5. **MitoCarta version sensitivity**. The reclassification of 70 proteins between MitoCarta 2.0 and 3.0 demonstrates that sub-mitochondrial localization annotations are still evolving. Results that depend on fine-grained compartment assignments should be interpreted with this uncertainty in mind.

## 8. Key Numbers Summary Table

| Parameter | Value | Source |
|---|---|---|
| **Dataset sizes** | | |
| GroEL substrates | 252 proteins (38 Class I, 126 Class II, 84 Class III) | Kerner et al. 2005 |
| HSP60 Tier 1 substrates | 266 proteins | Bruderer et al., filtered |
| E. coli background (size-matched) | 663 proteins | UniProt K-12 |
| Matrix background | 524 proteins (471 with CATH) | MitoCarta 3.0 |
| Full mitochondrial proteome | 1,132 proteins (898 with CATH) | MitoCarta 3.0 |
| **Cross-species orthology** | | |
| RBH pairs | 40 (median identity 35.8%) | Module C |
| Orthogroups with substrates | 143 (34 shared, 51 GroEL-only, 58 HSP60-only) | Module C |
| Substrate pairs from orthogroups | 62 | Module C |
| Merged Dataset 6 | 69 unique cross-species pairs | Module C |
| **Domain architecture (Goal 1)** | | |
| CATH coverage: GroEL | 247/252 (98.0%) | Module E |
| CATH coverage: HSP60 | 241/266 (90.6%) | Module E |
| Alpha-beta class (GroEL / HSP60) | 66.3% / 67.5% | Module E |
| Winged helix enrichment (GroEL) | OR = 42.80 [5.62, 325.83], p_BH = 2.35 x 10^-6 | Module H |
| TIM barrel enrichment (GroEL) | OR = 22.6, p_BH = 2.35 x 10^-6 | Module H |
| Mito carrier depletion (HSP60) | OR = 0.16 [0.05, 0.52], p_BH = 2.79 x 10^-2 | Module H |
| Foldseek shared clusters | 25 (of 16,242 total) | Module E |
| Homolog pairs with shared top SF | 55/69 (79.7%) | Module H |
| **N-vs-C asymmetry (Goal 2)** | | |
| Contact order N > C (GroEL) | W = 2,262, p_BH = 8.20 x 10^-4, r = 0.387 | Module H |
| Contact order N > C (HSP60) | W = 1,776, p_BH = 7.34 x 10^-6, r = 0.519 | Module H |
| Contact order N > C (matrix bg) | W = 3,195, p_BH = 3.59 x 10^-4, r = 0.388 | Module H |
| Contact order N > C (mito bg) | W = 1,086, p_BH = 9.00 x 10^-8, r = 0.644 | Module H |
| Substrate vs. background (CO asym.) | All p_BH > 0.57, all r < 0.11 | Module H |
| GroEL class comparison (CO asym.) | H = 0.32, p_BH = 0.964, eta^2 = -0.014 | Module H |
| Cross-species RCO correlation | r = 0.82 (Spearman, 69 homolog pairs) | Module I (Fig 5) |
| **Matrix targeting (Goal 3)** | | |
| HSP60 matrix enrichment | OR = 3.29 [2.46, 4.40], p = 1.60 x 10^-16 | Module H |
| MTS prevalence (HSP60 vs bg) | OR = 1.54 [1.05, 2.25], p = 0.029 | Module H |
| MTS pre-domain fraction | 368/436 (84.4%), binomial p = 3.42 x 10^-51 | Module H |
| MTS-to-domain gap | Median 18.0 residues (mean 37.4, range 0--579) | Module G |
| MTS-domain overlap extent | Mean 10.3 residues (max 48) | Module G |
| HSP60 substrates with TP | 168/266 (63.2%) | Module G |
| **Statistical framework** | | |
| Total tests | 62 | Module H |
| Significant (hierarchical BH) | 45 | Module H |
| Goal families passing BH gate | 3/3 | Module H |
| Alpha level | 0.05 (two-sided) | Module H |
