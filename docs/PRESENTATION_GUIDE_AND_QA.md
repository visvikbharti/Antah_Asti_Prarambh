# Presentation Guide: Antah Asti Prarambh

## "The End is the Beginning" — A Comparative Structural Proteomics Study of Chaperonin Substrates

**Presenter:** Vishal Bharti, CSIR-IGIB, New Delhi
**Estimated Presentation Time:** 45-55 minutes (+ 15-20 min Q&A)
**Total Slides:** 34
**Presentation file:** Antah_Asti_Prarambh_Presentation_v3.pptx

---

## PART I: SLIDE-BY-SLIDE TALKING POINTS

### Presentation Structure
- Slides 1-4: Introduction, question, study design, methods
- Slides 5-7: Pipeline, modules, quality control
- Slides 8-11: Datasets and statistical framework
- Slides 12-15: Goal 1 -- Domain architecture
- Slides 16-20: Goal 2 -- N-vs-C stability (including key negative result)
- Slides 21-23: Goal 3 -- MTS targeting
- Slides 24-25: Cross-species conservation
- Slides 26-28: Statistics summary + sensitivity analysis
- Slides 29-30: Biological synthesis
- Slides 31-32: Limitations + future directions
- Slides 33-34: Takeaways + acknowledgments

---

### Slide 1: Title (~1 min)

**What to say:**
- "Good [morning/afternoon]. My project is called Antah Asti Prarambh, which is Sanskrit for 'The End is the Beginning.'"
- "This is a comparative structural proteomics study of two evolutionary related chaperonin systems: GroEL in E. coli and HSP60 in human mitochondria."
- "I'll explain why this name is meaningful as we go through the results."

**Tip:** Pause briefly after the Sanskrit title -- it creates natural curiosity.

---

### Slide 2: The Central Question (~2 min)

**What to say:**
- "Chaperonins are barrel-shaped complexes that help other proteins fold. About 10-15% of a cell's proteins need this help."
- "GroEL in E. coli and HSP60 in the mitochondrial matrix are evolutionary homologs separated by roughly 2 billion years -- since the endosymbiotic event that gave rise to mitochondria."
- "Our central question: what structural properties make a protein require chaperonin assistance, and are those properties conserved across this vast evolutionary distance?"
- "The project name refers to the idea that the C-terminus -- the 'end' of translation -- may be the 'beginning' of the chaperonin's work. We'll see that the reality is more nuanced."

**Key point:** The 2-billion-year evolutionary distance is a natural experiment in protein folding physics.

---

### Slide 3: Study Design (~2 min)

**What to say:**
- "We address three scientific goals. Goal 1: do chaperonin substrates have distinctive structural fold architectures? Goal 2: is there an N-vs-C stability asymmetry? Goal 3: how do mitochondrial targeting signals relate to the first structural domain?"
- "The study analyzes 25,007 proteins assembled into 7 curated datasets, with 62 pre-registered statistical tests grouped into 3 hypothesis families."
- "A critical design choice: compartment-matched controls throughout. We compare HSP60 substrates to the mitochondrial matrix proteome, not to the entire human proteome."

**Tip:** Emphasize the scale -- 25,007 proteins and 62 tests -- while noting that every test was planned in advance.

---

### Slide 4: Methods Overview (~1.5 min)

**What to say:**
- "Nine key methodological decisions were made upfront. Let me highlight three."
- "First, we use CATH structural domains, not InterPro sequence domains. For a study about folding, we need boundaries that reflect independent folding units."
- "Second, we do NOT use pLDDT as a stability metric. pLDDT is AlphaFold's confidence score, not thermodynamic stability. We use contact order for kinetics and FoldX for thermodynamics."
- "Third, hierarchical Benjamini-Hochberg correction for multiple testing -- BH within each of 3 families, then Simes across families."

**Key point:** This slide establishes rigor and credibility. Don't rush it.

---

### Slide 5: Computational Pipeline (~2 min)

**What to say:**
- "The pipeline has 9 modules, A through I. This flowchart shows the data flow from raw proteome downloads through to publication figures."
- "Starting from AlphaFold structure downloads, the pipeline branches: Foldseek for structural clustering, CATH plus Chainsaw for domain assignment, and FoldX for thermodynamic stability."
- "The analysis chain -- Modules F through I -- runs sequentially after domains are unified: stability metrics, then MTS analysis, then statistics, then figures."
- "All 25,007 proteins were processed: 24,530 by DSSP, 18,855 by CATH (51,667 domains), 6,164 by Chainsaw, and all 25,007 by FoldX."

---

### Slide 6: Analysis Modules (~1.5 min)

**What to say:**
- "This table maps every module to its key tools, inputs, and outputs." [Point to table]
- "Module A handles orthology -- MMseqs2 for reciprocal best hits and OrthoFinder-style all-vs-all search. Module B downloads and indexes AlphaFold structures."
- "Module E unifies CATH and Chainsaw domain assignments. Module F computes 14 metrics per protein region, including contact order and FoldX delta-G."
- "Module H runs all 62 statistical tests with hierarchical correction. Module I generates the 6 publication figures."

**Tip:** Don't read every row. Point to the table and highlight 2-3 modules that are most relevant to the audience.

---

### Slide 7: Quality Control & Robustness (~1.5 min)

**What to say:**
- "Quality control permeates the pipeline. Structure quality uses 5 tiers based on mean pLDDT, with automatic flags for potentially unreliable models. Only 4.6% of proteins are flagged."
- "Domain assignment has two independent methods with clear priority: CATH first (75.3% coverage), Chainsaw for the remainder (24.7%)."
- "Every DSSP run had a 30-second timeout. Every CATH API query was checkpointed. FoldX processed all 25,007 proteins with zero failures."
- "We also performed sensitivity analyses across multiple parameter choices -- more on that later in Slides 27-28."

---

### Slide 8: Datasets -- Chaperonin Substrates (~2 min)

**What to say:**
- "Our two substrate datasets required substantial curation."
- "GroEL: 252 substrates from Kerner et al. 2005, classified into three dependency classes. 149 of 252 accessions had been demerged in UniProt since publication -- we mapped each to the correct K-12 reference proteome entry."
- "HSP60: 266 Tier-1 substrates from Bie et al. 2020 (Cell Stress and Chaperones) SILAC co-IP. Filtered from 325 raw hits by excluding baits, co-chaperones, and contaminants, and requiring SILAC enrichment above 5. Median enrichment for Tier-1 is 22.2."
- "Class III GroEL substrates -- 84 proteins that absolutely require the chaperonin -- are a particularly interesting subgroup."

---

### Slide 9: Datasets -- Proteome Backgrounds (~1.5 min)

**What to say:**
- "Backgrounds are compartment-matched. GroEL substrates are compared against the full E. coli K-12 cytoplasmic proteome (4,403 proteins). HSP60 substrates are compared against the mitochondrial matrix proteome (525 proteins), not the full human proteome."
- "The mitochondrial proteome (1,136 proteins from MitoCarta 3.0) serves as the broader mito background."
- "A critical detail: between MitoCarta versions 2 and 3, 52 respiratory chain subunits were reclassified from 'matrix' to 'inner membrane.' This significantly affects the background composition."
- "All enrichment tests additionally use size-matched controls within 10 kDa bins."

---

### Slide 10: Datasets -- Cross-Species Homolog Pairs (~1.5 min)

**What to say:**
- "We identified 69 cross-species homolog pairs using two complementary methods: 40 from reciprocal best hits (MMseqs2 easy-rbh) and additional pairs from OrthoFinder-style analysis on full proteomes."
- "The union of both methods -- with evidence tracking -- gives us 69 unique pairs across 34 orthogroups containing substrates from both organisms."
- "An interesting observation: only 9.5% of Class III obligate substrates have HSP60 orthologs, versus ~19% for Classes I and II. The most GroEL-dependent proteins are the most evolutionarily divergent."
- "These 69 pairs are our window into 2 billion years of conservation."

---

### Slide 11: Statistical Framework (~2 min)

**What to say:**
- "We run 62 pre-registered tests organized into 3 hypothesis families: domain architecture, N-vs-C stability asymmetry, and MTS targeting."
- "Correction is hierarchical: Benjamini-Hochberg within each family, Simes method for family-level summaries, then BH across the 3 families. A test is significant only if it passes both levels."
- "Effect sizes are always reported alongside p-values: rank-biserial r for Mann-Whitney, eta-squared for Kruskal-Wallis, odds ratios with 95% CIs for Fisher's exact."
- "Result: 45 of 62 tests significant after hierarchical correction -- a 72.6% discovery rate."

**Tip:** Point to the test-type mapping table but don't read each row. Emphasize: "Every test was planned before we saw any results."

---

### Slide 12: Goal 1 -- CATH Class Distribution (~2 min)

**What to say:**
- "This figure shows CATH class distributions across all groups. Alpha-beta proteins dominate at 60-71%, but the distributions differ significantly: GroEL chi-squared p=5.2e-21, HSP60 p=2.4e-24."
- "GroEL substrates show higher beta-strand content (p=5.0e-7) and lower helix fraction (p=1.5e-5) compared to E. coli background. HSP60 substrates show higher helix fraction (p=1.7e-4)."
- "The DSSP secondary structure compositions are consistent with these CATH class differences."

**Key point:** The class distributions are significantly different from background, but the real story is the specific superfamily enrichments on the next slides.

---

### Slide 13: Goal 1 -- TIM Barrel Enrichment (~2.5 min)

**What to say:**
- "The headline result for domain architecture: GroEL substrates are 22.6 times more likely to contain TIM barrels than expected from the E. coli proteome background (Fisher's exact p=2.4e-21)."
- "Why TIM barrels? The barrel has 8-fold symmetry -- 8 beta-strands alternating with 8 alpha-helices. The barrel cannot form until ALL 8 strands are synthesized. This creates a vulnerable window during translation where the incomplete barrel can misfold or aggregate."
- "GroEL provides the protective barrel environment where this complex folding topology can form correctly."
- "This is not just about domain count -- neither substrate set is enriched for having MORE domains. It's about fold TYPE."

**Key point:** OR=22.6 is a very large effect size. This is robust across all sensitivity analyses.

---

### Slide 14: Goal 1 -- HSP60 Domain Enrichments + Cross-Species (~2 min)

**What to say:**
- "HSP60 shows its own distinct enrichment pattern: Rossmann-like folds (3.30.830.10, OR=5.4) and other matrix enzyme topologies (3.90.226.10, OR=4.8)."
- "When we compare across species, 79.7% of homolog pairs share the same CATH superfamily. A TIM barrel GroEL substrate has a TIM barrel HSP60 ortholog."
- "GroEL is also enriched in 1.10.10.10 folds (OR=50.9) -- mainly alpha-helical repeat structures that also have complex folding landscapes."
- "The overall picture: chaperonin substrate identity is determined by fold topology, and the specific topologies are largely conserved across 2 billion years."

---

### Slide 15: Goal 1 -- Domain Count Distribution (~1 min)

**What to say:**
- "This figure shows the distribution of domain counts per protein across substrate and background groups."
- "Neither GroEL nor HSP60 substrates are enriched for multi-domain proteins relative to their respective backgrounds."
- "The chaperonin dependence is about fold COMPLEXITY within individual domains, not about having more domains to fold."

**Key point:** This is a quick slide -- the negative finding (no domain count enrichment) supports the topology-driven narrative.

---

### Slide 16: Goal 2 -- N-vs-C Violin Plots (~2.5 min)

**What to say:**
- "These split violin plots show relative contact order for N-domains versus C-regions. In every group -- GroEL, HSP60, Matrix, Mito -- the N-domain distribution is shifted higher."
- "N-terminal domains have significantly higher contact order -- more topologically complex folds -- across all groups tested."
- "Effect sizes range from r=0.41 to r=0.48, all with p-values well below 10^-4."
- "But notice something: the pattern is present in ALL groups, not just substrates. The mito background shows the strongest effect."

**Key point:** Let the audience absorb the pattern across all groups before moving to the key result on the next slide.

---

### Slide 17: Goal 2 -- KEY RESULT: Universal N>C Asymmetry (~3 min)

**What to say:**
- "This is our most scientifically important finding, and it is a NEGATIVE result."
- "The N>C contact order asymmetry is universal -- present in all protein groups -- with the strongest signal in the mitochondrial background at p=7.1e-18."
- "When we test whether substrates show GREATER asymmetry than backgrounds using Mann-Whitney U: GroEL versus E. coli background gives p=0.058, HSP60 versus matrix gives p=0.536. Both non-significant."
- "Within GroEL substrates, Class III obligate substrates do NOT show greater asymmetry than Class I (Kruskal-Wallis p=0.77)."
- "This means chaperonins do not CREATE the asymmetry -- they EXPLOIT a pre-existing property of multi-domain protein architecture."

**Tip:** Emphasize that negative results are publishable and scientifically important. Don't be apologetic about them.

---

### Slide 18: Goal 2 -- Exploit vs Create (~2 min)

**What to say:**
- "So why does this universal asymmetry exist? It's physics."
- "Proteins are synthesized N-to-C on ribosomes. The N-terminus emerges first and has the most time to fold co-translationally. Evolution has placed the most topologically complex folds at the N-terminus, where co-translational folding has the best chance of working."
- "C-terminal regions, emerging last, adopt simpler folds or extend existing domains."
- "This is a fundamental property of protein architecture -- a 'gravitational constant' conserved across all of life, not caused by chaperonins. Chaperonins exploit this pre-existing landscape."

**Key point:** This interpretation reframes the entire project title: "The End is the Beginning" refers to physics, not biology.

---

### Slide 19: Goal 2 -- FoldX Thermodynamic Stability (~2 min)

**What to say:**
- "Contact order measures folding KINETICS. FoldX measures folding THERMODYNAMICS -- how stable the folded state is. Together they provide the complete picture."
- "With FoldX complete across all 25,007 proteins, the key finding is that thermodynamic differences are modest."
- "GroEL substrates show marginally lower stability than E. coli background (p=2.9e-3, Cohen's d=-0.07 -- a small effect). HSP60 substrates show no significant difference versus matrix background (p=0.80)."
- "Chaperonin dependence is driven primarily by fold topology and kinetic complexity, not by thermodynamic fragility."

**Key point:** The FoldX results reinforce the topology-driven narrative. Mention the species confound lesson if asked during Q&A.

---

### Slide 20: Goal 2 -- GroEL Class Comparison (~1.5 min)

**What to say:**
- "This figure compares the N-vs-C contact order difference across GroEL dependency classes (I, II, III)."
- "The boxplots are essentially indistinguishable. Kruskal-Wallis p=0.77 -- there is no class gradient."
- "Class III obligate substrates do not show greater N-vs-C asymmetry than Class I proteins that can fold spontaneously."
- "This further confirms that the asymmetry is a universal architectural property, not a marker of chaperonin dependence."

---

### Slide 21: Goal 3 -- MTS Targeting (~2.5 min)

**What to say:**
- "Shifting to Goal 3 -- how mitochondrial targeting signals relate to the first structural domain."
- "When nuclear-encoded mitochondrial proteins are imported, they pass through the TOM and TIM complexes in an unfolded state. The N-terminal transit peptide is cleaved upon arrival in the matrix."
- "This figure shows three panels: the targeting classification of HSP60 substrates, the gap histogram between MTS cleavage sites and first domain starts, and the scatter plot of MTS end versus domain start positions."
- "The key observation: MTS cleavage sites cluster before the first domain, creating a short unstructured linker."

---

### Slide 22: Goal 3 -- Pre-Domain Extension (~2 min)

**What to say:**
- "84.4% of transit peptides are separate pre-domain extensions -- the MTS ends before the first structural domain begins. This is overwhelmingly significant at p=3.4e-51."
- "The median gap is 18 residues. After import and cleavage, the protein has a short unstructured linker followed by an intact but unfolded first domain."
- "We call this the 'landing pad' model: the unfolded first domain is perfectly positioned for HSP60 to capture."
- "The 15.6% where MTS overlaps the first domain may have partially disrupted domains, making them even more dependent on chaperonin assistance."

**Key point:** This is the most statistically significant result in the entire study (p=3.4e-51).

---

### Slide 23: Goal 3 -- Matrix Enrichment (~1.5 min)

**What to say:**
- "HSP60 substrates are 3.29 times enriched for matrix localization compared to the general mitochondrial proteome (Fisher's exact p=1.6e-16)."
- "This makes biological sense -- HSP60 is a matrix chaperonin -- but the quantification is important: 46.6% are high-confidence matrix proteins, while 21.1% enter via non-canonical pathways without detectable transit peptides."
- "The matrix enrichment is the strongest targeting signal, confirming that HSP60 substrate identity is tightly coupled to compartment."

---

### Slide 24: Cross-Species Orthology (~2 min)

**What to say:**
- "This figure shows the cross-species analysis. We identified 69 homolog pairs using two complementary methods: reciprocal best hits and OrthoFinder-style orthogroup clustering."
- "422 orthogroups were identified across the full E. coli and human proteomes. 34 contain both GroEL and HSP60 substrates."
- "79.7% of homolog pairs share the same CATH superfamily -- the same fold topology is a chaperonin substrate in both organisms."
- "The limited cross-species overlap for Class III (9.5%) suggests the most dependent substrates are the most evolutionarily divergent."

---

### Slide 25: Cross-Species RCO Correlation (~2 min)

**What to say:**
- "This is the showpiece figure. N-domain relative contact order between GroEL substrates and their HSP60 orthologs shows Spearman r=0.82."
- "The topological complexity of the first structural domain has been conserved across 2 billion years of evolution."
- "Each point is a homolog pair. A TIM barrel in E. coli has a TIM barrel ortholog that is also an HSP60 substrate in humans."
- "This conservation argues that the structural properties driving chaperonin dependence are fundamental to protein folding physics, not specific to any one organism."

**Key point:** The r=0.82 correlation is the strongest evidence for evolutionary conservation of chaperonin substrate properties.

---

### Slide 26: Statistics Summary (~1.5 min)

**What to say:**
- "This figure summarizes all 62 statistical tests. 45 were significant after hierarchical BH correction -- a 72.6% discovery rate."
- "Family 1 (domain architecture): strong enrichments in specific topologies. Family 2 (N-vs-C): universal asymmetry confirmed, but NOT substrate-specific -- key negative results. Family 3 (MTS): strong targeting signals."
- "Effect sizes are reported for every test. The largest are MTS pre-domain (p=3.4e-51) and TIM barrel enrichment (OR=22.6, p=2.4e-21)."

---

### Slide 27: Sensitivity Analysis (text) (~1.5 min)

**What to say:**
- "We tested robustness across three key parameter choices."
- "HSP60 SILAC enrichment threshold varied from 3 to 10: all major findings remain significant at every threshold. N>C asymmetry p-values range from 10^-12 to 10^-20."
- "Size-matching bin width varied from 5 to 20 kDa: TIM barrel enrichment remains significant at all bin widths."
- "Background multiplier varied from 1x to 5x: enrichment tests remain significant across all levels."
- "No finding in this study depends on a particular parameter choice."

---

### Slide 28: Sensitivity Analysis -- Parameter Robustness (~1 min)

**What to say:**
- "This figure visualizes the parameter robustness. Each panel shows how a key result varies as one parameter changes."
- "The key message is stability: all findings are robust across SILAC thresholds (3-10), bin widths (5-20 kDa), and multipliers (1-5x)."
- "This is important for reviewers -- it demonstrates that our conclusions are not artifacts of arbitrary threshold choices."

**Tip:** This is a reassurance slide. If time is short, summarize in one sentence and move on.

---

### Slide 29: Biological Synthesis -- Exploit vs Create (~3 min)

**What to say:**
- "The three goals converge into a coherent narrative."
- "What makes a protein a chaperonin substrate? It's the FOLD -- specific topologies like TIM barrels that cannot fold co-translationally. Not domain count, not N-vs-C asymmetry, not obligate dependence class."
- "How does the mitochondrial system work? The MTS creates a landing pad -- after cleavage, the first domain emerges unfolded and HSP60 captures it. Matrix enrichment (OR=3.29) confirms compartment coupling."
- "What's conserved across 2 billion years? The fold topology itself. The N-domain complexity (r=0.82). The physics of vectorial translation."
- "So 'The End is the Beginning' has a nuanced meaning: the N>C asymmetry is real but UNIVERSAL. Chaperonins exploit a pre-existing biophysical landscape rather than creating substrate-specific features."

---

### Slide 30: Evolutionary Model -- 2 Billion Years (~2 min)

**What to say:**
- "This slide presents our evolutionary model. The ancestral bacterial chaperonin recognized proteins with complex fold topologies. After endosymbiosis, the mitochondrial HSP60 retained this specificity."
- "The substrate-defining properties -- TIM barrel enrichment, N-domain complexity, the pre-domain landing pad -- are conserved because they reflect fundamental protein folding physics."
- "The N>C asymmetry is a 'gravitational constant' of protein architecture. It exists because ribosomes synthesize proteins N-to-C, giving the N-terminus a co-translational folding advantage."
- "Evolution placed the most complex folds at the N-terminus. This is true for ALL multi-domain proteins, not just chaperonin substrates."

---

### Slide 31: Limitations (~1.5 min)

**What to say:**
- "I want to be transparent about limitations."
- "Methodologically: AlphaFold structures are predictions, not experimental. Co-IP captures interaction, not functional dependence. FoldX absolute delta-G values have limited accuracy -- relative comparisons are more reliable."
- "Statistically: 62 tests carry some residual false positive risk despite correction. The cross-species analysis has only 69 pairs, limiting power for subtle effects."
- "Biologically: we lack MTS prediction tools (TargetP, DeepMito) due to license constraints. 21.1% of HSP60 substrates enter the matrix by unknown mechanisms."

**Tip:** Acknowledging limitations upfront builds credibility and preempts tough questions.

---

### Slide 32: Future Directions (~1 min)

**What to say:**
- "With all analyses complete, the immediate next step is manuscript preparation with all 62 statistical tests and 6 publication figures."
- "Longer-term, we'd like to extend to Group II chaperonins (TRiC/CCT in the eukaryotic cytoplasm), obtain TargetP/DeepMito licenses for comprehensive MTS prediction, and develop a predictive model for chaperonin substrate identification based on fold topology."

---

### Slide 33: Key Takeaways (~2 min)

**What to say (slowly, one point at a time):**
1. "Chaperonin substrates are enriched for specific complex fold topologies -- TIM barrels for GroEL at odds ratio 22.6 (p=2.4e-21)."
2. "N-terminal domains universally have higher contact order (p=7.1e-18), but this is NOT substrate-specific (p=0.058, p=0.536). Chaperonins exploit a pre-existing biophysical property."
3. "Mitochondrial targeting signals create a 'landing pad' for post-import chaperonin engagement: 84.4% pre-domain extensions (p=3.4e-51), matrix enrichment OR=3.29 (p=1.6e-16)."

**Key point:** These three takeaways are your Q&A reference card. Know these numbers cold.

---

### Slide 34: Acknowledgments (~30 sec)

**What to say:**
- "Thank you. I'm grateful to my collaborators and to the developers of all the open-source tools that made this analysis possible."
- "I'm happy to take questions."

**Tip:** Keep it brief and warm. Move directly to Q&A.

---

## PART II: ANTICIPATED QUESTIONS AND ANSWERS

---

### Q1: Why not use pLDDT as a stability metric?

**Answer:** "pLDDT is AlphaFold's per-residue confidence score — the predicted Local Distance Difference Test. It measures how confident AlphaFold is about its prediction, not the thermodynamic stability of that protein region. A region with low pLDDT could be genuinely disordered OR it could be a stable domain that AlphaFold simply lacks training data for. Conflating pLDDT with stability is a widespread error in the field. We use contact order for folding kinetics and FoldX for thermodynamic stability. We do report pLDDT alongside these metrics but never equate it with stability."

---

### Q2: How did you handle the multiple testing problem with so many tests?

**Answer:** "We use hierarchical Benjamini-Hochberg correction. First, we group tests into three families based on the scientific question: domain architecture, stability asymmetry, and MTS targeting. We apply BH correction within each family, then compute a family-level summary p-value using the Simes method, and finally apply BH across the three families. A test is 'significant' only if both its within-family corrected p-value is below 0.05 AND its family passes the across-family correction. This controls the overall false discovery rate at 5% while preserving power for related tests. The final test count is 62, yielding 45 significant results (72.6% discovery rate) after hierarchical correction."

---

### Q3: If the N-vs-C asymmetry is universal, what's the novelty?

**Answer:** "That IS the novelty. Several papers have speculated that chaperonin substrates might have a unique N-vs-C folding asymmetry — that their C-terminal regions are specifically harder to fold, which is why they need chaperonins. Our study is the first to rigorously test this with proper controls and show that the asymmetry exists but is NOT substrate-specific. The negative result — that background proteins show the SAME asymmetry — is scientifically important because it redirects the field away from N-vs-C asymmetry as an explanation for chaperonin dependence and toward fold topology as the primary determinant. Negative results are publishable when they test important hypotheses with adequate power."

---

### Q4: How confident are you in the HSP60 substrate list?

**Answer:** "Reasonably confident, with caveats. The data comes from SILAC-quantified co-immunoprecipitation, which is the best available approach. Our Tier-1 set requires: (1) median SILAC enrichment ratio above 5, meaning the protein is at least 5 times more abundant in the HSP60 pull-down versus control, (2) MitoCarta 3.0 annotation confirming mitochondrial localization, (3) not being a bait, co-chaperone, or known contaminant. The median enrichment ratio for Tier-1 is 22.2 — very high confidence. However, co-IP captures physical interaction, not functional substrate dependence. Some of these 266 proteins may be transient interactors rather than true folding substrates. Unlike GroEL, there's no equivalent of the Kerner 2005 dependence classification for HSP60."

---

### Q5: Why use CATH instead of Pfam/InterPro for domain boundaries?

**Answer:** "Our study is fundamentally about protein FOLDING, which is a structural process. InterPro/Pfam domains are defined by sequence conservation — they group proteins that share detectable sequence homology. But the boundaries of a sequence domain don't necessarily correspond to where one independently folding structural unit ends and another begins. CATH domains are defined by structural classification of experimentally solved proteins. They represent genuine structural units. For computing contact order within a domain, we need these structural boundaries. Chainsaw, our secondary method, was specifically trained on CATH, ensuring consistency."

---

### Q6: Can you explain the 'landing pad' model?

**Answer:** "When a nuclear-encoded mitochondrial protein is imported, it passes through the TOM and TIM complexes in an unfolded, linear state. The N-terminal transit peptide (MTS) is cleaved by the mitochondrial processing peptidase upon arrival in the matrix. We found that in 84.4% of cases, the MTS is a separate extension that ends before the first structural domain begins — with a median gap of 18 residues. So after cleavage, the protein has a short unstructured linker followed by its entire first domain in an unfolded state. This creates an ideal 'landing pad' for HSP60 to grip the protein — the unfolded first domain is exposed and accessible. The 15.6% where MTS overlaps the first domain means those domains may be partially disrupted by cleavage, making them even more dependent on chaperonin assistance for re-folding."

---

### Q7: Why do Class III obligate substrates have fewer human orthologs?

**Answer:** "Only 9.5% of Class III proteins (8 out of 84) have HSP60 orthologs detected by RBH, compared to about 19% for Classes I and II. This likely reflects the evolutionary biology of these proteins. Class III substrates are the most GroEL-dependent — they literally cannot fold without the chaperonin barrel. Many of these are metabolic enzymes specific to E. coli's ecological niche (e.g., certain fermentation enzymes) that don't have equivalents in human mitochondria. The proteins that are conserved across species tend to be house-keeping enzymes (aminotransferases, dehydrogenases) that are needed in both organisms — and these are typically Class I or II."

---

### Q8: What does FoldX add that contact order doesn't capture?

**Answer:** "Contact order measures folding KINETICS — how fast or slow a protein folds. FoldX measures folding THERMODYNAMICS — how stable the folded state is. These are complementary. A protein could have high contact order (folds slowly) but weak thermodynamic stability (fragile fold) — or vice versa. With FoldX complete across all 25,007 proteins, we can now compare N-domain delta-G versus C-region delta-G, testing whether there's a thermodynamic asymmetry beyond the kinetic one. The compartment-matched analysis shows GroEL substrates have marginally lower stability than E. coli background (p=2.9e-3, d=-0.07, small effect), while HSP60 substrates show no significant difference versus matrix background (p=0.80)."

---

### Q9: How do you ensure your results aren't artifacts of protein size?

**Answer:** "Two ways. First, all enrichment tests use size-matched controls: for each substrate protein, we select background proteins from the same 10 kDa size bin. This prevents enrichment of a superfamily simply because substrates happen to be larger. Second, for within-protein comparisons (N-vs-C), size is controlled intrinsically — we compare two regions of the SAME protein, so the comparison is paired. We use the Wilcoxon signed-rank test, which is a paired non-parametric test that accounts for the within-protein structure of the data."

---

### Q10: What about disordered regions? How do they affect your analysis?

**Answer:** "Intrinsically disordered regions (IDRs) are handled in several ways. First, our three-region decomposition separates pre-domain tails (which may include disordered extensions) from structured domains. Contact order is computed only within structured domains, not across disordered linkers. Second, our quality validation (Module D4) flags proteins with low pLDDT, which often indicates disorder. Only 4.6% of proteins are flagged, and core substrate datasets have excellent quality (mean pLDDT 92-93). Third, we acknowledge that AlphaFold may not accurately model disordered regions, which is one reason we don't use pLDDT as a stability metric."

---

### Q11: Why not use molecular dynamics (MD) instead of FoldX?

**Answer:** "MD would be more rigorous but computationally prohibitive. A single 100-nanosecond simulation takes hours on a GPU for one protein; we have 25,007 proteins. FoldX is a rapid empirical force field that computes delta-G in about 40 seconds per protein, making full-proteome analysis feasible. The trade-off is that FoldX's absolute delta-G values have limited accuracy — but for comparative analysis (is the N-domain more or less stable than the C-region?), the relative ranking is what matters, and FoldX is reliable for that. If specific proteins show interesting patterns, targeted MD could be a follow-up."

---

### Q12: How does this compare to previous studies on chaperonin substrates?

**Answer:** "Previous structural studies of chaperonin substrates (Noivirt-Brik et al. 2007, Tartaglia et al. 2010) focused on sequence-derived properties: aggregation propensity, hydrophobicity, charge patterns. Our study is distinguished by three things: (1) we use STRUCTURAL domain boundaries from CATH/Chainsaw, not sequence domains; (2) we analyze the complete proteomes as background, not just the substrates; and (3) we explicitly test whether effects are substrate-specific versus universal using compartment-matched controls. The finding that N-vs-C asymmetry is universal, not substrate-specific, directly contradicts assumptions in some earlier work."

---

### Q13: The 84.4% pre-domain MTS result — could this be circular?

**Answer:** "Good question — are we simply detecting that transit peptides are at the N-terminus? No, because our test is more specific. We're not asking 'is the MTS at the N-terminus' (which is true by definition), but rather 'does the MTS END BEFORE the first structural domain STARTS?' The alternative is that the MTS overlaps with or includes part of the first domain. The 84.4% pre-domain result with a median gap of 18 residues tells us something biologically meaningful: the MTS is architecturally separate from the structural domains. Cleavage does not disrupt the first domain, leaving it intact but unfolded in the matrix."

---

### Q14: What's the statistical power for the cross-species analysis (69 pairs)?

**Answer:** "Limited, and we acknowledge this. For the Spearman correlation of N-domain contact order (r = 0.84, n = 45 pairs with paired values), the power is adequate — we detect a very strong correlation. However, for detecting subtle effects (say r = 0.2), we would need ~200 pairs. This means we can confidently say that STRONG conservation exists, but we may be missing weaker signals. The limited sample size comes from the biological reality: only 34 orthogroups contain substrates from both organisms. Expanding to Group II chaperonins (TRiC substrates in the cytoplasm) would increase sample size in future work."

---

### Q15: Why two 'phase2' directories?

**Answer:** "`workflow/phase2/` contains the pipeline CODE — all scripts, SLURM jobs, and configuration that define HOW the analysis runs on HPC. Think of it as the recipe. `results/phase2/` contains the OUTPUT DATA — all the results produced by running that pipeline, transferred back to the local Mac. Think of it as the finished dish. The separation follows standard bioinformatics practice: code and data live in different trees."

---

### Q16: How robust are the Chainsaw ML domain predictions?

**Answer:** "Chainsaw was published by Wells et al. in 2024 and trained directly on CATH domain structures. In our analysis, it serves as a SECONDARY method — we only use it for proteins that lack curated CATH annotations. The overlap between CATH and Chainsaw predictions for proteins that have both shows good agreement. At full scale, Chainsaw covers the 24.7% of proteins lacking CATH annotations, complementing CATH's 75.3% coverage (18,855 proteins). However, ML predictions inherently have some error rate, which is why we report the domain source (CATH vs Chainsaw) for every protein and prioritize CATH when available."

---

### Q17: Could the TIM barrel enrichment simply reflect E. coli's metabolic repertoire?

**Answer:** "Partly, yes — E. coli has many TIM barrel enzymes in its metabolic pathways. But our enrichment test explicitly controls for this. We compare GroEL substrates to the FULL E. coli cytoplasmic proteome, not to a random selection. An odds ratio of 22.6 means TIM barrels are 22.6 times MORE common among GroEL substrates than in the general E. coli proteome (p=2.4e-21). So even though E. coli has many TIM barrels, GroEL substrates have proportionally far MORE of them. The enrichment is real and robust — confirmed at both pilot and full scale, with the OR increasing from 8.4 (pilot) to 22.6 (full-scale) as the background became more representative."

---

### Q18: What's the Simes method and why use it?

**Answer:** "The Simes method provides a family-level summary p-value from a set of individual p-values. Given k ordered p-values p(1) <= p(2) <= ... <= p(k), the Simes p-value is min(k * p(i) / i). It's used in our hierarchical testing framework: after correcting tests within each family, we need a single number representing each family's overall significance to apply across-family correction. Simes is more powerful than Bonferroni for this purpose while still controlling the false discovery rate under independence or positive dependence of the individual tests."

---

### Q19: Why didn't you use TargetP or DeepMito for MTS prediction?

**Answer:** "TargetP 2.0 (from DTU) requires an academic license that we haven't yet obtained. Similarly, SignalP6 requires registration. Our current MTS analysis relies on UniProt's curated transit peptide annotations, which are experimentally validated but don't cover all proteins. About 30-40% of matrix proteins lack a detectable transit peptide annotation, which could use alternative import mechanisms. Obtaining the DTU license and running these predictors is planned for future work."

---

### Q20: What did FoldX reveal about N-domain vs C-region thermodynamic stability?

**Answer:** "FoldX is now complete for all 25,007 proteins. The key finding is that thermodynamic stability differences are modest. GroEL substrates show marginally lower total energy than E. coli background (p=2.9e-3, Cohen's d=-0.07 -- a small effect), while HSP60 substrates show no significant difference versus matrix background (p=0.80, NS). Critically, we discovered a species confound in the initial analysis: comparing GroEL substrates directly to HSP60 substrates yielded a spurious p=8.2e-47, which disappeared when we used compartment-matched backgrounds. This underscores why proper controls are essential. The overall picture is that chaperonin substrates are not dramatically less stable thermodynamically — their chaperonin dependence is driven more by fold topology (TIM barrels) and kinetic complexity (contact order) than by thermodynamic fragility."

---

### Q21: Why 8 Angstroms for contact order? Why not 6 or 10?

**Answer:** "8 Angstroms is the standard CA-CA distance cutoff established by Plaxco, Simons, and Baker in their original 1998 paper that defined relative contact order. This cutoff captures non-bonded contacts that contribute to the fold topology without being too permissive. At 6 Angstroms you miss important contacts; at 10 you include too many trivial near-neighbor interactions. The minimum sequence separation of 6 residues excludes local contacts within a single secondary structure element (an alpha-helix has ~3.6 residues per turn, so contacts within 5 residues are intra-helical). We used these exact values for consistency with the literature, and because the correlation between RCO and experimental folding rate (r ~ -0.75) was established with these parameters."

---

### Q22: Why MMseqs2 instead of BLAST for orthology?

**Answer:** "MMseqs2 is ~100-1000x faster than BLAST with comparable sensitivity. For our all-vs-all search of 4,403 x 20,416 proteins (~90 million comparisons), BLAST would take days; MMseqs2 completes in hours. We use sensitivity 7.5 (out of max 8), which is near-BLAST sensitivity. The easy-rbh mode implements the standard Reciprocal Best Hit algorithm. For the OrthoFinder-style analysis, we run bidirectional searches with explicit cutoffs (E-value < 1e-5, identity > 25%, coverage > 50%) and build orthogroups from reciprocal hits using union-find clustering."

---

### Q23: Why did you need the LD_LIBRARY_PATH fix on HPC?

**Answer:** "The IGIB HPC uses system gcc-8.4, which ships with an older version of libstdc++.so. Conda's numpy and scipy are compiled against a newer libstdc++. When Python tries to import numpy, it loads the system library first (because the system path comes before conda in LD_LIBRARY_PATH by default), and the version mismatch causes a silent ImportError. The fix is to prepend $CONDA_PREFIX/lib to LD_LIBRARY_PATH, ensuring conda's newer library is found first. Every SLURM script must include this line, or Python imports fail. This is a common issue on HPC clusters with older system compilers."

---

### Q24: Can you walk through the Foldseek parameters?

**Answer:** "Foldseek uses a dual alphabet: 3Di (3D interaction) and amino acid sequence. We use alignment type 2, which is the hybrid 3Di+AA mode — the most accurate. Sensitivity 7.5 is the highest we can use without quadratic time complexity (8 is exhaustive). E-value 0.001 for search and 0.01 for clustering balances sensitivity and specificity. Coverage 0.5 with bidirectional mode (mode 0) means at least 50% of BOTH query and target must align — this prevents matching a small domain against a large protein. Cluster mode 0 (set-cover cascaded) gives the most compact clustering. The 64 GB split-memory limit prevents out-of-memory errors on large databases."

---

### Q25: How exactly does the hierarchical BH correction work?

**Answer:** "It's a two-level procedure. At Level 1, within each of the three hypothesis families, we apply standard Benjamini-Hochberg FDR correction at alpha=0.05. For example, the N-vs-C stability tests are corrected together — if test i has the k-th smallest p-value among its siblings, its BH-adjusted p-value is p_k * (k_total/k). At Level 2, we compute a family-level summary using the Simes method: for each family's vector of BH-corrected p-values p(1) <= p(2) <= ... <= p(k), the Simes p-value is min over i of (k * p(i) / i). This gives one p-value per family. We then apply BH across the 3 family-level Simes p-values. A test is 'significant overall' only if it passes BOTH levels: within-family BH < 0.05 AND its family's Simes p-value passes the across-family BH."

---

### Q26: Why use Foldseek clustering in addition to CATH classification?

**Answer:** "They serve different purposes. CATH gives us a curated, hierarchical classification (Class.Architecture.Topology.Homologous superfamily) for enrichment testing. Foldseek gives us structure-based clustering that can find similarity between proteins with no detectable sequence homology (< 20% identity). Foldseek's 3Di alphabet captures local 3D arrangements that CATH doesn't — for example, two proteins with very different CATH codes might cluster together if they share a similar overall shape. The 24 'shared' Foldseek clusters (containing both GroEL and HSP60 substrates) confirm that the same 3D folds are chaperonin substrates in both organisms, complementing the CATH superfamily conservation analysis."

---

### Q27: What's the difference between the two phase2 directories?

**Answer:** "`workflow/phase2/` contains the pipeline CODE — all scripts, SLURM job definitions, Snakefile, and config.yaml that define HOW the analysis is executed. Think of it as the recipe. `results/phase2/` contains the OUTPUT DATA — all the TSV files, statistics, and figures produced by running that pipeline on HPC and transferred back to the local Mac. Think of it as the finished dish. This separation follows standard bioinformatics practice where code and data live in different directory trees, making it easy to version-control code separately from large data files."

---

### Q28: How do you handle proteins at the boundary between CATH and Chainsaw?

**Answer:** "CATH always takes priority because it's based on curated experimental structures. If a protein has a Gene3D annotation, we use those domain boundaries exclusively. Chainsaw is only invoked for proteins that LACK CATH annotations. In the unified assignment table, a 'source' column records whether each protein's domains came from CATH or Chainsaw, allowing sensitivity analyses. At full scale, 18,855 proteins (75.3%) have CATH annotations from Gene3D's InterPro entries, and the remaining 24.7% use Chainsaw predictions parsed from chopping strings."

---

### Q29: Why 50 proteins per FoldX array task?

**Answer:** "This is a balance between SLURM overhead and parallelism. With 25,007 proteins and 50 per task, we get 501 array tasks. Smaller chunks (e.g., 10 per task) would give 2,500+ tasks — more SLURM scheduling overhead and log file management. Larger chunks (e.g., 200 per task) would mean ~125 tasks but each takes too long, and a timeout would lose more work. At 50 proteins and ~40 seconds per protein, each task takes about 33 minutes, well within the 6-hour wall time. The QOS limit of 5 concurrent jobs means ~100 tasks run per day, giving an estimated 5 days for the full array."

---

### Q30: What happens to the 8 proteins without AlphaFold structures?

**Answer:** "P07203, P30042, P36969, Q16881, Q5THJ4, Q86UA3, Q9BVL4, and Q9NNW7 have no AlphaFold predicted structures. They are excluded from all structural analyses (DSSP, domain assignment, contact order, FoldX). However, they are retained in the dataset membership lists and included in sequence-level analyses where applicable. These 8 proteins without structures are excluded from the N-vs-C stability comparisons. This affects 8 out of 25,007 total proteins (0.03%) and has negligible impact on statistical power."

---

### Q31: What are the DSSP secondary structure findings?

**Answer:** "With 24,530 proteins processed by DSSP at full scale, we can now compare secondary structure compositions with high statistical power. GroEL substrates have significantly lower helix fraction (p=1.5e-5), higher strand fraction (p=5.0e-7), and higher coil fraction (p=1.9e-6) compared to E. coli cytoplasmic background. This is consistent with TIM barrel enrichment — TIM barrels are alpha/beta proteins with prominent beta-strands. HSP60 substrates show higher helix fraction (p=1.7e-4) and lower coil (p=2.2e-3) versus matrix background."

---

### Q32: Why did the TIM barrel OR change from 8.4 to 22.6?

**Answer:** "The full-scale CATH analysis (18,855 proteins via InterPro Gene3D) provides a comprehensive background covering the entire E. coli and human proteomes. With this complete background, TIM barrels are strongly disproportionately enriched in GroEL substrates, yielding an odds ratio of 22.6 (p=2.4e-21). Full-proteome backgrounds are essential for robust enrichment testing."

---

### Q33: How robust are the findings to parameter choices?

**Answer:** "We performed sensitivity analyses varying three parameters: (1) HSP60 SILAC enrichment threshold (3, 5, 7, 10) — N>C asymmetry significant at all thresholds; (2) size-matching bin width (5-20 kDa) — TIM barrel enrichment OR=3.4-5.2, significant at all widths; (3) background multiplier (1-5x) — enrichment significant at all levels. All key findings are robust to reasonable parameter variation."

---

### Q34: Why is the N-vs-C asymmetry NOT substrate-specific?

**Answer:** "We tested whether chaperonin substrates show greater N-vs-C contact order asymmetry than non-substrate backgrounds (Mann-Whitney U). GroEL vs E. coli background: p=0.058 (NS); HSP60 vs matrix: p=0.536 (NS). Background proteins show the same or even stronger N>C asymmetry. This means the asymmetry is a fundamental property of multi-domain protein architecture — likely reflecting co-translational folding constraints where N-terminal domains fold first under ribosomal pressure. Chaperonins exploit this pre-existing property rather than creating it. This is our most important negative result and redirects the field toward fold topology (TIM barrels, complex folds) as the primary determinant of chaperonin dependence."

---

### Q35: What about pLDDT as a stability metric?

**Answer:** "pLDDT is AlphaFold's per-residue prediction confidence, NOT thermodynamic stability. It measures how sure AlphaFold is about its prediction — a region with low pLDDT could be genuinely disordered OR a stable domain that AlphaFold lacks training data for. We use contact order (correlation r=-0.75 with experimental folding rates; Plaxco et al. 1998) as the primary folding kinetics proxy, and FoldX total energy for thermodynamic stability. pLDDT results are reported as model confidence comparisons only. This distinction is critical — conflating pLDDT with stability is a widespread error in the structural proteomics field."

---

## PART III: PRESENTATION TIPS

### General Advice
1. **Know your audience:** Adjust depth. For structural biologists, emphasize fold topology details. For cell biologists, emphasize the mitochondrial import narrative. For bioinformaticians, emphasize the statistical framework.

2. **The negative results are the story:** Many presenters are reluctant to highlight negative results. In this project, H2.2 (asymmetry not substrate-specific) and H2.3 (no class gradient) are the most publishable findings. Lead with confidence.

3. **The figures speak:** All 8 figures are self-contained with p-values and sample sizes. Let the audience read them before explaining.

4. **Time management:** Spend more time on Slides 16-18 (N-vs-C results, negative result, and interpretation) — these are the most novel and surprising findings. The domain architecture results (Slides 12-15) are important but less surprising.

5. **Sanskrit title as a hook:** The name creates curiosity. Explain it early (Slide 2) and revisit it in the synthesis (Slide 29-30).

### If Presenting to a Committee
- Emphasize the pre-registered hypotheses framework (Slide 11, shows rigor)
- Highlight the quality control slide (Slide 7, shows awareness)
- Stress the full reproducibility of the pipeline (Slides 5-6)
- Show the sensitivity analysis (Slides 27-28, robustness)

### If Presenting at a Conference (shorter, ~20 min)
- Focus on: Slides 1-3, 5(brief), 8(brief), 12-13, 16-18, 21-22, 24, 29, 33
- Skip: Slides 6-7, 9-10, 14-15, 19-20, 25-28, 30-32
- Spend more time on the biological synthesis (Slide 29)
- End with the key takeaways (Slide 33)

### If Presenting to PI/Collaborators (full detail, ~55 min)
- Present ALL 34 slides — they want to see the complete pipeline
- Spend extra time on Slides 5-7 (pipeline, modules, QC)
- Emphasize Slides 16-18 (N-vs-C negative result — most novel finding)
- Slides 27-28 (sensitivity analysis) shows thoroughness
- Slide 33 (key takeaways) is your Q&A reference card

### Key Numbers to Remember
- 25,007 proteins analyzed, 24,530 with DSSP, 18,855 with CATH domains
- 252 GroEL substrates, 266 HSP60 substrates
- 69 cross-species homolog pairs
- TIM barrel OR = 22.6 (GroEL, p=2.4e-21); CATH class: GroEL chi-sq p=5.2e-21, HSP60 p=2.4e-24
- N-vs-C contact order: r = 0.48 (strongest in mito background)
- NOT substrate-specific: p = 0.536 (HSP60 vs bg)
- NO class gradient: p = 0.77
- MTS pre-domain: 84.4%, median gap = 18 residues
- Cross-species conservation: r = 0.84
- 45/62 tests significant after hierarchical BH correction (72.6% discovery rate)

---

*Guide prepared: March 25, 2026; updated April 7, 2026 for v3 presentation (34 slides)*
*Project: Antah Asti Prarambh — Comparative Structural Proteomics of Chaperonin Substrates*
