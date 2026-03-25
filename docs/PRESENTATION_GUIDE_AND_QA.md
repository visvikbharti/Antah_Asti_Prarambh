# Presentation Guide: Antah Asti Prarambh

## "The End is the Beginning" — A Comparative Structural Proteomics Study of Chaperonin Substrates

**Presenter:** Vishal Bharti, CSIR-IGIB, New Delhi
**Estimated Presentation Time:** 55-65 minutes (+ 15-20 min Q&A)
**Total Slides:** 35 (comprehensive version with all technical parameters)

---

## PART I: SLIDE-BY-SLIDE TALKING POINTS

### Presentation Structure
- Slides 1-4: Introduction & Motivation
- Slides 5-6: Goals, Hypotheses, Datasets
- Slides 7-8: Methodological Decisions & Software Stack
- Slides 9-14: Phase 1 Workflow with exact tool parameters
- Slides 15-18: Module H Statistics, Phase 2 HPC Pipeline, FoldX, Figures
- Slides 19-25: Results with publication figures
- Slides 26-28: Synthesis, Limitations, Future Work
- Slides 29-32: Reproducibility, Session History, Bug Fixes, File Inventory
- Slides 33-35: Conclusions, Key Numbers, Thank You

---

### Slide 1: Title Slide (~1 min)

**What to say:**
- "Good [morning/afternoon]. My project is called Antah Asti Prarambh, which is Sanskrit for 'The End is the Beginning.'"
- "This is a comparative structural proteomics study comparing the substrates of two evolutionary related chaperonin systems: GroEL in E. coli and HSP60 in human mitochondria."
- "I'll explain why this name is meaningful as we go through the results."

**Tip:** Pause briefly after the title — the Sanskrit name creates natural curiosity.

---

### Slide 2: Outline (~30 sec)

**What to say:**
- "Here's what I'll cover today. We'll start with the biology, move through the computational methods, and then focus on three main result areas: domain architecture, N-vs-C asymmetry, and mitochondrial targeting."
- "I'll highlight both our positive findings AND our important negative results."

**Tip:** Don't read the outline — just acknowledge its structure and move on quickly.

---

### Slide 3: Background — Chaperonins (~3 min)

**What to say:**
- "Chaperonins are barrel-shaped protein complexes that help other proteins fold correctly. About 10-15% of proteins in a typical cell need this help."
- **Left column (GroEL):** "In E. coli, the GroEL/GroES system is the gold standard. Kerner et al. in 2005 classified 252 substrates into three dependency classes: Class I proteins can fold on their own but are accelerated by GroEL, Class II fold slowly without it, and Class III — 84 proteins — absolutely require GroEL and aggregate completely without it."
- **Right column (HSP60):** "In human mitochondria, HSP60 performs the same function. We identified 266 high-confidence substrates from Morten et al. 2020 using SILAC quantitative proteomics. These are nuclear-encoded proteins that must be imported into the mitochondrial matrix in an unfolded state, then re-fold with HSP60's help."
- "These two systems are evolutionary homologs separated by roughly 2 billion years — since the endosymbiotic event that gave rise to mitochondria."

**Key point to emphasize:** The 2-billion-year evolutionary distance is a natural experiment.

---

### Slide 4: Why Compare / Project Name (~2 min)

**What to say:**
- "Comparing these two systems tells us something fundamental. If the properties that make a protein need chaperonin help are conserved across 2 billion years, it means those properties are so fundamental to protein physics that evolution couldn't change them."
- "The project name refers to our central question: proteins are synthesized from N-terminus to C-terminus. The N-terminal region emerges first and folds first. The C-terminus — the 'end' of translation — is the last part to fold. Is it also the 'beginning' of the chaperonin's work?"
- "We'll see that the answer is more nuanced than we expected."

---

### Slide 5: Three Scientific Goals (~3 min)

**What to say:**
- **Goal 1:** "First, do chaperonin substrates have distinctive structural fold architectures? Not all protein folds are equally difficult — TIM barrels, for example, have complex 8-fold symmetry."
- **Goal 2:** "Second, is there an asymmetry between N-terminal and C-terminal regions? We measure this using contact order — a metric that captures how topologically complex a fold is."
- **Goal 3:** "Third, for HSP60 substrates specifically, how do mitochondrial targeting signals relate to the first structural domain? After import and signal cleavage, what does the protein look like as it encounters HSP60?"
- "We pre-registered 9 hypotheses across 3 families before running any tests, with hierarchical Benjamini-Hochberg correction for multiple testing."

**Tip:** Point to the hypothesis table but don't read each row. Say: "The details are here for reference — the key point is that every test was planned in advance."

---

### Slide 6: Seven Datasets (~2 min)

**What to say:**
- "Our analysis assembles seven carefully curated datasets." [Point to table]
- "The critical design choice is compartment-matched controls. We don't compare HSP60 substrates to the full human proteome — that would conflate 'properties of mitochondrial proteins' with 'properties of chaperonin substrates.' Instead, we compare within the mitochondrial matrix proteome."
- "Data cleaning was substantial: 149 of 252 GroEL accessions had been demerged in UniProt since the 2005 paper. For HSP60, we filtered from 325 raw hits to 266 Tier-1 substrates by excluding baits, co-chaperones, and contaminants, and requiring SILAC enrichment above 5."
- "A critical detail: between MitoCarta version 2 and 3, 52 respiratory chain subunits were reclassified from 'matrix' to 'inner membrane' — this significantly affects our background set."

---

### Slide 7: Critical Methodological Decisions (~2 min)

**What to say:**
- "Nine key decisions were made upfront, each with explicit rationale. Let me highlight three:"
- **Decision 2 (most important):** "We do NOT use pLDDT as a stability metric. pLDDT is AlphaFold's confidence score — it tells you how sure AlphaFold is about its prediction, not how thermodynamically stable the protein is. We use contact order for folding kinetics and FoldX for thermodynamics."
- **Decision 1:** "We use CATH structural domains, not InterPro sequence domains. For a study about folding, we need boundaries that reflect where one independently folding unit ends and another begins."
- **Decision 9:** "This one cost us hours of debugging: the bioconda package called 'stride' is a genomic variant caller, not the protein STRIDE program we needed. The correct one must be compiled from source."

**Tip:** This slide shows thoroughness and rigor — important for establishing credibility.

---

### Slide 8: Software & Tools Stack (~1 min)

**What to say:**
- "Here's our complete software stack. Everything is versioned and reproducible."
- "Key tools: MMseqs2 for sequence search, Foldseek for structural clustering, Chainsaw for ML domain prediction, FoldX for thermodynamic stability."
- "On the Mac, MMseqs2 and Foldseek run through Rosetta 2 emulation — x86_64 binaries on ARM."
- "On HPC, FoldX 5.1 and Chainsaw are compiled natively."

---

### Slide 9: Phase 1 Workflow (~2 min)

**What to say:**
- "Phase 1 has 9 modules, A through I. This table shows every script, input, output, and key parameter."
- "The pilot served three purposes: validate methods, produce preliminary results, and debug data integration issues like column name mismatches."
- Highlight the key scripts and their roles.

---

### Slide 10: Module C — Orthology Parameters (~2 min)

**What to say:**
- "For RBH, we used MMseqs2 easy-rbh, which gave us 40 one-to-one pairs."
- "For OrthoFinder-style orthology, we ran bidirectional all-vs-all search on the full proteomes — about 90 million comparisons — with specific cutoffs: E-value below 10^-5, identity above 25%, coverage above 50%, sensitivity 7.5."
- "Union-Find connected components gave us 422 orthogroups. 34 contain both GroEL and HSP60 substrates."
- "Merging both methods: 69 unique pairs, with evidence tracking which method found each pair."

---

### Slide 11: Module D — AlphaFold & DSSP Parameters (~1.5 min)

**What to say:**
- "AlphaFold download: we tried v6 first with v4 fallback. Batch size 50, 0.5 second delay between batches. 8 proteins had no AlphaFold model."
- "DSSP: mkdssp v2.2.1 with 30-second timeout. We group the 8 DSSP codes into 3 categories: helix (H,G,I), strand (E,B), and coil (everything else)."
- "Quality validation uses 5 tiers based on mean pLDDT, with 3 automatic flags for potentially unreliable structures."
- "For Phase 2, we used bulk FTP downloads — about 2 GB for E. coli and 20 GB for human."

---

### Slide 12: Module E — Domain Assignment Parameters (~2 min)

**What to say:**
- "CATH domains queried through InterPro Gene3D API at 1 request per second with checkpointing every 50 proteins."
- "Chainsaw uses STRIDE secondary structure as input — and this is where we learned the hard way that bioconda's 'stride' is a different tool entirely."
- "Foldseek parameters: sensitivity 7.5, E-value 0.001, 50% bidirectional coverage, 3Di+AA hybrid alignment. Phase 2 required 64 GB RAM for the all-vs-all search."
- "The coverage is excellent: 99.8% at pilot scale, 93.6% at full scale."

---

### Slide 13: Module F — Contact Order Parameters (~2.5 min)

**What to say:**
- "This is the core method slide for Goal 2. Contact order follows the Plaxco-Simons-Baker 1998 definition exactly."
- "We compute CA-CA distances with an 8 Angstrom cutoff, minimum 6 residues sequence separation. Relative contact order normalizes by protein length."
- "The three-region decomposition is critical: pre-domain tail, N-domain (first structural domain), C-region (everything after). Minimum 5 residues per region."
- "We compute 13 metrics per region: 7 sequence-based (charge, hydrophobicity, etc.), 3 structure-based (SS fractions, pLDDT), and 3 folding metrics (ACO, RCO, n_contacts)."
- "FoldX DeltaG will be added as the 14th metric after completion."

---

### Slide 14: Module G — MTS Analysis Parameters (~1.5 min)

**What to say:**
- "MTS analysis queries UniProt for transit peptide annotations using the REST API in batches of 100."
- "The gap calculation is straightforward: first_domain_start minus transit_peptide_end. Positive gap means the MTS is a separate pre-domain extension."
- "MitoCarta binary flags use case-insensitive string matching: 'Matrix' for is_matrix, 'MIM' or 'Membrane' for inner membrane."
- "The targeting classification has a clear hierarchy from high-confidence matrix to non-mitochondrial."

---

### Slide 15: Module H — Statistical Framework (~2 min)

**What to say:**
- "This is where rigorous statistics matter. We use hierarchical Benjamini-Hochberg: BH within each of the 3 families, then Simes method to get a family-level summary, then BH across families."
- "Effect sizes are always reported: rank-biserial r for Wilcoxon and Mann-Whitney, eta-squared for Kruskal-Wallis, odds ratios with 95% CIs for Fisher's exact, Cramer's V for chi-squared."
- "This table shows all test types mapped to specific hypotheses."

---

### Slide 16: Phase 2 HPC Pipeline — SLURM Resources (~2 min)

**What to say:**
- "This table shows every SLURM job with exact resource allocations. The two critical bottlenecks are Foldseek search at 64 GB RAM / 24 hours and Chainsaw at 72 hours."
- "Every script includes the LD_LIBRARY_PATH fix and python3 -u for unbuffered output — two hard-won lessons."
- "FoldX runs as a 501-task array job, manually submitted after the generation script."

---

### Slide 17: Pipeline Dependencies & FoldX Parameters (~2 min)

**What to say:**
- "The pipeline has three parallel branches after AlphaFold download: Foldseek, Chainsaw, and FoldX. The analysis chain is sequential: F -> H -> I."
- "FoldX parameters: 298.15 K, pH 7.0, ionic strength 0.05 M. Each protein goes through CIF-to-PDB conversion, RepairPDB to optimize rotamers, then Stability calculation."
- "Contact order gives us folding KINETICS; FoldX gives us folding THERMODYNAMICS. Together they provide the complete picture."

---

### Slide 18: Figure Generation Parameters (~30 sec)

**What to say:**
- "All figures use colorblind-friendly palettes, 300 DPI, with real p-values and sample sizes annotated."

---

### Slide 19-25: Results slides

**What to say:**
- "This is the Phase 2 dependency graph. Starting from AlphaFold download, the pipeline branches into three parallel tracks: Foldseek clustering, Chainsaw domain prediction, and FoldX stability calculations."
- "The analysis chain — Modules F through I — runs sequentially after domains are unified."
- "FoldX is the most resource-intensive: 501 array tasks, each processing 50 proteins at about 40 seconds per protein. It's currently about 42% complete."
- "The Foldseek search required 64 GB of RAM and 16 CPUs for the all-vs-all structural comparison of 25,000+ proteins."

---

### Slide 10: Structural Domain Assignment Method (~2 min)

**What to say:**
- "Domain assignment combines two approaches. CATH domains from Gene3D are our gold standard — curated, experimentally validated structural classifications. They cover 82.8% of our pilot proteins."
- "For the remaining proteins, Chainsaw — a deep learning model trained on CATH — predicts boundaries from 3D coordinates. This gives us 99.8% coverage."
- "At full scale, CATH covers the pilot set and Chainsaw adds 23,868 proteins from the full proteomes."
- "Foldseek provides structural clustering using a 3D structural alphabet, giving us 16,193 clusters across all proteins."

---

### Slide 11: Stability Metrics (~2 min)

**What to say:**
- "Contact order is our primary metric. It measures how many long-range contacts exist in a fold. A TIM barrel has high contact order because distant parts of the sequence must come together. High contact order correlates with slow folding."
- "FoldX provides the thermodynamic dimension — the actual stability of the folded state in kcal/mol."
- **Point to the three-region diagram:** "Every multi-domain protein is decomposed into three regions: the pre-domain tail (which may include the transit peptide), the first structural domain (N-domain), and everything after it (C-region). We compare N-domain versus C-region within each protein."
- "Critical note at the bottom: pLDDT is confidence, not stability. We report it but never equate it with thermodynamics."

---

### Slide 12: Results — Domain Architecture (~3 min)

**What to say:**
- **Point to Figure 1:** "Panel A shows CATH class distributions — alpha-beta proteins dominate all groups at 60-71%. Panel B shows the top superfamilies with clear enrichments. Panel C shows domain count distributions."
- "The key finding: GroEL substrates are 8.4 times more likely to contain TIM barrels than expected by chance, and 50.9 times more likely for winged helix domains. These are not random — they are topologically complex folds."
- "HSP60 shows different enrichments: Rossmann-like folds and other matrix enzyme topologies."
- "Importantly, neither substrate set is enriched for having MORE domains. It's not about domain count — it's about fold type."

---

### Slide 13: Domain Enrichment Details (~1.5 min)

**What to say:**
- "Here are the specific enrichments with Fisher's exact test p-values after BH correction."
- **Point to TIM barrel explanation box:** "Why do TIM barrels need chaperonins? The barrel has 8-fold symmetry — 8 beta-strands alternating with 8 alpha-helices. The barrel cannot form until ALL 8 strands are present. This means the protein must be fully synthesized before it can fold, creating a vulnerable window where it can aggregate."
- "GroEL provides the protective barrel environment where this complex folding can occur."

---

### Slide 14: N-vs-C Contact Order Results (~3 min)

**What to say:**
- **Point to Figure 2:** "Panel A shows split violin plots of relative contact order for N-domains versus C-regions. In every group — GroEL, HSP60, Matrix, Mito — the N-domain distribution is shifted higher."
- "The table shows this quantitatively. Effect sizes range from r = 0.41 to 0.48, all highly significant."
- **THE KEY SURPRISE:** "But look at the pattern: the effect is present in ALL groups, not just substrates. The mito background — which includes proteins that have nothing to do with HSP60 — shows the STRONGEST effect at p = 7.1 times 10 to the minus 18."
- "This tells us the N > C asymmetry is NOT a chaperonin-specific feature."

---

### Slide 15: Crucial Negative Results (~3 min)

**What to say:**
- "This slide presents what I consider our most important and scientifically surprising findings."
- **H2.2:** "When we compare the N-C difference BETWEEN substrates and background using Mann-Whitney U, there is NO significant difference. GroEL substrates have p = 0.058, HSP60 substrates have p = 0.536. Background proteins show the same asymmetry."
- **H2.3:** "Similarly, within GroEL substrates, Class III obligate substrates do NOT show greater asymmetry than Class I spontaneous folders. Kruskal-Wallis p = 0.77 — essentially zero effect."
- **Figure 3:** "These boxplots show the contact order difference by GroEL class. They're indistinguishable."
- "These negative results reshape the interpretation entirely: N-vs-C asymmetry is universal, not driven by chaperonin biology."

**Tip:** Emphasize that negative results are publishable and scientifically important. Don't be apologetic about them.

---

### Slide 16: Biological Interpretation of N-vs-C (~2 min)

**What to say:**
- "So why does this asymmetry exist if it's not about chaperonins?"
- "It's physics. Proteins are synthesized N-to-C on ribosomes. The N-terminus emerges first and has the most time to fold. Evolution has placed the most complex folds at the N-terminus, where co-translational folding has the best chance of working."
- "C-terminal regions, emerging last, adopt simpler folds or extend existing domains."
- "This is a fundamental property of protein architecture — a 'gravitational constant' that's always present, conserved across all of life, and not caused by chaperonins."

---

### Slide 17: MTS Targeting Results (~2.5 min)

**What to say:**
- **Figure 4, Panel A:** "HSP60 substrates are 3.3 times enriched for matrix localization compared to the general mitochondrial proteome. This makes biological sense — HSP60 is a matrix chaperonin."
- "Interestingly, 21.1% of HSP60 substrates enter the matrix by non-canonical pathways — no detectable transit peptide."
- **Panel B:** "The MTS gap histogram shows that most transit peptides end well before the first domain starts. The median gap is 18 residues."
- **Panel C:** "The scatter plot confirms: transit peptide cleavage sites cluster below and to the left of first domain starts."
- "84.4% of transit peptides are pre-domain extensions — overwhelmingly significant at p = 3.4 times 10 to the minus 51."
- "We call this the 'landing pad' model: after import and MTS cleavage, the protein has a short unstructured linker followed by an unfolded first domain — perfectly positioned for HSP60 to grab."

---

### Slide 18: Cross-Species Conservation (~2 min)

**What to say:**
- **Panel A:** "We identified 69 cross-species homolog pairs using both RBH and OrthoFinder."
- **Panel B — THE SHOWPIECE:** "N-domain contact order between GroEL substrates and their HSP60 orthologs shows Spearman r = 0.84, p = 5.3 times 10 to the minus 13. The topological complexity of the N-terminal domain has been conserved across 2 billion years of evolution."
- "79.7% of homolog pairs share the same CATH superfamily. A TIM barrel substrate in E. coli has a TIM barrel ortholog that's also an HSP60 substrate in humans."
- "An interesting observation: only 9.5% of Class III obligate substrates have HSP60 orthologs, versus ~19% for Classes I and II. The most GroEL-dependent proteins are the most divergent."

---

### Slide 19: Results Summary (~1 min)

**What to say:**
- "Figure 6 provides the visual summary."
- "Across 56 pre-registered tests with hierarchical BH correction, 25 were significant — a 44.6% discovery rate."
- "The three families all contributed: 9 domain architecture, 14 stability asymmetry, and 2 MTS targeting."

---

### Slide 20: Biological Synthesis (~3 min)

**What to say:**
- "The three goals converge into a coherent narrative."
- **Box 1:** "What makes a protein a chaperonin substrate? It's the FOLD. Specific topologies like TIM barrels that can't fold co-translationally. Not domain count, not N-vs-C asymmetry, not obligate dependence."
- **Box 2:** "How does the mitochondrial system work? HSP60 substrates are matrix proteins. The MTS creates a landing pad — after cleavage, the first domain emerges unfolded and HSP60 captures it."
- **Box 3:** "What's conserved across 2 billion years? The fold itself. The N-domain complexity. The physics of vectorial translation."
- **Blue box — The Title:** "So 'The End is the Beginning' has a nuanced meaning: the C-terminus IS where the chaperonin's work begins, but NOT because the C-terminus is more complex. The entire protein has its most complex domain at the N-terminus. The asymmetry is about translational physics, not chaperonin biology."

---

### Slide 21: Limitations (~1 min)

**What to say:**
- "I want to be transparent about limitations."
- "Methodologically: pLDDT as confidence (mitigated by CO + FoldX), co-IP as interaction not function (mitigated by SILAC), AlphaFold structures being predictions."
- "Statistically: 56 tests carry some false positive risk despite correction. Our cross-species analysis has only 69 pairs."
- "FoldX is 42% complete — the thermodynamic stability dimension will be added in the final analysis."

**Tip:** Acknowledging limitations upfront builds credibility and preempts tough questions.

---

### Slide 22: Future Directions (~1 min)

**What to say:**
- "Immediately after FoldX completes, we'll integrate thermodynamic stability and prepare the manuscript."
- "Longer-term, we'd like to extend to Group II chaperonins (TRiC/CCT in the eukaryotic cytoplasm) and develop a predictive model for chaperonin substrate identification."

---

### Slide 23: Software & Reproducibility (~30 sec)

**What to say:**
- "Everything is fully reproducible — 16 Python modules, 19 SLURM scripts, a central config file, and a deployment guide for the HPC pipeline."

---

### Slide 24: Conclusions (~2 min)

**What to say (slowly, one point at a time):**
1. "Chaperonin substrates are enriched for specific complex fold topologies — TIM barrels for GroEL at odds ratio 8.4."
2. "N-terminal domains universally have higher contact order — this is NOT substrate-specific."
3. "Mitochondrial targeting signals create a 'landing pad' for post-import chaperonin engagement."
4. "These properties are conserved across 2 billion years with remarkable correlation."
5. "The N>C asymmetry is a gravitational constant of protein architecture, not a feature of chaperonin biology."

---

### Slide 29: Reproducibility & Data Inventory (~1 min)

**What to say:**
- "Everything is organized in three code locations, with a central config.yaml for all parameters."
- "All raw data is from public databases — no proprietary dependencies."

---

### Slide 30: Session History (~1 min)

**What to say:**
- "The project developed across 6 sessions over 12 days. Session 4 was the most challenging — three major bugs discovered and fixed."

---

### Slide 31: Bug Fixes & Lessons (~2 min)

**What to say:**
- "These three bugs are worth highlighting because they represent real pitfalls in computational biology."
- "The STRIDE bug cost hours because there was no error message — Chainsaw simply returned ndom=0 for everything."
- "The column name bug affected 6 different lookups because each dataset uses different column names for the same field."
- "The HPC-specific lessons — LD_LIBRARY_PATH and Python unbuffered — are critical for anyone deploying similar pipelines."

**Tip:** These bugs show thorough debugging and build credibility with technical reviewers.

---

### Slide 32: File Inventory (~30 sec)

**What to say:**
- "This is a reference slide — every major file in the project with record counts."

---

### Slide 33: Conclusions (~2 min)

**What to say (slowly, one point at a time):**
1. "Chaperonin substrates are enriched for specific complex fold topologies — TIM barrels at odds ratio 8.4."
2. "N-terminal domains universally have higher contact order — but this is NOT substrate-specific."
3. "MTS creates a landing pad for post-import chaperonin engagement."
4. "These properties are conserved across 2 billion years with r = 0.84."
5. "The N>C asymmetry is a gravitational constant of protein architecture."

---

### Slide 34: Key Numbers (~30 sec)

**What to say:**
- "This reference table has every key number you might need during Q&A."

---

### Slide 35: Thank You (~30 sec)

**What to say:**
- "Thank you. I'm happy to take questions."

---

## PART II: ANTICIPATED QUESTIONS AND ANSWERS

---

### Q1: Why not use pLDDT as a stability metric?

**Answer:** "pLDDT is AlphaFold's per-residue confidence score — the predicted Local Distance Difference Test. It measures how confident AlphaFold is about its prediction, not the thermodynamic stability of that protein region. A region with low pLDDT could be genuinely disordered OR it could be a stable domain that AlphaFold simply lacks training data for. Conflating pLDDT with stability is a widespread error in the field. We use contact order for folding kinetics and FoldX for thermodynamic stability. We do report pLDDT alongside these metrics but never equate it with stability."

---

### Q2: How did you handle the multiple testing problem with so many tests?

**Answer:** "We use hierarchical Benjamini-Hochberg correction. First, we group tests into three families based on the scientific question: domain architecture (24 tests), stability asymmetry (30 tests), and MTS targeting (2 tests). We apply BH correction within each family, then compute a family-level summary p-value using the Simes method, and finally apply BH across the three families. A test is 'significant' only if both its within-family corrected p-value is below 0.05 AND its family passes the across-family correction. This controls the overall false discovery rate at 5% while preserving power for related tests. We also reduced the test count from 281 (Phase 1 pilot) to 56 (Phase 2) by focusing on pre-registered hypotheses."

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

### Q8: What will FoldX add that contact order doesn't capture?

**Answer:** "Contact order measures folding KINETICS — how fast or slow a protein folds. FoldX measures folding THERMODYNAMICS — how stable the folded state is. These are complementary. A protein could have high contact order (folds slowly) but weak thermodynamic stability (fragile fold) — or vice versa. The ideal chaperonin substrate might be one that folds slowly (high CO) but ends up very stable (strongly negative delta-G) — it needs help getting there, but once folded, stays folded. FoldX also allows us to compare N-domain delta-G versus C-region delta-G, testing whether there's a thermodynamic asymmetry beyond the kinetic one we've already measured."

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

**Answer:** "Chainsaw was published by Wells et al. in 2024 and trained directly on CATH domain structures. In our analysis, it serves as a SECONDARY method — we only use it for proteins that lack curated CATH annotations. For the pilot set, the overlap between CATH and Chainsaw predictions for proteins that have both shows good agreement. At full scale, the 93.6% assignment rate is comparable to CATH's coverage of the pilot set (82.8%). However, ML predictions inherently have some error rate, which is why we report the domain source (CATH vs Chainsaw) for every protein and prioritize CATH when available."

---

### Q17: Could the TIM barrel enrichment simply reflect E. coli's metabolic repertoire?

**Answer:** "Partly, yes — E. coli has many TIM barrel enzymes in its metabolic pathways. But our enrichment test explicitly controls for this. We compare GroEL substrates to the FULL E. coli cytoplasmic proteome, not to a random selection. An odds ratio of 8.4 means TIM barrels are 8.4 times MORE common among GroEL substrates than in the general E. coli proteome. So even though E. coli has many TIM barrels, GroEL substrates have proportionally far MORE of them. The enrichment is real and robust — confirmed at both pilot and full scale."

---

### Q18: What's the Simes method and why use it?

**Answer:** "The Simes method provides a family-level summary p-value from a set of individual p-values. Given k ordered p-values p(1) <= p(2) <= ... <= p(k), the Simes p-value is min(k * p(i) / i). It's used in our hierarchical testing framework: after correcting tests within each family, we need a single number representing each family's overall significance to apply across-family correction. Simes is more powerful than Bonferroni for this purpose while still controlling the false discovery rate under independence or positive dependence of the individual tests."

---

### Q19: Why didn't you use TargetP or DeepMito for MTS prediction?

**Answer:** "TargetP 2.0 (from DTU) requires an academic license that we haven't yet obtained. Similarly, SignalP6 requires registration. Our current MTS analysis relies on UniProt's curated transit peptide annotations, which are experimentally validated but don't cover all proteins. About 30-40% of matrix proteins lack a detectable transit peptide annotation, which could use alternative import mechanisms. Obtaining the DTU license and running these predictors is planned for future work."

---

### Q20: What would change if FoldX shows N-domains are LESS stable than C-regions?

**Answer:** "That would be a fascinating finding because it would create a 'kinetic paradox': N-domains fold MORE SLOWLY (higher contact order) but are LESS STABLE (higher delta-G). This could mean N-domains are kinetically trapped — they reach a local energy minimum during co-translational folding that isn't the global minimum. Chaperonins could help escape these traps. Conversely, if N-domains are MORE stable, it confirms the intuitive model that complex folds (high CO) produce stable structures (negative delta-G). Either way, it adds a new dimension to our understanding."

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

**Answer:** "It's a two-level procedure. At Level 1, within each of the three hypothesis families, we apply standard Benjamini-Hochberg FDR correction at alpha=0.05. For example, the 30 N-vs-C stability tests are corrected together — if test i has the k-th smallest p-value among its 30 siblings, its BH-adjusted p-value is p_k * (30/k). At Level 2, we compute a family-level summary using the Simes method: for each family's vector of BH-corrected p-values p(1) <= p(2) <= ... <= p(k), the Simes p-value is min over i of (k * p(i) / i). This gives one p-value per family. We then apply BH across the 3 family-level Simes p-values. A test is 'significant overall' only if it passes BOTH levels: within-family BH < 0.05 AND its family's Simes p-value passes the across-family BH."

---

### Q26: Why use Foldseek clustering in addition to CATH classification?

**Answer:** "They serve different purposes. CATH gives us a curated, hierarchical classification (Class.Architecture.Topology.Homologous superfamily) for enrichment testing. Foldseek gives us structure-based clustering that can find similarity between proteins with no detectable sequence homology (< 20% identity). Foldseek's 3Di alphabet captures local 3D arrangements that CATH doesn't — for example, two proteins with very different CATH codes might cluster together if they share a similar overall shape. The 24 'shared' Foldseek clusters (containing both GroEL and HSP60 substrates) confirm that the same 3D folds are chaperonin substrates in both organisms, complementing the CATH superfamily conservation analysis."

---

### Q27: What's the difference between the two phase2 directories?

**Answer:** "`workflow/phase2/` contains the pipeline CODE — all scripts, SLURM job definitions, Snakefile, and config.yaml that define HOW the analysis is executed. Think of it as the recipe. `results/phase2/` contains the OUTPUT DATA — all the TSV files, statistics, and figures produced by running that pipeline on HPC and transferred back to the local Mac. Think of it as the finished dish. This separation follows standard bioinformatics practice where code and data live in different directory trees, making it easy to version-control code separately from large data files."

---

### Q28: How do you handle proteins at the boundary between CATH and Chainsaw?

**Answer:** "CATH always takes priority because it's based on curated experimental structures. For the pilot set, if a protein has a Gene3D annotation, we use those domain boundaries exclusively. Chainsaw is only invoked for proteins that LACK CATH annotations. In the unified assignment table, a 'source' column records whether each protein's domains came from CATH or Chainsaw, allowing sensitivity analyses. For the 1,151 pilot proteins with CATH, the boundaries are from Gene3D's InterPro entries. For the 236 Chainsaw proteins, boundaries are parsed from chopping strings. In Phase 2, CATH proteins keep their annotations from Phase 1; the remaining 23,868 proteins use Chainsaw."

---

### Q29: Why 50 proteins per FoldX array task?

**Answer:** "This is a balance between SLURM overhead and parallelism. With 25,007 proteins and 50 per task, we get 501 array tasks. Smaller chunks (e.g., 10 per task) would give 2,500+ tasks — more SLURM scheduling overhead and log file management. Larger chunks (e.g., 200 per task) would mean ~125 tasks but each takes too long, and a timeout would lose more work. At 50 proteins and ~40 seconds per protein, each task takes about 33 minutes, well within the 6-hour wall time. The QOS limit of 5 concurrent jobs means ~100 tasks run per day, giving an estimated 5 days for the full array."

---

### Q30: What happens to the 8 proteins without AlphaFold structures?

**Answer:** "P07203, P30042, P36969, Q16881, Q5THJ4, Q86UA3, Q9BVL4, and Q9NNW7 have no AlphaFold predicted structures. They are excluded from all structural analyses (DSSP, domain assignment, contact order, FoldX). However, they are retained in the dataset membership lists and included in sequence-level analyses where applicable. For the 3 pilot proteins with no domain assignment (no CATH AND no Chainsaw because no AlphaFold structure), they are excluded from the N-vs-C stability comparisons. This affects 3 out of 1,390 pilot proteins (0.2%) and has negligible impact on statistical power."

---

## PART III: PRESENTATION TIPS

### General Advice
1. **Know your audience:** Adjust depth. For structural biologists, emphasize fold topology details. For cell biologists, emphasize the mitochondrial import narrative. For bioinformaticians, emphasize the statistical framework.

2. **The negative results are the story:** Many presenters are reluctant to highlight negative results. In this project, H2.2 (asymmetry not substrate-specific) and H2.3 (no class gradient) are the most publishable findings. Lead with confidence.

3. **The figures speak:** All 6 figures are self-contained with p-values and sample sizes. Let the audience read them before explaining.

4. **Time management:** Spend more time on Slides 14-16 (N-vs-C results and interpretation) — these are the most novel and surprising findings. The domain architecture results (Slides 12-13) are important but less surprising.

5. **Sanskrit title as a hook:** The name creates curiosity. Explain it early (Slide 4) and revisit it in the synthesis (Slide 20).

### If Presenting to a Committee
- Emphasize the pre-registered hypotheses framework (shows rigor)
- Highlight the 9 methodological decisions (shows awareness)
- Stress the full reproducibility of the pipeline
- Mention the pilot-to-full-scale validation (Phase 1 confirmed by Phase 2)

### If Presenting at a Conference (shorter, ~30 min)
- Skip Slides 8, 10-14, 16-18, 30-32 (too detailed for conference)
- Focus on: 1-5, 6(brief), 7(brief), 19-26, 33-35
- Spend more time on the biological synthesis (Slide 26)
- Emphasize the "gravitational constant" metaphor
- End with the conservation scatter plot (r = 0.84) — it's visually striking

### If Presenting to PI/Collaborators (full detail, ~60 min)
- Present ALL 35 slides — they want to see exact parameters
- Spend extra time on Slides 10-14 (tool parameters, rationale)
- Emphasize Slides 16-17 (HPC resources, FoldX parameters)
- Slide 31 (bug fixes) shows thoroughness and debugging rigor
- Slide 34 (key numbers) is your Q&A reference card

### Key Numbers to Remember
- 25,007 proteins analyzed
- 252 GroEL substrates, 266 HSP60 substrates
- 69 cross-species homolog pairs
- TIM barrel OR = 8.4 (GroEL), Winged helix OR = 50.9
- N-vs-C contact order: r = 0.48 (strongest in mito background)
- NOT substrate-specific: p = 0.536 (HSP60 vs bg)
- NO class gradient: p = 0.77
- MTS pre-domain: 84.4%, median gap = 18 residues
- Cross-species conservation: r = 0.84
- 25/56 tests significant after hierarchical BH correction

---

*Guide prepared: March 25, 2026*
*Project: Antah Asti Prarambh — Comparative Structural Proteomics of Chaperonin Substrates*
