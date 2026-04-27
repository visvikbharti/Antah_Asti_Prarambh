# Collaborator Onboarding Guide

**Welcome to the Antah Asti Prarambh project repository.**

This guide tells you exactly where to start, what to read, and in what order. The project analyzes 25,007 proteins across 7 datasets to answer three questions about chaperonin substrate recognition.

---

## Quick Start (5 minutes)

1. **Read the README.md** (project root) — overview, key findings, datasets, pipeline diagram
2. **Open the PowerPoint** — `Antah_Asti_Prarambh_Presentation_v3.pptx` (34 slides, ~15 min presentation)
   - Generate it: `python3 create_presentation_v3.py`

That's enough for a first meeting.

---

## Recommended Reading Order

### Level 1: The Story (30 minutes)

| Order | Document | What You'll Learn |
|:-----:|----------|-------------------|
| 1 | **[README.md](../README.md)** | Project overview, 3 goals, key findings, pipeline diagram |
| 2 | **[RESEARCH_FINDINGS_PHASE2.md](RESEARCH_FINDINGS_PHASE2.md)** | Complete manuscript-style findings (640 lines) — abstract, methods, all results with statistics, discussion |
| 3 | **[PRESENTATION_GUIDE_AND_QA.md](PRESENTATION_GUIDE_AND_QA.md)** | 35 anticipated questions with detailed answers |

### Level 2: The Methods (1 hour)

| Order | Document | What You'll Learn |
|:-----:|----------|-------------------|
| 4 | **[PRIMARY_HYPOTHESES.md](PRIMARY_HYPOTHESES.md)** | 9 pre-registered hypotheses, statistical framework |
| 5 | **[METHODS_AND_PROTOCOLS.md](METHODS_AND_PROTOCOLS.md)** | Exact commands, tool versions, reproducibility protocol |
| 6 | **[COLLABORATOR_PRESENTATION.md](COLLABORATOR_PRESENTATION.md)** | Slide-by-slide talking points and figure explanations |

### Level 3: The Data (as needed)

| Order | Document | What You'll Learn |
|:-----:|----------|-------------------|
| 7 | **[DATA_HANDOVER_INDEX.md](DATA_HANDOVER_INDEX.md)** | File-by-file guide to every result file, column descriptions |
| 8 | **[PHASE2_RESULTS_REPORT.md](PHASE2_RESULTS_REPORT.md)** | Detailed results with all statistical tests |
| 9 | **[COMPREHENSIVE_PROJECT_DOCUMENT.md](COMPREHENSIVE_PROJECT_DOCUMENT.md)** | Complete technical reference (~1,600 lines) |

### Level 4: Full Technical Deep-Dive (for reproducing)

| Order | Document | What You'll Learn |
|:-----:|----------|-------------------|
| 10 | **[DOCUMENTATION.md](DOCUMENTATION.md)** | Master technical docs, all modules A-I, HPC deployment |
| 11 | **[INSTALLATION.md](INSTALLATION.md)** | Setup guide for all tools |
| 12 | **[HPC_DEPLOYMENT_GUIDE.md](HPC_DEPLOYMENT_GUIDE.md)** | SLURM job scripts, HPC configuration |

---

## Where Are the Key Results?

### Statistical Results
```
results/phase2/stats/
├── corrected_pvalues_full.tsv        # All 62 tests with p-values, effect sizes, BH correction
├── statistics_summary_full.txt       # Human-readable summary report
├── supplementary_table_S1.tsv        # Publication-ready table of all tests
└── sensitivity_analysis.tsv          # Robustness across parameter choices
```

### Figures
```
results/phase2/figures/
├── fig1_domain_architecture.pdf      # CATH class distribution
├── fig2_n_vs_c_stability_full.pdf    # Contact order violin plots (N vs C)
├── fig3_foldx_deltag_comparison.pdf  # FoldX thermodynamic stability
├── fig3_groel_class_comparison.pdf   # GroEL class I/II/III effects
├── fig4_mts_targeting.pdf            # MTS-domain spatial relationship
├── fig5_orthology.pdf                # Cross-species conservation
├── fig6_summary.pdf                  # Combined key findings
└── fig8_sensitivity.pdf              # Parameter robustness analysis
```
All figures available in both PDF (vector) and PNG (300 DPI raster).

### Core Data Files
```
data/processed/
├── groel_substrates_standardized.tsv # 252 GroEL substrates (Kerner 2005)
├── hsp60_tier1_substrates.tsv        # 266 HSP60 substrates (Bie et al. 2020)
├── human_mito_proteome.tsv           # 1,136 mitochondrial proteins
├── human_matrix_proteome.tsv         # 525 matrix-only proteins
└── groel_hsp60_homologs.tsv          # 69 cross-species homolog pairs
```

---

## Three Key Findings to Know

### 1. Chaperonin substrates prefer specific folds
GroEL substrates are enriched in TIM barrels (OR=22.6, p=2.4e-21) and have distinctive secondary structure composition: higher beta-strand (p=5.0e-7) and coil (p=1.9e-6) fractions vs E. coli background. 24/24 domain architecture tests are significant.

### 2. N-terminal complexity is universal, NOT substrate-specific
N-terminal domains have higher contact order than C-terminal regions across ALL multi-domain proteins (p=7.1e-18). This is NOT specific to chaperonin substrates (substrate vs background: p>0.05). Chaperonins exploit this pre-existing architectural property.

### 3. Transit peptides are separate from structural domains
84.4% of mitochondrial transit peptides are pre-domain extensions (p=3.4e-51), with a median gap of 18 residues between MTS cleavage site and first structural domain.

---

## Important Scientific Notes

- **pLDDT is confidence, not stability.** AlphaFold's pLDDT score reflects prediction confidence. Contact order (Plaxco et al. 1998) is our folding kinetics proxy. FoldX provides thermodynamic stability.
- **FoldX was parameterized on experimental structures**, applied to AlphaFold predictions. Relative comparisons within the same pipeline are valid; absolute values should be interpreted with caution.
- **HSP60 data is co-IP based** (physical interaction), unlike GroEL data which has functional dependence classification (Kerner 2005). This limits direct "functional conservation" claims.
- **Sensitivity analysis confirms robustness** — all key findings hold across SILAC thresholds (3-10), size-matching bin widths (5-20 kDa), and background multipliers (1-5x).

---

## How to Reproduce

```bash
# Clone the repository
git clone https://github.com/visvikbharti/Antah_Asti_Prarambh.git
cd Antah_Asti_Prarambh

# Setup environment
make setup

# Run core analysis
make phase1

# Generate presentation
python3 create_presentation_v3.py
```

For full-scale HPC pipeline (25,007 proteins), see [HPC_DEPLOYMENT_GUIDE.md](HPC_DEPLOYMENT_GUIDE.md).

---

## Questions?

See [PRESENTATION_GUIDE_AND_QA.md](PRESENTATION_GUIDE_AND_QA.md) for 35 detailed Q&As covering methodology, statistics, biological interpretation, and limitations.

---

*Repository: https://github.com/visvikbharti/Antah_Asti_Prarambh (private)*
*Investigator: Vishal Bharti, CSIR-IGIB*
*Last updated: April 7, 2026*
