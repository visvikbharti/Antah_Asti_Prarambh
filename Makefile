# =============================================================================
# Antah Asti Prarambh — Master Makefile
# =============================================================================
# Reproduces the full analysis pipeline from scratch.
#
# Quick start:
#   make setup          # Install dependencies + download data
#   make phase1         # Run full Phase 1 pilot analysis
#   make phase2-local   # Run Phase 2 locally (needs AlphaFold structures)
#
# On HPC:
#   make phase2-hpc     # Submit Phase 2 SLURM jobs
# =============================================================================

SHELL := /bin/bash
.DEFAULT_GOAL := help

# Project root (auto-detected)
PROJECT_DIR := $(shell pwd)

# Python (use conda env if available)
PYTHON := python3

# Phase 1 scripts
SCRIPTS := scripts
WORKFLOW := workflow/scripts

# Phase 2
PHASE2 := workflow/phase2
SLURM := workflow/phase2/slurm_jobs

# =============================================================================
# HELP
# =============================================================================

.PHONY: help
help:  ## Show this help message
	@echo ""
	@echo "Antah Asti Prarambh — Chaperonin Substrate Proteomics Pipeline"
	@echo "================================================================"
	@echo ""
	@echo "Setup:"
	@echo "  make setup             Install conda env + download external data"
	@echo "  make env               Create conda environment only"
	@echo "  make download-data     Download external datasets only"
	@echo ""
	@echo "Phase 1 (pilot, ~1,390 proteins):"
	@echo "  make phase1            Run full Phase 1 pipeline"
	@echo "  make phase0            Phase 0: data cleaning and curation"
	@echo "  make orthology         Module C: orthology and homology"
	@echo "  make structures        Module D: AlphaFold download + DSSP"
	@echo "  make domains           Module E: CATH + Chainsaw domains"
	@echo "  make stability         Module F: N-vs-C stability analysis"
	@echo "  make targeting         Module G: MTS targeting analysis"
	@echo "  make statistics        Module H: comparative statistics"
	@echo "  make figures           Module I: publication figures"
	@echo ""
	@echo "Phase 2 (full-scale, ~25,000 proteins):"
	@echo "  make phase2-hpc        Submit SLURM pipeline on HPC"
	@echo "  make phase2-local      Run Phase 2 locally (requires AlphaFold data)"
	@echo ""
	@echo "Utilities:"
	@echo "  make validate          Validate all Phase 1 outputs exist"
	@echo "  make presentation      Generate PowerPoint presentation"
	@echo "  make clean             Remove intermediate files"
	@echo ""

# =============================================================================
# SETUP
# =============================================================================

.PHONY: setup env download-data
setup: env download-data  ## Full setup: conda env + data download

env:  ## Create conda environment
	@echo "Creating conda environment from environment.yml..."
	conda env create -f environment.yml || conda env update -f environment.yml
	@echo ""
	@echo "Activate with: conda activate proteomics"

download-data:  ## Download external datasets
	bash scripts/download_external_data.sh

# =============================================================================
# PHASE 0: Data Cleaning and Curation
# =============================================================================

.PHONY: phase0
phase0: data/processed/groel_substrates_standardized.tsv data/processed/hsp60_tier1_substrates.tsv data/processed/human_mito_proteome.tsv  ## Phase 0: clean and curate all datasets

data/processed/groel_substrates_standardized.tsv: data/raw/custom/kerner_2005_groel_interactors_table_s3.csv
	$(PYTHON) $(SCRIPTS)/validate_uniprot_accessions.py

data/processed/hsp60_tier1_substrates.tsv: data/raw/custom/hsp60_interactome_clean.tsv
	$(PYTHON) $(SCRIPTS)/filter_hsp60_interactome.py

data/processed/human_mito_proteome.tsv: data/raw/mitocarta/Human.MitoCarta3.0.xls
	$(PYTHON) $(WORKFLOW)/parse_mitocarta.py

# =============================================================================
# MODULE C: Orthology / Homology
# =============================================================================

.PHONY: orthology
orthology: data/processed/groel_hsp60_homologs.tsv  ## Module C: RBH + orthogroups

data/processed/groel_substrates.fasta: data/processed/groel_substrates_standardized.tsv
	$(PYTHON) $(SCRIPTS)/module_c_extract_fasta.py

results/homology/rbh_groel_hsp60_annotated.tsv: data/processed/groel_substrates.fasta data/processed/hsp60_tier1_substrates.fasta
	$(PYTHON) $(WORKFLOW)/run_orthology.py

results/homology/rbh_groel_hsp60.tsv: results/homology/rbh_groel_hsp60_annotated.tsv
	$(PYTHON) $(SCRIPTS)/module_c_analyze_rbh.py

data/processed/groel_hsp60_homologs.tsv: results/homology/rbh_groel_hsp60_annotated.tsv results/homology/orthogroups_ecoli_human.tsv
	$(PYTHON) $(WORKFLOW)/build_dataset6_homologs.py

# =============================================================================
# MODULE D: Structure Acquisition + Quality
# =============================================================================

.PHONY: structures
structures: results/structures/dssp_summary.tsv results/structures/structure_quality_validation.tsv  ## Module D: AlphaFold + DSSP + quality

results/structures/structure_index.tsv: data/processed/groel_substrates_standardized.tsv data/processed/hsp60_tier1_substrates.tsv
	$(PYTHON) $(WORKFLOW)/download_alphafold_pilot.py

results/structures/dssp_summary.tsv: results/structures/structure_index.tsv
	$(PYTHON) $(WORKFLOW)/run_dssp.py

results/structures/structure_quality_validation.tsv: results/structures/structure_index.tsv results/structures/dssp_summary.tsv
	$(PYTHON) $(WORKFLOW)/validate_structure_quality.py

# =============================================================================
# MODULE E: Domain Architecture
# =============================================================================

.PHONY: domains
domains: results/domains/ml_domain_assignments.tsv results/domains/foldseek_clusters.tsv  ## Module E: CATH + Chainsaw + Foldseek

results/domains/cath_domain_assignments.tsv: results/structures/structure_index.tsv
	$(PYTHON) $(WORKFLOW)/get_cath_domains.py

results/domains/domain_structural_metrics.tsv: results/domains/cath_domain_assignments.tsv
	$(PYTHON) $(WORKFLOW)/compute_domain_structural_metrics.py

results/domains/chainsaw_domain_predictions.tsv: results/domains/cath_protein_summary.tsv
	$(PYTHON) $(WORKFLOW)/run_chainsaw_e2.py

results/domains/ml_domain_assignments.tsv: results/domains/cath_domain_assignments.tsv results/domains/chainsaw_domain_predictions.tsv
	@echo "Unified domain assignments built from CATH + Chainsaw"

results/domains/foldseek_clusters.tsv: results/structures/structure_index.tsv
	$(PYTHON) $(WORKFLOW)/analyze_foldseek.py

results/domains/domain_distribution_summary.tsv: results/domains/ml_domain_assignments.tsv
	$(PYTHON) $(WORKFLOW)/domain_distribution_summary.py

# =============================================================================
# MODULE F: N-vs-C Stability Analysis
# =============================================================================

.PHONY: stability
stability: results/termini/n_vs_c_paired_extended.tsv  ## Module F: N-vs-C paired comparisons

results/termini/n_vs_c_paired.tsv: results/domains/cath_domain_assignments.tsv results/structures/dssp_summary.tsv
	$(PYTHON) $(WORKFLOW)/module_f_n_vs_c_analysis.py

results/termini/n_vs_c_paired_extended.tsv: results/domains/ml_domain_assignments.tsv
	$(PYTHON) $(WORKFLOW)/module_f_extension_chainsaw.py

# =============================================================================
# MODULE G: Mitochondrial Targeting
# =============================================================================

.PHONY: targeting
targeting: results/mts/combined_targeting.tsv  ## Module G: MTS analysis

results/mts/combined_targeting.tsv: data/processed/hsp60_tier1_substrates.tsv results/domains/ml_domain_assignments.tsv
	$(PYTHON) $(WORKFLOW)/module_g_mts_analysis.py

# =============================================================================
# MODULE H: Comparative Statistics
# =============================================================================

.PHONY: statistics
statistics: results/stats/corrected_pvalues.tsv  ## Module H: hierarchical hypothesis tests

results/stats/corrected_pvalues.tsv: results/termini/n_vs_c_paired_extended.tsv results/mts/combined_targeting.tsv results/domains/domain_distribution_summary.tsv
	$(PYTHON) $(WORKFLOW)/module_h_comparative_stats.py

# =============================================================================
# MODULE I: Figures
# =============================================================================

.PHONY: figures
figures: results/stats/corrected_pvalues.tsv  ## Module I: publication figures
	$(PYTHON) $(WORKFLOW)/generate_figures.py

# =============================================================================
# FULL PHASE 1
# =============================================================================

.PHONY: phase1
phase1: phase0 orthology structures domains stability targeting statistics figures  ## Run complete Phase 1 pipeline
	@echo ""
	@echo "============================================================"
	@echo "  Phase 1 complete!"
	@echo "  Results in: results/"
	@echo "  Figures in: results/figures/"
	@echo "============================================================"

# =============================================================================
# PHASE 2 (HPC)
# =============================================================================

.PHONY: phase2-hpc phase2-local
phase2-hpc:  ## Submit Phase 2 SLURM pipeline on HPC
	@echo "Submitting Phase 2 pipeline via SLURM..."
	@echo "Edit workflow/phase2/config.yaml first to set your HPC paths."
	bash $(SLURM)/submit_pipeline.sh

phase2-local:  ## Run Phase 2 locally via Snakemake
	@echo "Running Phase 2 locally (requires AlphaFold structures downloaded)..."
	cd $(PHASE2) && snakemake --configfile config.yaml --cores 8

# =============================================================================
# UTILITIES
# =============================================================================

.PHONY: validate presentation clean

validate:  ## Validate all Phase 1 output files exist
	@echo "Checking Phase 1 outputs..."
	@PASS=0; FAIL=0; \
	for f in \
	  results/structures/structure_index.tsv \
	  results/structures/dssp_summary.tsv \
	  results/structures/structure_quality_validation.tsv \
	  results/domains/cath_domain_assignments.tsv \
	  results/domains/ml_domain_assignments.tsv \
	  results/domains/foldseek_clusters.tsv \
	  results/homology/rbh_groel_hsp60_annotated.tsv \
	  results/homology/orthogroups_ecoli_human.tsv \
	  data/processed/groel_hsp60_homologs.tsv \
	  results/termini/n_vs_c_paired.tsv \
	  results/termini/n_vs_c_paired_extended.tsv \
	  results/mts/combined_targeting.tsv \
	  results/stats/corrected_pvalues.tsv \
	  results/figures/fig1_domain_architecture.png; do \
	  if [ -f "$$f" ]; then \
	    echo "  [OK]   $$f"; \
	    PASS=$$((PASS + 1)); \
	  else \
	    echo "  [MISS] $$f"; \
	    FAIL=$$((FAIL + 1)); \
	  fi; \
	done; \
	echo ""; \
	echo "$$PASS passed, $$FAIL missing"

presentation:  ## Generate PowerPoint presentation
	$(PYTHON) create_presentation_v2.py

clean:  ## Remove intermediate/regenerable files
	@echo "Removing intermediate files..."
	rm -rf results/domains/_cath_checkpoint.json
	rm -rf results/homology/_mmseqs_ortho_work/
	rm -rf results/domains/foldseek/
	rm -rf __pycache__ workflow/scripts/__pycache__
	@echo "Done. Raw data, processed data, and final results preserved."
