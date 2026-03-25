# Antah Asti Prarambh: Methods & Protocols

**Reproducible Protocol for Comparative Structural Proteomics of Chaperonin Substrates (GroEL/HSP60)**

Document version: 1.0
Generated: 2026-03-12
Project root: `/Users/vishalbharti/Downloads/Antah_Asti_Prarambh`

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Protocol for Each Module](#2-protocol-for-each-module)
   - [Module A: Data Acquisition & Cleaning](#module-a-data-acquisition--cleaning)
   - [Module B: Dataset Construction & Reference Data](#module-b-dataset-construction--reference-data)
   - [Module C: Orthology / Homology Layer](#module-c-orthology--homology-layer)
   - [Module D: Structure Acquisition & Indexing](#module-d-structure-acquisition--indexing)
   - [Module E: Structural Domain Architecture](#module-e-structural-domain-architecture)
   - [Module F: N-domain vs C-region Stability Analysis](#module-f-n-domain-vs-c-region-stability-analysis)
   - [Module G: Mitochondrial Targeting Analysis](#module-g-mitochondrial-targeting-analysis)
   - [Module H: Comparative Statistics](#module-h-comparative-statistics)
   - [Module I: Visualization & Figures](#module-i-visualization--figures)
3. [Data Flow Diagram](#3-data-flow-diagram)
4. [Statistical Analysis Details](#4-statistical-analysis-details)
5. [Figure Generation](#5-figure-generation)
6. [Quality Control Checklist](#6-quality-control-checklist)

---

## 1. Prerequisites

### 1.1 Hardware

| Resource | Minimum | Used in this study |
|----------|---------|-------------------|
| CPU | Apple M1 or x86_64 multi-core | Apple M1 (arm64) |
| RAM | 8 GB | 8 GB |
| Free disk | 20 GB | ~18 GB consumed |
| OS | macOS or Linux | macOS Darwin 25.2.0 |

### 1.2 Software Requirements (Exact Versions)

| Software | Version | Installation | Purpose |
|----------|---------|-------------|---------|
| Python | 3.9.16 | Anaconda base | All scripts |
| MMseqs2 | 18.8cc5c | conda (proteomics env) | Reciprocal best hits, orthology |
| Foldseek | 10.941cd33 | conda (proteomics env) | Structural clustering |
| mkdssp | (system) | conda base | Secondary structure assignment |
| Anaconda | 3 (2023+) | Installer | Environment management |

**Binary paths:**
```
Python:   /Users/vishalbharti/opt/anaconda3/bin/python3
MMseqs2:  /Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/mmseqs
Foldseek: /Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/foldseek
mkdssp:   /Users/vishalbharti/opt/anaconda3/bin/mkdssp
```

### 1.3 Python Packages

```
pandas>=2.2.2
scipy>=1.9.2
numpy
matplotlib
seaborn
biopython>=1.78
openpyxl
xlrd
requests
statsmodels
```

Install all at once:
```bash
pip install pandas scipy numpy matplotlib seaborn biopython openpyxl xlrd requests statsmodels
```

### 1.4 Conda Environments

Create the `proteomics` environment for MMseqs2 and Foldseek:
```bash
conda create -n proteomics python=3.9
conda activate proteomics
conda install -c bioconda mmseqs2 foldseek
```

### 1.5 Data Downloads Needed

| Data | Source | Download URL / Method |
|------|--------|----------------------|
| GroEL substrates (Kerner 2005) | Supplementary Table S3 | Already in `data/raw/custom/kerner_2005_groel_interactors_table_s3.csv` |
| HSP60 interactome (Bruderer 2020) | MOESM4 supplement | Already in `data/raw/custom/hsp60_interactome_clean.tsv` |
| Bruderer 2020 Excel supplement | Publisher | `data/raw/custom/12192_2020_1080_MOESM4_ESM.xlsx` |
| MitoCarta 3.0 | Broad Institute | `https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls` -> `data/raw/mitocarta/Human.MitoCarta3.0.xls` |
| E. coli K-12 proteome (FASTA) | UniProt UP000000625 | `data/raw/uniprot/ecoli_k12_proteome.fasta` |
| E. coli K-12 proteome (TSV) | UniProt UP000000625 | `data/raw/uniprot/ecoli_k12_proteome.tsv` |
| Human reference proteome (FASTA) | UniProt UP000005640 | `data/raw/uniprot/human_proteome.fasta` |
| Human reference proteome (TSV) | UniProt UP000005640 | `data/raw/uniprot/human_proteome.tsv` |
| AlphaFold structures | AlphaFold DB v6 (EBI) | Downloaded programmatically by Module D script |

### 1.6 Disk Space Estimate

| Component | Approximate Size |
|-----------|-----------------|
| Raw data (FASTA, TSV, XLS) | ~500 MB |
| AlphaFold CIF files (1,382 structures) | ~3 GB |
| DSSP output files (1,382 files) | ~200 MB |
| MMseqs2 working directories | ~2 GB |
| Foldseek databases and results | ~1.5 GB |
| Results (TSV, figures) | ~100 MB |
| **Total** | **~8 GB** |

---

## 2. Protocol for Each Module

---

### Module A: Data Acquisition & Cleaning

**Purpose:** Obtain and standardize all raw input datasets. Validate identifiers, resolve obsolete UniProt accessions, and produce clean TSV files with consistent column naming.

#### Step A1: Standardize GroEL substrates (Kerner 2005)

**Script:** `scripts/validate_uniprot_accessions.py`

**What it does:** Reads the 252 GroEL substrate entries from Kerner 2005 Table S3, resolves all UniProt accessions against the current database via REST API, follows demerge chains to find E. coli K-12 entries, and writes a standardized TSV.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 scripts/validate_uniprot_accessions.py
```

**Input files:**
- `data/raw/custom/kerner_2005_groel_interactors_clean.csv`
- `data/raw/custom/kerner_2005_groel_interactors_table_s3.csv`

**Output files:**
- `data/processed/groel_substrates_standardized.tsv` (252 rows)

**Key columns in output:** `original_accession`, `current_accession`, `accession_status`, `entry_name`, `gene_symbol`, `protein_name`, `organism`, `length`, `reviewed`, `groel_class` (I/II/III), `mass_kDa`, `subcellular_location`, `location_category`, `scop_folds`, `description_clean`

**Validation checks:**
- Verify 252 rows in output file
- All `current_accession` values should be valid 6-character UniProt IDs
- GroEL class distribution: Class I = 38, Class I or II = 4, Class II = 126, Class III = 84
- No duplicate accessions

**Troubleshooting:**
- If UniProt REST API returns 429 (rate limit), the script has built-in retry logic with backoff.
- If accessions have been further demerged since the last run, manual resolution may be needed for a handful of entries.

---

#### Step A2: Filter HSP60 interactome (Bruderer 2020)

**Script:** `scripts/filter_hsp60_interactome.py`

**What it does:** Reads the 325-protein HSP60 interactome from Bruderer 2020, imputes NDIC (Not Detected In Control) values at 2x the 95th percentile per SILAC replicate, assigns quality tiers, excludes bait/co-chaperone/contaminant proteins, and produces Tier 1 substrates.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 scripts/filter_hsp60_interactome.py
```

**Input files:**
- `data/raw/custom/hsp60_interactome_clean.tsv` (325 proteins)

**Output files:**
- `data/processed/hsp60_interactome_standardized.tsv` (325 rows, all proteins with flags)
- `data/processed/hsp60_tier1_substrates.tsv` (266 rows, high-confidence substrates)
- `data/processed/hsp60_filtering_report.txt`

**NDIC imputation values (2x 95th percentile):**
- `silac_ratio_coIP1`: 80.3386
- `silac_ratio_coIP2`: 80.6358
- `silac_ratio_coIP3`: 91.0113

**Tier criteria:**
- **Tier 1** (266): MitoCarta-positive AND median imputed SILAC ratio > 5
- **Tier 2** (49): SILAC > 2 but not MitoCarta+ or SILAC <= 5
- **Excluded** (10): Bait (HSPD1, HSPE1), co-chaperones (TRAP1, HSPA9, GRPEL1, DNAJA3), contaminants (IGHG2, TUBA1B, TUBB4B, TUBB)

**Validation checks:**
- 325 total input proteins
- 266 Tier 1 proteins in output
- 49 Tier 2 proteins
- 10 excluded proteins (2 bait + 4 co-chaperone + 4 contaminant)
- 100% of Tier 1 proteins are MitoCarta-annotated

---

#### Step A3: Parse MitoCarta 3.0

**Script:** `workflow/scripts/parse_mitocarta.py`

**What it does:** Parses the MitoCarta 3.0 Excel file and produces TSV files for the full human mitochondrial proteome and the matrix-only subset. Cross-references with the HSP60 interactome and documents MitoCarta 2.0 vs 3.0 discrepancies.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/parse_mitocarta.py
# Optionally, to download MitoCarta if not already present:
python3 workflow/scripts/parse_mitocarta.py --download
```

**Input files:**
- `data/raw/mitocarta/Human.MitoCarta3.0.xls`

**Output files:**
- `data/processed/human_mito_proteome.tsv` (1,136 proteins)
- `data/processed/human_matrix_proteome.tsv` (525 proteins)
- `data/processed/mitocarta_summary_report.txt`

**Key numbers:**
- Total MitoCarta 3.0 proteins: 1,136
- Matrix-localized: 523 (primary) / 525 (binary flag, includes multi-compartment)
- Inner membrane (MIM): 359
- Outer membrane (MOM): 110
- Intermembrane space (IMS): 51
- Unknown: 56

**Validation checks:**
- 1,136 rows in `human_mito_proteome.tsv`
- 525 rows in `human_matrix_proteome.tsv`
- 274 of 325 HSP60 interactors present in MitoCarta 3.0 (84.3%)
- 70 localization changes documented between MitoCarta 2.0 and 3.0
- 6 HSP60 interactors lost from MitoCarta 3.0 (C21ORF33, GAPDH, P4HB, RPS15A, RPS18, ATP5J2-PTCD1)

---

### Module B: Dataset Construction & Reference Data

**Purpose:** Download reference proteomes from UniProt and extract FASTA sequences for all substrate protein sets.

#### Step B1: Download UniProt reference proteomes

Download E. coli K-12 (UP000000625) and Human (UP000005640) reference proteomes from UniProt in both FASTA and TSV formats. Place them at:
- `data/raw/uniprot/ecoli_k12_proteome.fasta`
- `data/raw/uniprot/ecoli_k12_proteome.tsv`
- `data/raw/uniprot/human_proteome.fasta`
- `data/raw/uniprot/human_proteome.tsv`

#### Step B2: Extract FASTA sequences for substrate sets

**Script:** `scripts/module_c_extract_fasta.py`

**What it does:** Extracts FASTA sequences for GroEL substrates (from E. coli proteome) and HSP60 Tier 1 substrates (from Human proteome). Falls back to UniProt REST API for any accessions not found in the reference FASTA files.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 scripts/module_c_extract_fasta.py
```

**Input files:**
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `data/raw/uniprot/ecoli_k12_proteome.fasta`
- `data/raw/uniprot/human_proteome.fasta`

**Output files:**
- `data/processed/groel_substrates.fasta` (252 sequences)
- `data/processed/hsp60_tier1_substrates.fasta` (266 sequences)

**Validation checks:**
- 252 sequences in GroEL FASTA
- 266 sequences in HSP60 FASTA
- No duplicate headers in either file

**Troubleshooting:**
- 4 GroEL accessions (P69408, Q99390, P62593, P29368) are expected to be missing from the K-12 reference proteome and are fetched directly from UniProt REST API.

---

### Module C: Orthology / Homology Layer

**Purpose:** Identify orthologous protein pairs between E. coli GroEL substrates and human HSP60 substrates using both reciprocal best hits (RBH) and orthogroup analysis.

#### Step C1: MMseqs2 Reciprocal Best Hit (RBH) analysis

**What it does:** Finds pairs of proteins where protein A's best match in set B is protein B, and protein B's best match in set A is protein A. These are high-confidence 1:1 ortholog pairs.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh

# Run MMseqs2 easy-rbh
/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/mmseqs easy-rbh \
  data/processed/groel_substrates.fasta \
  data/processed/hsp60_tier1_substrates.fasta \
  results/homology/rbh_groel_hsp60.tsv \
  /tmp/mmseqs_tmp \
  --format-output "query,target,evalue,bits,alnlen,qcov,tcov,pident,qlen,tlen"
```

**Parameters:**
- MMseqs2 v18.8cc5c
- Default sensitivity
- E-value cutoff: 0.001 (MMseqs2 default for easy-rbh)
- Output format: query, target, evalue, bits, alignment length, query coverage, target coverage, percent identity, query length, target length

**Output file:**
- `results/homology/rbh_groel_hsp60.tsv` (40 pairs, headerless TSV)

**Key results:**
- 40 RBH pairs found
- 40/252 GroEL proteins matched (15.9%)
- 40/266 HSP60 proteins matched (15.0%)
- Median percent identity: 35.8%
- Median e-value: 3.39e-58

---

#### Step C2: Annotate RBH results

**Script:** `scripts/module_c_analyze_rbh.py`

**What it does:** Joins RBH pairs with GroEL class metadata and HSP60 gene names, computes statistics by GroEL class, and generates the summary report.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 scripts/module_c_analyze_rbh.py
```

**Input files:**
- `results/homology/rbh_groel_hsp60.tsv`
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `data/processed/groel_substrates.fasta`
- `data/processed/hsp60_tier1_substrates.fasta`

**Output files:**
- `results/homology/rbh_groel_hsp60_annotated.tsv` (40 rows)
- `results/homology/rbh_summary_report.txt`

**RBH pairs per GroEL class:**
- Class I: 7 (18.4% of 38)
- Class I or II: 1 (25.0% of 4)
- Class II: 24 (19.0% of 126)
- Class III: 8 (9.5% of 84)

**Validation checks:**
- 40 rows in annotated output
- All 40 pairs have non-null GroEL class and HSP60 gene annotations
- 8 Class III RBH pairs (the key finding: obligate GroEL substrates with conserved HSP60 dependence)

---

#### Step C3: Orthology group analysis (many-to-many)

**Script:** `workflow/scripts/run_orthology.py`

**What it does:** Runs MMseqs2 bidirectional all-vs-all search between the full E. coli and Human proteomes, filters by e-value/identity/coverage, identifies reciprocal pairs, then clusters them into orthogroups using union-find connected components. Intersects orthogroups with GroEL and HSP60 substrate lists.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/run_orthology.py
```

**Parameters:**
- E-value cutoff: 1e-05
- Minimum percent identity: 25.0%
- Minimum coverage (query or target): 50%
- Clustering: Union-find connected components on reciprocal pairs

**Input files:**
- `data/raw/uniprot/ecoli_k12_proteome.fasta`
- `data/raw/uniprot/human_proteome.fasta`
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `results/homology/rbh_groel_hsp60_annotated.tsv`

**Output files:**
- `results/homology/orthogroups_ecoli_human.tsv` (422 orthogroups)
- `results/homology/substrate_orthogroups.tsv` (143 substrate-containing orthogroups)
- `results/homology/orthology_comparison.tsv`
- `results/homology/orthology_summary_report.txt`
- `results/homology/_mmseqs_ortho_work/` (working directory)

**Key results:**
- 422 total orthogroups (both species)
- 199 one-to-one, 114 one-to-many (EC:HS), 40 many-to-one, 69 many-to-many
- 2,895 total reciprocal pairs, 453 strict RBH pairs
- 143 orthogroups containing chaperonin substrates
- 34 shared orthogroups (containing both GroEL and HSP60 substrates)
- 62 substrate pairs in shared orthogroups (vs 40 from simple RBH)
- 29 additional substrate pairs found beyond simple RBH

**Validation checks:**
- 422 orthogroups in output
- 34 shared substrate orthogroups
- 62 substrate pairs total
- 33 of the 40 original RBH pairs recovered in orthogroups

---

#### Step C4: Build unified homolog table (Dataset 6)

**Script:** `workflow/scripts/build_dataset6_homologs.py`

**What it does:** Merges evidence from RBH (40 pairs) and orthogroup analysis (62 substrate pairs) into a single unified homolog table with evidence source annotations.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/build_dataset6_homologs.py
```

**Output file:**
- `data/processed/groel_hsp60_homologs.tsv` (69 rows)

**Validation checks:**
- 69 unique homolog pairs in the final table
- Each pair has an `evidence` column indicating `rbh`, `orthogroup`, or `both`

---

### Module D: Structure Acquisition & Indexing

**Purpose:** Download AlphaFold predicted structures for all proteins in the analysis and compute secondary structure assignments using DSSP.

#### Step D1-D2: Download AlphaFold structures and build index

**Script:** `workflow/scripts/download_alphafold_pilot.py`

**What it does:** Collects unique UniProt accessions from all datasets (GroEL, HSP60, MitoCarta mito, matrix), downloads AlphaFold DB v6 CIF files (falling back to v4 if v6 unavailable), and builds a structure index with pLDDT summary statistics.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/download_alphafold_pilot.py
```

**Parameters:**
- AlphaFold DB URL pattern: `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif`
- Fallback: `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.cif`
- Batch size: 50
- Delay between batches: 0.5 seconds
- Per-request timeout: 30 seconds

**Input files:**
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `data/processed/human_mito_proteome.tsv`

**Output files:**
- `data/raw/alphafold/pilot/AF-{ACCESSION}-F1-model_v6.cif` (1,382 CIF files)
- `results/structures/structure_index.tsv` (1,382 rows; columns: `uniprot_accession`, `source_dataset`, `model_path`, `mean_plddt`, `min_plddt`, `max_plddt`, `n_residues`)

**Validation checks:**
- 1,382 CIF files in `data/raw/alphafold/pilot/`
- 1,382 rows in `results/structures/structure_index.tsv`
- All CIF files are non-empty and parseable

**Troubleshooting:**
- The EBI AlphaFold server may throttle downloads. The script handles retries automatically.
- A small number of accessions may not have AlphaFold models (e.g., very recently added or small peptides). These will be logged as missing.

---

#### Step D3: Run DSSP secondary structure assignment

**Script:** `workflow/scripts/run_dssp.py`

**What it does:** Runs `mkdssp` on every AlphaFold CIF file, parses the output to assign secondary structure codes (H/G/I = helix, E/B = strand, T/S/- = coil) to each residue, and produces per-protein summary statistics and per-residue tables.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/run_dssp.py
```

**Binary used:**
```
/Users/vishalbharti/opt/anaconda3/bin/mkdssp
```

**DSSP code grouping:**
- **Helix:** H (alpha-helix), G (3-10 helix), I (pi-helix)
- **Strand:** E (extended strand), B (isolated beta-bridge)
- **Coil:** T (turn), S (bend), `-` (loop/unassigned)

**Input files:**
- All CIF files in `data/raw/alphafold/pilot/`

**Output files:**
- `results/structures/dssp/{ACCESSION}.dssp` (1,382 individual DSSP files)
- `results/structures/dssp_summary.tsv` (per-protein SS fractions, pLDDT stats)
- `results/structures/dssp_per_residue.tsv` (per-residue SS assignments)

**Validation checks:**
- 1,382 DSSP files in `results/structures/dssp/`
- `dssp_summary.tsv` and `dssp_per_residue.tsv` both present and non-empty

---

### Module E: Structural Domain Architecture

**Purpose:** Assign CATH structural domains to all proteins, compute per-domain structural metrics, run Foldseek structural clustering, and generate domain distribution summaries.

#### Step E1: Obtain CATH domain assignments

**Script:** `workflow/scripts/get_cath_domains.py`

**What it does:** Queries the InterPro Gene3D REST API for each UniProt accession to retrieve CATH superfamily domain assignments. Uses checkpointing for resilience against network failures.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/get_cath_domains.py
```

**API endpoint:**
```
https://www.ebi.ac.uk/interpro/api/entry/cathgene3d/protein/uniprot/{accession}?format=json
```

**Rate limiting:** 1 request/second, with retry backoff of 5 seconds on failure.

**Input files:**
- `results/structures/structure_index.tsv` (1,382 accessions)

**Output files:**
- `results/domains/cath_domain_assignments.tsv` (2,141 domain assignments)
- `results/domains/cath_protein_summary.tsv` (1,390 rows = 1,382 queried)
- `results/domains/_cath_checkpoint.json` (checkpoint for resume)

**Key numbers by dataset:**
| Dataset | Proteins | With CATH domains | Total domains | Mean domains/protein |
|---------|----------|-------------------|---------------|---------------------|
| GroEL | 252 | 247 (98.0%) | 499 | 2.02 |
| HSP60 | 266 | 241 (90.6%) | 471 | 1.95 |
| Matrix | 524 | 471 (89.9%) | 945 | 2.01 |
| Mito | 1,132 | 898 (79.3%) | 1,628 | 1.81 |

**Validation checks:**
- 2,141 rows in `cath_domain_assignments.tsv` (header + 2,141 data rows = 2,142 lines)
- 1,390 rows in `cath_protein_summary.tsv` (header + 1,390 data rows = 1,391 lines)
- CATH superfamily codes follow the `X.XX.XX.XX` pattern

**Troubleshooting:**
- The InterPro API can be slow. The checkpoint file allows resuming from where it left off.
- If the API returns 503, wait and re-run; the checkpoint ensures no duplicate queries.

---

#### Step E2: Compute per-domain structural metrics

**Script:** `workflow/scripts/compute_domain_structural_metrics.py`

**What it does:** For each CATH domain assignment, extracts the corresponding residue range from the DSSP per-residue data and AlphaFold CIF files. Computes secondary structure fractions (helix, strand, coil) and pLDDT statistics (mean, min, fraction > 70, fraction > 90).

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/compute_domain_structural_metrics.py
```

**Input files:**
- `results/domains/cath_domain_assignments.tsv`
- `results/structures/dssp_per_residue.tsv`
- `results/structures/structure_index.tsv`
- CIF files in `data/raw/alphafold/pilot/`

**Output file:**
- `results/domains/domain_structural_metrics.tsv`

**Metrics computed per domain:**
- `mean_plddt`, `min_plddt`, `frac_plddt>70`, `frac_plddt>90`
- `frac_helix`, `frac_strand`, `frac_coil`

**Validation checks:**
- Output file has one row per domain (matching `cath_domain_assignments.tsv` row count)
- No NaN values in pLDDT columns for domains with valid CIF files

---

#### Step E3: Domain distribution summary

**Script:** `workflow/scripts/domain_distribution_summary.py`

**What it does:** Computes domain architecture distributions (domain count, top CATH superfamilies, CATH class proportions, most common domain architectures) for each of the four datasets (GroEL, HSP60, Matrix, Mito).

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/domain_distribution_summary.py
```

**Output files:**
- `results/domains/domain_distribution_summary.tsv`
- `results/domains/domain_distribution_report.txt`

**Key CATH class distributions:**
| Dataset | Alpha-beta | Mainly-alpha | Mainly-beta | Few-SS | Class-6 |
|---------|-----------|-------------|-------------|--------|---------|
| GroEL | 66.3% | 19.8% | 11.6% | 1.2% | 1.0% |
| HSP60 | 67.5% | 21.9% | 8.5% | 1.1% | 1.1% |
| Mito | 60.1% | 28.2% | 8.8% | 1.4% | 1.5% |

---

#### Step E4: Foldseek structural clustering

**What it does:** Creates a Foldseek structure database from all 1,382 AlphaFold CIF files, runs cascaded set-cover clustering, and an all-vs-all structural search. Analyzes shared clusters between GroEL and HSP60 substrates.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh

# Step 1: Create Foldseek database from CIF files
/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/foldseek createdb \
  data/raw/alphafold/pilot/ \
  results/domains/foldseek/pilot_db

# Step 2: Run cascaded clustering
/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/foldseek cluster \
  results/domains/foldseek/pilot_db \
  results/domains/foldseek/cluster_result \
  /tmp/foldseek_tmp \
  --cluster-mode 1 \
  --min-seq-id 0.3 \
  -c 0.5 \
  -e 0.01

# Step 3: Convert clustering to TSV
/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/foldseek createtsv \
  results/domains/foldseek/pilot_db \
  results/domains/foldseek/pilot_db \
  results/domains/foldseek/cluster_result \
  results/domains/foldseek/cluster_membership.tsv

# Step 4: Run all-vs-all search
/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/foldseek easy-search \
  results/domains/foldseek/pilot_db \
  results/domains/foldseek/pilot_db \
  results/domains/foldseek/search_results.tsv \
  /tmp/foldseek_search_tmp \
  -e 0.001 \
  --format-output "query,target,evalue,bits,alnlen,qcov,tcov,pident,qlen,tlen"
```

**Clustering parameters:**
- Method: Cascaded, set-cover (`--cluster-mode 1`)
- Minimum sequence identity: 0.3 (30%)
- Coverage threshold: 0.5 (50%)
- E-value: 0.01 (cluster), 0.001 (search)

**Step 5: Analyze clustering results**

**Script:** `workflow/scripts/analyze_foldseek.py`

```bash
python3 workflow/scripts/analyze_foldseek.py
```

**Output files:**
- `results/domains/foldseek/pilot_db*` (database files)
- `results/domains/foldseek/cluster_membership.tsv`
- `results/domains/foldseek/search_results.tsv`
- `results/domains/foldseek_clusters.tsv` (per-protein cluster table)
- `results/domains/foldseek_cluster_summary.txt`

**Key results:**
- 1,382 total proteins clustered
- 1,155 clusters total
- 999 singletons (86.5%), 149 small (2-5), 7 medium (6-20), 0 large (>20)
- 24 shared structural clusters (containing both GroEL and HSP60 members)
- 217 GroEL-only clusters, 218 HSP60-only clusters

**Validation checks:**
- 1,155 total clusters
- 24 shared clusters between GroEL and HSP60
- All-vs-all search: 11,664 total hits (including self), 10,282 non-self hits

---

### Module F: N-domain vs C-region Stability Analysis

**Purpose:** Decompose each multi-domain protein into three structural regions (pre-domain tail, first CATH domain = N-domain, remainder = C-region), compute sequence-derived and structure-derived metrics per region, and perform within-protein paired comparisons.

#### Step F1-F5: N-vs-C analysis

**Script:** `workflow/scripts/module_f_n_vs_c_analysis.py`

**What it does:**
1. Defines three regions per protein: pre-domain tail (residues before first CATH domain), N-domain (first CATH domain), C-region (everything after).
2. For single-domain proteins: compares N-half vs C-half of the domain.
3. Computes per-region metrics: length, net charge, hydrophobic fraction, SS fractions, pLDDT stats, relative contact order.
4. Contact order is computed from CA-CA distances in AlphaFold CIF files (distance cutoff for contacts, normalized by region length).

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/module_f_n_vs_c_analysis.py
```

**Input files:**
- `results/domains/cath_domain_assignments.tsv`
- `results/domains/cath_protein_summary.tsv`
- `results/structures/dssp_per_residue.tsv`
- `results/structures/structure_index.tsv`
- `data/processed/groel_substrates_standardized.tsv`
- CIF files in `data/raw/alphafold/pilot/`

**Output files:**
- `results/termini/n_vs_c_paired.tsv` (567 proteins with N-vs-C comparisons)
- `results/termini/contact_order.tsv`
- `results/termini/region_boundaries.tsv`
- `results/termini/sequence_metrics.tsv`
- `results/termini/structure_metrics.tsv`

**Amino acid property definitions used:**
- Hydrophobic: A, V, I, L, F, W, M
- Charged: K, R, D, E
- Polar: S, T, N, Q, Y, H, C
- Aromatic: F, W, Y
- Small: G, A, S

**Validation checks:**
- 567 proteins in `n_vs_c_paired.tsv` (= 568 lines including header)
- All paired comparisons have both N-domain and C-region metrics populated
- Contact order values are non-negative

---

### Module G: Mitochondrial Targeting Analysis

**Purpose:** Analyze mitochondrial targeting sequences (MTS), sub-mitochondrial localization, and the relationship between MTS and first structural domains.

#### Steps G1-G5: Targeting analysis

**Script:** `workflow/scripts/module_g_mts_analysis.py`

**What it does:**
1. **G1-G2:** Queries UniProt REST API for transit peptide and signal peptide annotations for all HSP60 and MitoCarta proteins.
2. **G3:** Cross-references with MitoCarta 3.0 sub-mitochondrial localization annotations.
3. **G4:** Integrates all targeting evidence into a combined classification table.
4. **G5:** Analyzes the relationship between MTS endpoints and first CATH domain boundaries.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/module_g_mts_analysis.py
```

**Input files:**
- `data/processed/hsp60_tier1_substrates.tsv`
- `data/processed/human_mito_proteome.tsv`
- `data/processed/human_matrix_proteome.tsv`
- `data/raw/uniprot/human_proteome.tsv`
- `results/domains/cath_domain_assignments.tsv`

**Output files:**
- `results/mts/uniprot_transit_signal_cache.tsv` (cached UniProt API results)
- `results/mts/combined_targeting.tsv` (1,139 proteins with targeting classification)
- `results/mts/mts_domain_relationship.tsv` (436 proteins with both TP and CATH domain)
- `results/mts/targeting_summary_report.txt`

**Targeting classification categories:**
| Category | Count | Percent |
|----------|-------|---------|
| Inner membrane (MIM) | 391 | 34.3% |
| High-confidence matrix | 325 | 28.5% |
| Non-canonical matrix import | 192 | 16.9% |
| Outer membrane (MOM) | 110 | 9.7% |
| Intermembrane space (IMS) | 52 | 4.6% |
| Mitochondrial (other/unspecified) | 50 | 4.4% |
| Other categories | 19 | 1.7% |

**MTS vs first domain relationship:**
- 436 proteins with both transit peptide and CATH domain
- MTS is pre-domain (non-overlapping): 368 (84.4%)
- MTS overlaps first domain: 68 (15.6%)
- Mean gap (domain_start - tp_end) for non-overlapping: 37.4 residues (median 18.0)
- Mean overlap extent: 10.3 residues (max 48)

**HSP60 substrate targeting:**
- HSP60 substrates with transit peptide: 168/266 (63.2%)
- HSP60 matrix substrates: 181 total, 124 with TP (68.5%)
- All matrix proteins: 524 total, 325 with TP (62.0%)

**Validation checks:**
- 1,139 proteins in `combined_targeting.tsv`
- 436 proteins in `mts_domain_relationship.tsv`
- 84.4% of MTS are pre-domain (non-overlapping)
- 260 HSP60 substrates overlap with MitoCarta

---

### Module H: Comparative Statistics

**Purpose:** Perform all statistical tests for the three primary goals (domain architecture enrichment, N-vs-C stability asymmetry, matrix targeting/MTS) with hierarchical multiple testing correction.

#### Step H1: Run comparative statistics

**Script:** `workflow/scripts/module_h_comparative_stats.py`

**What it does:** Brings together all previous analyses and runs the full battery of statistical tests organized by three goals:
1. **Goal 1 (Domain Architecture):** Fisher's exact tests for CATH superfamily enrichment, chi-squared tests for CATH class distribution.
2. **Goal 2 (N-vs-C Stability):** Wilcoxon signed-rank tests (within-protein paired), Mann-Whitney U (cross-dataset), Kruskal-Wallis (GroEL class comparison), Hotelling's T-squared (multivariate).
3. **Goal 3 (MTS/Targeting):** Fisher's exact test for matrix enrichment, binomial test for MTS pre-domain position.

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/module_h_comparative_stats.py
```

**Input files:**
- `results/domains/cath_protein_summary.tsv`
- `results/domains/cath_domain_assignments.tsv`
- `results/domains/domain_structural_metrics.tsv`
- `results/domains/foldseek_clusters.tsv`
- `results/domains/domain_distribution_summary.tsv`
- `results/termini/n_vs_c_paired.tsv`
- `results/termini/contact_order.tsv`
- `results/termini/region_boundaries.tsv`
- `results/mts/combined_targeting.tsv`
- `results/mts/mts_domain_relationship.tsv`
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `data/processed/groel_hsp60_homologs.tsv`
- `results/structures/structure_index.tsv`

**Output files:**
- `results/stats/domain_enrichment.tsv` (per-superfamily Fisher test results)
- `results/stats/stability_comparisons.tsv` (paired and between-group tests)
- `results/stats/targeting_stats.tsv` (targeting enrichment tests)
- `results/stats/corrected_pvalues.tsv` (all tests with hierarchical BH correction)
- `results/stats/statistics_summary_report.txt`

**Validation checks:**
- 281 total tests in `corrected_pvalues.tsv`
- 22 significant tests after hierarchical correction
- 3 goal families all significant at family level
- Domain architecture: 8 tests significant
- N-vs-C stability: 11 tests significant
- Matrix targeting: 3 tests significant

(See Section 4 for full statistical details.)

---

### Module I: Visualization & Figures

**Purpose:** Generate all publication-quality figures summarizing the project findings.

#### Step I1: Generate all figures

**Script:** `workflow/scripts/generate_figures.py`

**Commands:**
```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/generate_figures.py
```

**Output files** (each in PDF + PNG at 300 DPI):
- `results/figures/fig1_domain_architecture.{pdf,png}`
- `results/figures/fig2_n_vs_c_stability.{pdf,png}`
- `results/figures/fig3_groel_class_comparison.{pdf,png}`
- `results/figures/fig4_mts_targeting.{pdf,png}`
- `results/figures/fig5_orthology.{pdf,png}`
- `results/figures/fig6_summary.{pdf,png}`

**Validation checks:**
- 6 figures generated (12 files total: 6 PDF + 6 PNG)
- All PNG files are non-empty
- All PDF files contain Type 42 fonts (for journal compatibility)

(See Section 5 for style and color details.)

---

## 3. Data Flow Diagram

```
                         RAW INPUTS
                             |
     +--------------------------------------------------+
     |                       |                           |
 Kerner 2005         Bruderer 2020              MitoCarta 3.0
 GroEL Table S3      HSP60 interactome          Human.MitoCarta3.0.xls
     |                       |                           |
     v                       v                           v
 +-----------+      +-----------------+        +-------------------+
 | Module A1 |      |    Module A2    |        |    Module A3      |
 | validate_ |      | filter_hsp60_  |        | parse_mitocarta   |
 | uniprot   |      | interactome    |        |                   |
 +-----------+      +-----------------+        +-------------------+
     |                       |                     |           |
     v                       v                     v           v
 groel_subs_       hsp60_tier1_              mito_proteome  matrix_proteome
 standardized      substrates               (1,136)        (525)
 (252)             (266)                         |           |
     |                       |                   |           |
     +-------+-------+      +-------+-----------+-----------+
             |               |       |
             v               v       v
     +------------------+   +-----------------------+
     |    Module B       |   | UniProt proteomes     |
     | Extract FASTA     |   | E. coli + Human FASTA |
     | sequences         |   |                       |
     +------------------+   +-----------------------+
         |       |               |           |
         v       v               v           v
    groel.fasta  hsp60.fasta  ecoli.fasta  human.fasta
         |       |               |           |
         +---+---+               +-----+-----+
             |                         |
             v                         v
     +-------------------+    +-------------------+
     | Module C (C1-C2)  |    | Module C (C3)     |
     | MMseqs2 easy-rbh  |    | Orthology groups  |
     | 40 RBH pairs      |    | 422 orthogroups   |
     +-------------------+    +-------------------+
             |                         |
             +------+    +-------------+
                    |    |
                    v    v
            +------------------+
            |   Module C (C4)  |
            |   Unified        |
            |   homolog table  |
            |   69 pairs       |
            +------------------+
                    |
                    v
     +-----------------------------------------+
     |  All 1,382 unique UniProt accessions    |
     +-----------------------------------------+
                    |
                    v
     +----------------------------+
     |       Module D (D1-D2)     |
     | Download AlphaFold v6 CIF  |
     | 1,382 structures           |
     +----------------------------+
                    |
                    v
     +----------------------------+
     |       Module D (D3)        |
     | Run DSSP (mkdssp)          |
     | 1,382 DSSP files           |
     +----------------------------+
                    |
          +---------+----------+
          |                    |
          v                    v
 +-------------------+  +-------------------+
 | Module E (E1-E3)  |  | Module E (E4)     |
 | CATH domains       |  | Foldseek cluster  |
 | 2,141 assignments |  | 1,155 clusters    |
 | Domain metrics    |  | 24 shared         |
 +-------------------+  +-------------------+
          |                    |
          +---+---+---+--------+
              |   |   |
              v   v   v
 +------------------+  +-------------------+  +-------------------+
 | Module F          |  | Module G          |  |                   |
 | N-vs-C analysis  |  | MTS targeting     |  |                   |
 | 567 paired       |  | 1,139 classified  |  |                   |
 | comparisons      |  | 436 MTS-domain    |  |                   |
 +------------------+  +-------------------+  |                   |
          |                    |               |                   |
          +----+---------------+               |                   |
               |                               |                   |
               v                               |                   |
     +----------------------------+            |                   |
     |       Module H             |<-----------+                   |
     | Comparative statistics     |                                |
     | 281 tests, 22 significant  |                                |
     +----------------------------+                                |
               |                                                   |
               v                                                   |
     +----------------------------+                                |
     |       Module I             |<-------------------------------+
     | Generate 6 figures         |   (reads all upstream results)
     | PDF + PNG at 300 DPI       |
     +----------------------------+
               |
               v
         FINAL OUTPUTS:
         - 6 publication figures
         - Statistical reports
         - Homolog table
         - Domain/targeting annotations
```

---

## 4. Statistical Analysis Details

### 4.1 Pre-registered Hypotheses

All hypotheses are documented in `docs/PRIMARY_HYPOTHESES.md`. Nine primary hypotheses organized into three goal families:

| ID | Hypothesis | Goal Family |
|----|-----------|-------------|
| H1.1 | GroEL substrates enriched for specific CATH superfamilies vs size-matched cytoplasmic background | Domain architecture |
| H1.2 | HSP60 substrates enriched for specific fold architectures vs size-matched matrix background | Domain architecture |
| H1.3 | Structural fold enrichment conserved between GroEL and HSP60 systems | Domain architecture |
| H2.1 | N-terminal domains have different contact order than C-terminal regions (within-protein paired test) | Stability asymmetry |
| H2.2 | N-vs-C asymmetry greater in chaperonin substrates than in matched background | Stability asymmetry |
| H2.3 | Class III GroEL substrates show greater N-vs-C asymmetry than Class I | Stability asymmetry |
| H3.1 | HSP60 substrates enriched for matrix localization vs general mito proteome | Matrix targeting |
| H3.2 | MTS-bearing HSP60 substrates have distinct first-domain properties | Matrix targeting |
| H3.3 | MTS is predominantly a pre-domain extension, not part of the first domain | Matrix targeting |

### 4.2 Statistical Framework

- **Alpha level:** 0.05 (two-sided)
- **Multiple testing correction:** Hierarchical Benjamini-Hochberg (BH)
  - **Level 1:** Three goal families tested simultaneously with BH
  - **Level 2:** Within each significant family, individual tests corrected with BH
- **Controls:** Compartment-matched AND size-matched (10 kDa bins)

### 4.3 Statistical Tests Used

#### Goal 1: Domain Architecture Enrichment

| Test | Scipy Function | Application | Parameters |
|------|---------------|-------------|------------|
| Fisher's exact test | `scipy.stats.fisher_exact(table, alternative='two-sided')` | Per-superfamily enrichment (GroEL vs background, HSP60 vs background) | 2x2 contingency table |
| Chi-squared test | `scipy.stats.chi2_contingency(table)` | CATH class distribution comparison | k x 2 contingency table |
| Hypergeometric test | `scipy.stats.hypergeom.sf(...)` | Shared superfamily overlap significance | |
| Odds ratio 95% CI | Manual: `log(OR) +/- 1.96 * SE` | Effect size for Fisher tests | Haldane correction (+0.5) for zero cells |

**Specifics from results:**
- H1.1: 123 superfamilies tested, 5 significant after BH correction
  - Top: 1.10.10.10 (Winged helix) OR=42.80, p_BH=2.35e-06
  - Top: 3.20.20.70 (Aldolase class I) OR=9.16, p_BH=2.35e-06
  - CATH class chi-squared: chi2=24.24, dof=4, p=7.16e-05, Cramer's V=0.120
- H1.2: 119 superfamilies tested, 1 significant after BH correction
  - 1.50.40.10 (Mitochondrial carrier domain) OR=0.16, p_BH=2.79e-02
  - CATH class chi-squared: chi2=16.79, dof=4, p=2.13e-03, Cramer's V=0.101
- H1.3: 85 shared superfamilies, Jaccard index=0.202, 55/69 homolog pairs share top SF

#### Goal 2: N-vs-C Stability Asymmetry

| Test | Scipy Function | Application | Parameters |
|------|---------------|-------------|------------|
| Wilcoxon signed-rank | `scipy.stats.wilcoxon(n_values, c_values, alternative='two-sided')` | Within-protein paired N-vs-C comparison | Paired observations |
| Mann-Whitney U | `scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')` | Cross-dataset asymmetry comparison | Independent samples |
| Kruskal-Wallis H | `scipy.stats.kruskal(*groups)` | GroEL class comparison (I vs II vs III) | k independent groups |
| Hotelling's T-squared | Custom implementation: T2 = n * mean' * inv(S) * mean | Multivariate omnibus test across all metrics | F-approximation reported |

**Effect size measures:**
- Wilcoxon: rank-biserial correlation `r = 1 - (2W) / (n*(n+1)/2)` where W is the test statistic
- Mann-Whitney U: rank-biserial `r = 1 - (2U) / (n1*n2)`
- Kruskal-Wallis: eta-squared `eta2 = (H - k + 1) / (N - k)`
- Chi-squared: Cramer's V = `sqrt(chi2 / (N * min(r-1, c-1)))`

**Metrics tested (5 per dataset, within-protein):**
1. `relative_contact_order` (folding kinetics proxy)
2. `mean_plddt` (AlphaFold confidence)
3. `frac_helix` (helix content)
4. `frac_strand` (strand content)
5. `frac_plddt_gt70` (fraction well-modeled)

**Key results from within-protein paired tests (Wilcoxon signed-rank):**

| Dataset | Metric | Direction | W | p_BH | r |
|---------|--------|-----------|---|------|---|
| GroEL | relative_contact_order | N > C | 2262.0 | 8.20e-04 | 0.387 |
| GroEL | mean_plddt | N > C | 3750.0 | 6.13e-03 | 0.282 |
| GroEL | frac_plddt>70 | N > C | 1623.0 | 4.16e-04 | 0.438 |
| HSP60 | relative_contact_order | N > C | 1776.0 | 7.34e-06 | 0.519 |
| Matrix BG | relative_contact_order | N > C | 3195.0 | 3.59e-04 | 0.388 |
| Matrix BG | mean_plddt | N > C | 4815.0 | 8.56e-04 | 0.322 |
| Matrix BG | frac_plddt>70 | N > C | 3543.0 | 8.20e-04 | 0.349 |
| Mito BG | relative_contact_order | N > C | 1086.0 | 9.00e-08 | 0.644 |
| Mito BG | mean_plddt | N > C | 2268.0 | 1.70e-03 | 0.354 |
| Mito BG | frac_strand | N > C | 1667.0 | 1.26e-03 | 0.389 |
| Mito BG | frac_plddt>70 | N > C | 1900.5 | 9.53e-04 | 0.389 |

**Multivariate omnibus test (Hotelling's T-squared):**
| Dataset | T2 | F | p |
|---------|----|----|---|
| GroEL | 18.08 | F(4,117)=4.41 | 2.36e-03 |
| HSP60 | 27.92 | F(4,117)=6.80 | 5.84e-05 |
| Matrix BG | 31.10 | F(4,140)=7.61 | 1.41e-05 |
| Mito BG | 60.19 | F(4,106)=14.63 | 1.51e-09 |

#### Goal 3: Matrix Targeting and MTS

| Test | Scipy Function | Application | Parameters |
|------|---------------|-------------|------------|
| Fisher's exact test | `scipy.stats.fisher_exact(table)` | Matrix enrichment (HSP60 vs non-HSP60) | 2x2 table |
| Fisher's exact test | `scipy.stats.fisher_exact(table)` | MTS prevalence comparison | 2x2 table |
| Binomial test | `scipy.stats.binom_test(k, n, p=0.5)` | MTS pre-domain vs overlap position | H0: p=0.5 |

**Key results:**
- H3.1 Matrix enrichment: OR=3.29 [2.46, 4.40], p=1.60e-16 (181/266 HSP60 matrix vs 343/873 non-HSP60 matrix)
- H3.2 MTS prevalence: OR=1.54 [1.05, 2.25], p=2.95e-02
- H3.3 MTS pre-domain: 368/436 (84.4%), binomial p=3.42e-51

### 4.4 Multiple Testing Correction Details

**Implementation:** `statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='fdr_bh')`

**Hierarchical BH procedure:**
1. Compute raw p-values for all 281 tests.
2. Group tests into three families: domain_architecture, stability_asymmetry, matrix_targeting.
3. Take the minimum raw p-value per family.
4. Apply BH correction to the 3 family-level p-values.
5. If a family is significant at level 0.05, apply BH correction to all tests within that family.
6. Tests in non-significant families are not examined further.

**Family-level results:**
| Family | Min raw p | Family BH p | Significant? |
|--------|-----------|-------------|-------------|
| domain_architecture | 2.82e-08 | 2.82e-08 | Yes |
| matrix_targeting | 3.42e-51 | 1.02e-50 | Yes |
| stability_asymmetry | 4.50e-09 | 6.75e-09 | Yes |

**Summary:** 281 total tests, 22 significant after hierarchical correction.

### 4.5 Fisher's Exact Test Implementation Detail

From `module_h_comparative_stats.py`:
```python
def safe_fisher(a, b, c, d):
    """Fisher's exact test with odds ratio and 95% CI."""
    table = np.array([[a, b], [c, d]])
    odds_ratio, pvalue = stats.fisher_exact(table, alternative='two-sided')
    # Compute 95% CI using log(OR) +/- 1.96*SE
    if a == 0 or b == 0 or c == 0 or d == 0:
        # Haldane correction
        a_c, b_c, c_c, d_c = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        a_c, b_c, c_c, d_c = a, b, c, d
    log_or = np.log(a_c * d_c / (b_c * c_c))
    se = np.sqrt(1/a_c + 1/b_c + 1/c_c + 1/d_c)
    ci_lower = np.exp(log_or - 1.96 * se)
    ci_upper = np.exp(log_or + 1.96 * se)
    return odds_ratio, ci_lower, ci_upper, pvalue
```

---

## 5. Figure Generation

### 5.1 Running Figure Generation

```bash
cd /Users/vishalbharti/Downloads/Antah_Asti_Prarambh
python3 workflow/scripts/generate_figures.py
```

All figures are saved to `results/figures/` in both PDF and PNG (300 DPI) formats.

### 5.2 Figures Produced

| Figure | File | Content |
|--------|------|---------|
| Figure 1 | `fig1_domain_architecture` | CATH class distribution comparison across GroEL, HSP60, Matrix, Mito datasets |
| Figure 2 | `fig2_n_vs_c_stability` | N-domain vs C-region stability asymmetry (contact order, pLDDT, SS fractions) |
| Figure 3 | `fig3_groel_class_comparison` | GroEL class (I/II/III) comparison of N-vs-C asymmetry metrics |
| Figure 4 | `fig4_mts_targeting` | MTS targeting classification and MTS-domain boundary relationship |
| Figure 5 | `fig5_orthology` | Orthology network and RBH statistics between GroEL and HSP60 substrates |
| Figure 6 | `fig6_summary` | Integrated summary of all three goals |

### 5.3 Style Specifications

**Matplotlib style:** `seaborn-v0_8-whitegrid`

**Global rcParams:**
```python
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,    # TrueType fonts (journal requirement)
    "ps.fonttype": 42,
})
```

### 5.4 Color Palette

**Primary palette:** Seaborn `Set2` (colorblind-friendly), 8 colors.

**Dataset color mapping:**
```python
CB_PALETTE = sns.color_palette("Set2", 8)

DATASET_COLORS = {
    "GroEL":  CB_PALETTE[0],   # green-teal
    "HSP60":  CB_PALETTE[1],   # orange
    "Mito":   CB_PALETTE[2],   # blue-gray
    "Matrix": CB_PALETTE[3],   # pink
}

CATH_CLASS_COLORS = {
    "Mainly Alpha": CB_PALETTE[4],   # sage green
    "Mainly Beta":  CB_PALETTE[5],   # golden yellow
    "Alpha Beta":   CB_PALETTE[6],   # salmon pink
    "Few SS":       CB_PALETTE[7],   # gray
}

NC_COLORS = {
    "N-domain": CB_PALETTE[0],   # green-teal
    "C-region": CB_PALETTE[1],   # orange
}
```

### 5.5 Output Format Details

- **PNG:** 300 DPI, tight bounding box
- **PDF:** Type 42 (TrueType) fonts for journal compatibility, tight bounding box
- **Backend:** `Agg` (non-interactive, suitable for scripts/servers)

### 5.6 Input Files Required for Figure Generation

All of the following must exist before running:
- `results/domains/cath_protein_summary.tsv`
- `results/domains/cath_domain_assignments.tsv`
- `results/domains/domain_structural_metrics.tsv`
- `results/domains/foldseek_clusters.tsv`
- `results/termini/n_vs_c_paired.tsv`
- `results/termini/contact_order.tsv`
- `results/termini/region_boundaries.tsv`
- `results/termini/sequence_metrics.tsv`
- `results/termini/structure_metrics.tsv`
- `results/mts/combined_targeting.tsv`
- `results/mts/mts_domain_relationship.tsv`
- `results/homology/rbh_groel_hsp60_annotated.tsv`
- `results/homology/substrate_orthogroups.tsv`
- `data/processed/groel_substrates_standardized.tsv`
- `data/processed/hsp60_tier1_substrates.tsv`
- `results/structures/structure_index.tsv`

---

## 6. Quality Control Checklist

Use this checklist to verify that the full pipeline ran correctly. All numbers are taken directly from the generated reports.

### Module A: Data Acquisition & Cleaning

- [ ] 252 GroEL proteins in `data/processed/groel_substrates_standardized.tsv`
- [ ] GroEL class distribution: Class I = 38, Class I or II = 4, Class II = 126, Class III = 84
- [ ] 325 total HSP60 interactors in `data/processed/hsp60_interactome_standardized.tsv`
- [ ] 266 HSP60 Tier 1 proteins in `data/processed/hsp60_tier1_substrates.tsv`
- [ ] 49 Tier 2 proteins, 10 excluded proteins
- [ ] 100% of Tier 1 proteins are MitoCarta-annotated (266/266)
- [ ] NDIC imputation values: coIP1=80.34, coIP2=80.64, coIP3=91.01
- [ ] 1,136 MitoCarta 3.0 proteins in `data/processed/human_mito_proteome.tsv`
- [ ] 525 matrix-localized proteins in `data/processed/human_matrix_proteome.tsv`
- [ ] 274 of 325 HSP60 interactors in MitoCarta 3.0 (84.3%)
- [ ] 70 localization changes documented between MitoCarta 2.0 and 3.0

### Module B: Dataset Construction

- [ ] 252 sequences in `data/processed/groel_substrates.fasta`
- [ ] 266 sequences in `data/processed/hsp60_tier1_substrates.fasta`

### Module C: Orthology / Homology

- [ ] 40 RBH pairs in `results/homology/rbh_groel_hsp60.tsv`
- [ ] 40 rows in `results/homology/rbh_groel_hsp60_annotated.tsv`
- [ ] RBH by class: Class I = 7, Class I or II = 1, Class II = 24, Class III = 8
- [ ] 8 Class III RBH pairs (9.5% of Class III substrates)
- [ ] Median RBH percent identity: 35.8%
- [ ] Median RBH e-value: 3.39e-58
- [ ] 422 orthogroups in `results/homology/orthogroups_ecoli_human.tsv`
- [ ] 143 substrate-containing orthogroups
- [ ] 34 shared orthogroups (GroEL + HSP60 substrates)
- [ ] 62 substrate pairs from orthogroup method
- [ ] 29 additional pairs beyond simple RBH
- [ ] 69 homolog pairs in `data/processed/groel_hsp60_homologs.tsv`

### Module D: Structures

- [ ] 1,382 AlphaFold CIF files in `data/raw/alphafold/pilot/`
- [ ] 1,382 rows in `results/structures/structure_index.tsv`
- [ ] 1,382 DSSP files in `results/structures/dssp/`
- [ ] `results/structures/dssp_summary.tsv` exists and is non-empty
- [ ] `results/structures/dssp_per_residue.tsv` exists and is non-empty

### Module E: Domain Architecture

- [ ] 2,141 domain assignments in `results/domains/cath_domain_assignments.tsv` (2,142 lines incl. header)
- [ ] 1,390 protein summaries in `results/domains/cath_protein_summary.tsv` (1,391 lines incl. header)
- [ ] GroEL: 247/252 proteins with CATH domains (98.0%), 499 total domains
- [ ] HSP60: 241/266 proteins with CATH domains (90.6%), 471 total domains
- [ ] Matrix: 471/524 proteins with CATH domains (89.9%), 945 total domains
- [ ] Mito: 898/1,132 proteins with CATH domains (79.3%), 1,628 total domains
- [ ] 1,151 total proteins with CATH domains (247 + 241 + 471 + 898 unique across overlapping sets, counted per the union)
- [ ] `results/domains/domain_structural_metrics.tsv` exists and is non-empty
- [ ] `results/domains/domain_distribution_report.txt` exists
- [ ] 1,155 Foldseek clusters in `results/domains/foldseek_cluster_summary.txt`
- [ ] 999 singleton clusters (86.5%)
- [ ] 24 shared structural clusters between GroEL and HSP60
- [ ] `results/domains/foldseek_clusters.tsv` exists

### Module F: N-vs-C Stability

- [ ] 567 multi-domain N-vs-C comparisons in `results/termini/n_vs_c_paired.tsv` (568 lines incl. header)
- [ ] `results/termini/contact_order.tsv` exists and is non-empty
- [ ] `results/termini/region_boundaries.tsv` exists and is non-empty
- [ ] `results/termini/sequence_metrics.tsv` exists and is non-empty
- [ ] `results/termini/structure_metrics.tsv` exists and is non-empty

### Module G: MTS Targeting

- [ ] 1,139 proteins in `results/mts/combined_targeting.tsv`
- [ ] 436 proteins in `results/mts/mts_domain_relationship.tsv`
- [ ] 494 proteins with transit peptide (43.4% of all analyzed)
- [ ] 168 HSP60 substrates with transit peptide (63.2%)
- [ ] 181 HSP60 matrix substrates total
- [ ] 368 MTS pre-domain (84.4%), 68 MTS overlapping (15.6%)
- [ ] Median gap (MTS end to first domain): 18.0 residues
- [ ] `results/mts/targeting_summary_report.txt` exists

### Module H: Statistics

- [ ] 281 total statistical tests in `results/stats/corrected_pvalues.tsv`
- [ ] 22 tests significant after hierarchical BH correction
- [ ] 3 goal families all significant at family level
- [ ] Domain architecture: 8 significant tests (5 SF enrichment + 2 chi-squared + 1 Jaccard)
- [ ] Stability asymmetry: 11 significant tests
- [ ] Matrix targeting: 3 significant tests
- [ ] `results/stats/domain_enrichment.tsv` exists
- [ ] `results/stats/stability_comparisons.tsv` exists
- [ ] `results/stats/targeting_stats.tsv` exists
- [ ] `results/stats/statistics_summary_report.txt` exists

### Module I: Figures

- [ ] `results/figures/fig1_domain_architecture.pdf` and `.png` exist
- [ ] `results/figures/fig2_n_vs_c_stability.pdf` and `.png` exist
- [ ] `results/figures/fig3_groel_class_comparison.pdf` and `.png` exist
- [ ] `results/figures/fig4_mts_targeting.pdf` and `.png` exist
- [ ] `results/figures/fig5_orthology.pdf` and `.png` exist
- [ ] `results/figures/fig6_summary.pdf` and `.png` exist
- [ ] 6 figures generated (12 files total: 6 PDF + 6 PNG)
- [ ] All PNG files are > 0 bytes
- [ ] All PDF files are > 0 bytes

### Cross-Module Consistency

- [ ] Number of RBH pairs (40) consistent between `rbh_summary_report.txt` and `orthology_summary_report.txt`
- [ ] Number of orthogroups (422) consistent between report and output TSV
- [ ] Number of structures (1,382) consistent between structure index, CIF file count, and DSSP file count
- [ ] Domain assignment count (2,141) consistent between report and TSV line count
- [ ] All 69 homolog pairs have both a GroEL and HSP60 accession
- [ ] Class III substrates total 84 across all reports (RBH: 84, domain dist: 84)

---

## References

1. Kerner MJ et al. (2005) "Proteome-wide analysis of chaperonin-dependent protein folding in Escherichia coli." *Cell* 122(2):209-220.
2. Bruderer R et al. (2020) "Analysis of 1508 plasma samples by capillary-flow data-independent acquisition profiles proteomics of weight loss and maintenance." *Molecular & Cellular Proteomics* 19(2):368-384. (HSP60 interactome supplement)
3. Rath S et al. (2021) "MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations." *Nucleic Acids Research* 49(D1):D1541-D1547. PMID: 33174596.
4. Steinegger M, Soding J (2017) "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets." *Nature Biotechnology* 35:1026-1028.
5. van Kempen M et al. (2024) "Fast and accurate protein structure search with Foldseek." *Nature Biotechnology* 42:243-246.
6. Jumper J et al. (2021) "Highly accurate protein structure prediction with AlphaFold." *Nature* 596:583-589.
7. Kabsch W, Sander C (1983) "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features." *Biopolymers* 22(12):2577-2637.
