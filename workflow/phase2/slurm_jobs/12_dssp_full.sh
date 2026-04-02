#!/bin/bash
#SBATCH --job-name=aap_dssp_full
#SBATCH --output=logs/12_dssp_full_%j.out
#SBATCH --error=logs/12_dssp_full_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00

# =============================================================================
# Antah Asti Prarambh — Step 12: Full-scale DSSP
# Run mkdssp on ALL 25,007 AlphaFold structures (not just pilot 1,382)
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 12: Full-scale DSSP — $(date)"
echo "============================================================"

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

python3 -u << 'PYTHON_SCRIPT'
import os
import sys
import glob
import subprocess
import re
import time
import pandas as pd
from collections import Counter
from pathlib import Path

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"
STRUCT_DIR_ECOLI = f"{PROJECT_DIR}/data/raw/alphafold/full/ecoli"
STRUCT_DIR_HUMAN = f"{PROJECT_DIR}/data/raw/alphafold/full/human"

os.makedirs(f"{RESULTS}/structures/dssp", exist_ok=True)

print("=" * 70)
print("Full-Scale DSSP Secondary Structure Assignment")
print("=" * 70)

# Find all CIF files
cif_files = []
for d in [STRUCT_DIR_ECOLI, STRUCT_DIR_HUMAN]:
    if os.path.exists(d):
        cif_files.extend(glob.glob(f"{d}/*.cif"))
        cif_files.extend(glob.glob(f"{d}/*.cif.gz"))  # in case still gzipped

print(f"Found {len(cif_files)} structure files")

# Also check pilot directory
pilot_dir = f"{PROJECT_DIR}/data/raw/alphafold/pilot"
if os.path.exists(pilot_dir):
    pilot_cifs = glob.glob(f"{pilot_dir}/*.cif")
    print(f"Found {len(pilot_cifs)} pilot structure files")
    cif_files.extend(pilot_cifs)

# Extract accession from filename
def get_accession(path):
    m = re.search(r'AF-([A-Z0-9]+)-F', os.path.basename(path))
    return m.group(1) if m else None

# Check which already have DSSP
existing_dssp = set()
dssp_dir = f"{RESULTS}/structures/dssp"
for f in os.listdir(dssp_dir) if os.path.exists(dssp_dir) else []:
    if f.endswith(".dssp"):
        existing_dssp.add(f.replace(".dssp", ""))

# Also check pilot DSSP
pilot_dssp_dir = f"{PROJECT_DIR}/results/structures/dssp"
if os.path.exists(pilot_dssp_dir):
    for f in os.listdir(pilot_dssp_dir):
        if f.endswith(".dssp"):
            existing_dssp.add(f.replace(".dssp", ""))

print(f"Already have DSSP for {len(existing_dssp)} proteins")

# Build work list
work = []
for cif_path in cif_files:
    acc = get_accession(cif_path)
    if acc and acc not in existing_dssp:
        work.append((acc, cif_path))

# Deduplicate
seen = set()
unique_work = []
for acc, path in work:
    if acc not in seen:
        seen.add(acc)
        unique_work.append((acc, path))
work = unique_work

print(f"Need to process: {len(work)} proteins")

# Find mkdssp
mkdssp = None
for candidate in ["mkdssp", f"{os.environ.get('CONDA_PREFIX', '')}/bin/mkdssp"]:
    try:
        subprocess.run([candidate, "--version"], capture_output=True, timeout=5)
        mkdssp = candidate
        break
    except (FileNotFoundError, subprocess.TimeoutExpired):
        continue

if not mkdssp:
    print("ERROR: mkdssp not found")
    sys.exit(1)
print(f"Using: {mkdssp}")

# Process
success = 0
failed = 0
per_protein_records = []

for i, (acc, cif_path) in enumerate(work):
    if (i + 1) % 500 == 0:
        print(f"  Progress: {i+1}/{len(work)} ({100*(i+1)/len(work):.1f}%)")

    dssp_out = f"{dssp_dir}/{acc}.dssp"
    try:
        result = subprocess.run(
            [mkdssp, cif_path, dssp_out],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0 and os.path.exists(dssp_out):
            # Parse DSSP output
            ss_counts = Counter()
            total_residues = 0
            with open(dssp_out) as f:
                in_data = False
                for line in f:
                    if line.startswith("  #  RESIDUE"):
                        in_data = True
                        continue
                    if in_data and len(line) > 16:
                        ss = line[16] if len(line) > 16 else ' '
                        if ss in ('H', 'G', 'I'):
                            ss_counts['helix'] += 1
                        elif ss in ('E', 'B'):
                            ss_counts['strand'] += 1
                        else:
                            ss_counts['coil'] += 1
                        total_residues += 1

            if total_residues > 0:
                per_protein_records.append({
                    "accession": acc,
                    "n_residues": total_residues,
                    "frac_helix": ss_counts['helix'] / total_residues,
                    "frac_strand": ss_counts['strand'] / total_residues,
                    "frac_coil": ss_counts['coil'] / total_residues,
                })
            success += 1
        else:
            failed += 1
    except (subprocess.TimeoutExpired, Exception) as e:
        failed += 1

print(f"\nCompleted: {success} success, {failed} failed")

# Save summary
if per_protein_records:
    df = pd.DataFrame(per_protein_records)

    # Also load pilot DSSP summary and merge
    pilot_dssp_path = f"{PROJECT_DIR}/results/structures/dssp_summary.tsv"
    if os.path.exists(pilot_dssp_path):
        pilot = pd.read_csv(pilot_dssp_path, sep="\t")
        # Harmonize columns
        acc_col = "uniprot_accession" if "uniprot_accession" in pilot.columns else "accession"
        if acc_col != "accession":
            pilot = pilot.rename(columns={acc_col: "accession"})
        # Keep columns that match
        common_cols = [c for c in df.columns if c in pilot.columns]
        pilot_slim = pilot[common_cols]
        df = pd.concat([pilot_slim, df], ignore_index=True)
        df = df.drop_duplicates(subset="accession", keep="last")
        print(f"Merged with pilot: {len(df)} total proteins")

    out_path = f"{RESULTS}/structures/dssp_summary_full.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path} ({len(df)} proteins)")

print("\n" + "=" * 70)
print("Full-scale DSSP complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "DSSP complete — $(date)"
echo "============================================================"
