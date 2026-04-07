#!/bin/bash
#SBATCH --job-name=aap_dssp_cont
#SBATCH --output=logs/12_dssp_full_continue_%j.out
#SBATCH --error=logs/12_dssp_full_continue_%j.err
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00

# =============================================================================
# Antah Asti Prarambh — Step 12 CONTINUATION: Full-scale DSSP
# Processes remaining ~1,143 proteins, then generates summary TSV from ALL .dssp files
# =============================================================================

set -euo pipefail

export PROJECT_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc"
CONDA_ENV="proteomics"
RESULTS="${PROJECT_DIR}/results/phase2"

echo "============================================================"
echo "  Step 12 CONTINUATION: Full-scale DSSP — $(date)"
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
import pandas as pd
from collections import Counter
from pathlib import Path

PROJECT_DIR = os.environ["PROJECT_DIR"]
RESULTS = f"{PROJECT_DIR}/results/phase2"
STRUCT_DIR_ECOLI = f"{PROJECT_DIR}/data/raw/alphafold/full/ecoli"
STRUCT_DIR_HUMAN = f"{PROJECT_DIR}/data/raw/alphafold/full/human"
dssp_dir = f"{RESULTS}/structures/dssp"
pilot_dssp_dir = f"{PROJECT_DIR}/results/structures/dssp"

os.makedirs(dssp_dir, exist_ok=True)

print("=" * 70)
print("PHASE 1: Process remaining proteins that don't have .dssp files yet")
print("=" * 70)

# Find all CIF files
cif_files = []
for d in [STRUCT_DIR_ECOLI, STRUCT_DIR_HUMAN]:
    if os.path.exists(d):
        cif_files.extend(glob.glob(f"{d}/*.cif"))

pilot_dir = f"{PROJECT_DIR}/data/raw/alphafold/pilot"
if os.path.exists(pilot_dir):
    cif_files.extend(glob.glob(f"{pilot_dir}/*.cif"))

print(f"Found {len(cif_files)} total structure files")

def get_accession(path):
    m = re.search(r'AF-([A-Z0-9]+)-F', os.path.basename(path))
    return m.group(1) if m else None

# Check which already have DSSP (from both phase2 dssp dir AND pilot dssp dir)
existing_dssp = set()
for d in [dssp_dir, pilot_dssp_dir]:
    if os.path.exists(d):
        for f in os.listdir(d):
            if f.endswith(".dssp"):
                existing_dssp.add(f.replace(".dssp", ""))

print(f"Already have DSSP for {len(existing_dssp)} proteins")

# Build work list of remaining proteins
work = []
seen = set()
for cif_path in cif_files:
    acc = get_accession(cif_path)
    if acc and acc not in existing_dssp and acc not in seen:
        seen.add(acc)
        work.append((acc, cif_path))

print(f"Need to process: {len(work)} remaining proteins")

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

# Process remaining proteins
success = 0
failed = 0
for i, (acc, cif_path) in enumerate(work):
    if (i + 1) % 100 == 0:
        print(f"  Progress: {i+1}/{len(work)} ({100*(i+1)/len(work):.1f}%)")

    dssp_out = f"{dssp_dir}/{acc}.dssp"
    try:
        result = subprocess.run(
            [mkdssp, cif_path, dssp_out],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0 and os.path.exists(dssp_out):
            success += 1
        else:
            failed += 1
    except (subprocess.TimeoutExpired, Exception):
        failed += 1

print(f"\nRemaining proteins: {success} success, {failed} failed")

# =========================================================================
print("\n" + "=" * 70)
print("PHASE 2: Generate summary TSV from ALL .dssp files")
print("=" * 70)

def parse_dssp_file(dssp_path):
    """Parse a single .dssp file and return SS fractions."""
    ss_counts = Counter()
    total_residues = 0
    try:
        with open(dssp_path) as f:
            in_data = False
            for line in f:
                if line.startswith("  #  RESIDUE"):
                    in_data = True
                    continue
                if in_data and len(line) > 16:
                    ss = line[16]
                    if ss in ('H', 'G', 'I'):
                        ss_counts['helix'] += 1
                    elif ss in ('E', 'B'):
                        ss_counts['strand'] += 1
                    else:
                        ss_counts['coil'] += 1
                    total_residues += 1
    except Exception:
        return None

    if total_residues == 0:
        return None

    return {
        "n_residues": total_residues,
        "frac_helix": ss_counts['helix'] / total_residues,
        "frac_strand": ss_counts['strand'] / total_residues,
        "frac_coil": ss_counts['coil'] / total_residues,
    }

# Collect ALL .dssp files from both phase2 and pilot directories
all_dssp_files = {}

# Phase 2 .dssp files (takes priority)
for f in sorted(os.listdir(dssp_dir)):
    if f.endswith(".dssp"):
        acc = f.replace(".dssp", "")
        all_dssp_files[acc] = os.path.join(dssp_dir, f)

# Pilot .dssp files (only if not already in phase2)
if os.path.exists(pilot_dssp_dir):
    for f in sorted(os.listdir(pilot_dssp_dir)):
        if f.endswith(".dssp"):
            acc = f.replace(".dssp", "")
            if acc not in all_dssp_files:
                all_dssp_files[acc] = os.path.join(pilot_dssp_dir, f)

print(f"Total .dssp files to parse: {len(all_dssp_files)}")

# Parse all files
records = []
parse_failed = 0
for i, (acc, path) in enumerate(all_dssp_files.items()):
    if (i + 1) % 2000 == 0:
        print(f"  Parsing: {i+1}/{len(all_dssp_files)} ({100*(i+1)/len(all_dssp_files):.1f}%)")

    result = parse_dssp_file(path)
    if result:
        result["accession"] = acc
        records.append(result)
    else:
        parse_failed += 1

print(f"Parsed: {len(records)} proteins, {parse_failed} failed")

# Save summary TSV
if records:
    df = pd.DataFrame(records)
    # Reorder columns
    df = df[["accession", "n_residues", "frac_helix", "frac_strand", "frac_coil"]]
    df = df.sort_values("accession").reset_index(drop=True)

    out_path = f"{RESULTS}/structures/dssp_summary_full.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved: {out_path}")
    print(f"Total proteins in summary: {len(df)}")
    print(f"  Mean helix: {df['frac_helix'].mean():.3f}")
    print(f"  Mean strand: {df['frac_strand'].mean():.3f}")
    print(f"  Mean coil: {df['frac_coil'].mean():.3f}")
else:
    print("ERROR: No records to save!")
    sys.exit(1)

print("\n" + "=" * 70)
print("Full-scale DSSP continuation complete.")
print("=" * 70)
PYTHON_SCRIPT

echo ""
echo "DSSP continuation complete — $(date)"
echo "============================================================"
