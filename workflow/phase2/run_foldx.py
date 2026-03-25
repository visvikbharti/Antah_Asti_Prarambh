#!/usr/bin/env python3
"""
Phase 2: FoldX stability calculations on full proteome AlphaFold structures.

Antah Asti Prarambh — Module F (Full-Scale)

Computes thermodynamic stability estimates (DeltaG) for all AlphaFold structures
using FoldX RepairPDB + Stability commands. Supports SLURM array job
parallelization for processing ~25,000 proteins.

Pipeline per protein:
  1. Convert CIF -> PDB (FoldX requires PDB format)
  2. RepairPDB: optimize side chain rotamers and fix clashes
  3. Stability: compute folding free energy (DeltaG)
  4. (Optional) AnalyseComplex: domain-domain interface energy

Estimated resources per protein:
  - RepairPDB: ~30-120s, 1-2 GB RAM
  - Stability:  ~10-30s, <1 GB RAM
  - Total for 25,000 proteins at 50/job: ~500 SLURM array tasks
  - Walltime per task: ~4 hours

Usage:
  # Full run (serial, for small batches)
  python run_foldx.py --config config.yaml --organism ecoli

  # SLURM array mode: generate job scripts
  python run_foldx.py --config config.yaml --generate-slurm

  # Run a specific chunk (called by SLURM array)
  python run_foldx.py --config config.yaml --chunk-id 0 --chunk-size 50

  # Collect results from all chunks
  python run_foldx.py --config config.yaml --collect

SLURM array submission:
  # After --generate-slurm:
  sbatch foldx_array.slurm
"""

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml


# =============================================================================
# CONFIGURATION
# =============================================================================

def load_config(config_path):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    project_dir = os.path.expandvars(cfg["project_dir"])
    phase1_dir = os.path.expandvars(cfg.get("phase1_dir", project_dir + "/phase1"))

    def resolve(d):
        if isinstance(d, dict):
            return {k: resolve(v) for k, v in d.items()}
        elif isinstance(d, str):
            return d.replace("{project_dir}", project_dir).replace("{phase1_dir}", phase1_dir)
        elif isinstance(d, list):
            return [resolve(v) for v in d]
        return d

    return resolve(cfg)


# =============================================================================
# CIF TO PDB CONVERSION
# =============================================================================

def cif_to_pdb(cif_path, pdb_path):
    """
    Convert mmCIF to PDB format using gemmi or BioPython.
    FoldX requires PDB format input.
    """
    # Try gemmi first (fast, reliable)
    try:
        import gemmi
        structure = gemmi.read_structure(cif_path)
        structure.write_pdb(pdb_path)
        return True
    except ImportError:
        pass

    # Try BioPython
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("protein", cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        return True
    except ImportError:
        pass

    # Try pdbx/mmcif -> PDB with command line tools
    try:
        result = subprocess.run(
            ["gemmi", "convert", cif_path, pdb_path],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode == 0:
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    print(f"    WARNING: Could not convert {cif_path} to PDB. "
          f"Install gemmi: pip install gemmi")
    return False


# =============================================================================
# FOLDX OPERATIONS
# =============================================================================

def run_foldx_repair(foldx_bin, pdb_path, work_dir, rotabase=None, timeout=300):
    """
    Run FoldX RepairPDB to optimize side chain rotamers and fix steric clashes.
    This should always be run before Stability calculations on AlphaFold models.

    Returns path to repaired PDB file, or None on failure.
    """
    pdb_name = os.path.basename(pdb_path)
    repaired_name = pdb_name.replace(".pdb", "_Repair.pdb")
    repaired_path = os.path.join(work_dir, repaired_name)

    # Skip if already repaired
    if os.path.exists(repaired_path) and os.path.getsize(repaired_path) > 100:
        return repaired_path

    # Copy PDB to working directory
    work_pdb = os.path.join(work_dir, pdb_name)
    if not os.path.exists(work_pdb):
        shutil.copy2(pdb_path, work_pdb)

    # Build FoldX command
    cmd = [
        foldx_bin,
        "--command=RepairPDB",
        f"--pdb={pdb_name}",
        f"--output-dir={work_dir}",
    ]

    if rotabase and os.path.exists(rotabase):
        cmd.append(f"--rotabaseLocation={rotabase}")

    env = os.environ.copy()
    # Some FoldX versions need this
    if rotabase:
        env["FOLDX_ROTABASE"] = rotabase

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=work_dir,
            timeout=timeout,
            env=env,
        )

        if result.returncode != 0:
            # Check if repair file was still created (FoldX sometimes exits non-zero)
            if os.path.exists(repaired_path) and os.path.getsize(repaired_path) > 100:
                return repaired_path
            return None

        if os.path.exists(repaired_path):
            return repaired_path
        else:
            # FoldX sometimes uses different naming convention
            alt_name = f"RepairPDB_{pdb_name}"
            alt_path = os.path.join(work_dir, alt_name)
            if os.path.exists(alt_path):
                os.rename(alt_path, repaired_path)
                return repaired_path
            return None

    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT: RepairPDB for {pdb_name} ({timeout}s)")
        return None
    except Exception as e:
        print(f"    ERROR: RepairPDB for {pdb_name}: {e}")
        return None


def run_foldx_stability(foldx_bin, pdb_path, work_dir, rotabase=None,
                        temperature=298.15, ph=7.0, ionic_strength=0.05,
                        timeout=120):
    """
    Run FoldX Stability command to compute folding free energy (DeltaG).

    Returns dict with stability metrics, or None on failure.
    FoldX Stability output includes:
      - total energy (DeltaG)
      - backbone H-bond energy
      - sidechain-sidechain H-bond energy
      - van der Waals clashes
      - electrostatics
      - solvation (polar + hydrophobic)
      - entropy (mainchain + sidechain)
    """
    pdb_name = os.path.basename(pdb_path)

    # Copy PDB to working directory if needed
    work_pdb = os.path.join(work_dir, pdb_name)
    if not os.path.exists(work_pdb):
        shutil.copy2(pdb_path, work_pdb)

    cmd = [
        foldx_bin,
        "--command=Stability",
        f"--pdb={pdb_name}",
        f"--output-dir={work_dir}",
        f"--temperature={temperature}",
        f"--pH={ph}",
        f"--ionStrength={ionic_strength}",
    ]

    if rotabase and os.path.exists(rotabase):
        cmd.append(f"--rotabaseLocation={rotabase}")

    env = os.environ.copy()
    if rotabase:
        env["FOLDX_ROTABASE"] = rotabase

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=work_dir,
            timeout=timeout,
            env=env,
        )

        # Parse Stability output
        # FoldX writes to work_dir/<pdb_name>_0_ST.fxout or similar
        stability_data = parse_foldx_stability_output(work_dir, pdb_name)

        if stability_data:
            return stability_data
        else:
            # Try parsing from stdout as fallback
            return parse_foldx_stdout(result.stdout, pdb_name)

    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT: Stability for {pdb_name} ({timeout}s)")
        return None
    except Exception as e:
        print(f"    ERROR: Stability for {pdb_name}: {e}")
        return None


def parse_foldx_stability_output(work_dir, pdb_name):
    """
    Parse FoldX Stability output files.
    Looks for *_ST.fxout or Stability_*.fxout files.
    """
    # Search for output files
    stem = pdb_name.replace(".pdb", "")
    patterns = [
        os.path.join(work_dir, f"{stem}_0_ST.fxout"),
        os.path.join(work_dir, f"Stability_{stem}.fxout"),
        os.path.join(work_dir, f"*_ST.fxout"),
    ]

    for pattern in patterns:
        matches = glob.glob(pattern)
        for fpath in matches:
            try:
                with open(fpath) as f:
                    lines = f.readlines()

                # FoldX stability output format:
                # Header line with field names, then data line(s)
                for i, line in enumerate(lines):
                    if line.strip().startswith(stem) or line.strip().startswith(pdb_name):
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            return {
                                "total_energy": safe_float(parts[1]) if len(parts) > 1 else None,
                                "backbone_hbond": safe_float(parts[2]) if len(parts) > 2 else None,
                                "sidechain_hbond": safe_float(parts[3]) if len(parts) > 3 else None,
                                "vdw_clashes": safe_float(parts[4]) if len(parts) > 4 else None,
                                "electrostatics": safe_float(parts[5]) if len(parts) > 5 else None,
                                "solvation_polar": safe_float(parts[6]) if len(parts) > 6 else None,
                                "solvation_hydrophobic": safe_float(parts[7]) if len(parts) > 7 else None,
                                "entropy_mainchain": safe_float(parts[8]) if len(parts) > 8 else None,
                                "entropy_sidechain": safe_float(parts[9]) if len(parts) > 9 else None,
                            }
            except Exception:
                continue

    return None


def parse_foldx_stdout(stdout, pdb_name):
    """Fallback: parse stability values from FoldX stdout."""
    if not stdout:
        return None

    result = {}
    for line in stdout.split("\n"):
        if "Total" in line and "=" in line:
            m = re.search(r"Total\s*=\s*([-\d.]+)", line)
            if m:
                result["total_energy"] = float(m.group(1))
        elif "total energy" in line.lower():
            m = re.search(r"([-\d.]+)\s*$", line)
            if m:
                result["total_energy"] = float(m.group(1))

    return result if result else None


def safe_float(s):
    """Convert string to float, returning None on failure."""
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


# =============================================================================
# PROCESS A BATCH OF PROTEINS
# =============================================================================

def process_protein(accession, cif_path, output_dir, cfg):
    """
    Full FoldX pipeline for a single protein:
    1. Convert CIF -> PDB
    2. RepairPDB
    3. Stability calculation
    """
    foldx_cfg = cfg.get("foldx", {})
    foldx_bin = foldx_cfg.get("binary", "foldx")
    rotabase = foldx_cfg.get("rotabase", "")
    temperature = foldx_cfg.get("temperature", 298.15)
    ph = foldx_cfg.get("ph", 7.0)
    ionic_strength = foldx_cfg.get("ionic_strength", 0.05)
    timeout = foldx_cfg.get("timeout_per_protein", 300)

    result = {
        "accession": accession,
        "status": "pending",
        "total_energy": None,
        "backbone_hbond": None,
        "sidechain_hbond": None,
        "vdw_clashes": None,
        "electrostatics": None,
        "solvation_polar": None,
        "solvation_hydrophobic": None,
        "entropy_mainchain": None,
        "entropy_sidechain": None,
        "error": None,
    }

    # Create per-protein working directory
    work_dir = os.path.join(output_dir, "work", accession)
    os.makedirs(work_dir, exist_ok=True)

    # Check if result already exists
    result_file = os.path.join(output_dir, "per_protein", f"{accession}.json")
    if os.path.exists(result_file):
        try:
            with open(result_file) as f:
                cached = json.load(f)
            if cached.get("status") == "success":
                return cached
        except Exception:
            pass

    try:
        # Step 1: CIF -> PDB
        pdb_path = os.path.join(work_dir, f"{accession}.pdb")
        if not os.path.exists(pdb_path):
            ok = cif_to_pdb(cif_path, pdb_path)
            if not ok:
                result["status"] = "failed"
                result["error"] = "CIF_to_PDB_conversion_failed"
                save_result(result, result_file)
                return result

        # Step 2: RepairPDB
        if foldx_cfg.get("repair_pdb", True):
            repaired = run_foldx_repair(
                foldx_bin, pdb_path, work_dir, rotabase, timeout
            )
            if repaired:
                pdb_for_stability = repaired
            else:
                # Use unrepaired PDB as fallback
                pdb_for_stability = pdb_path
                result["error"] = "RepairPDB_failed_using_unrepaired"
        else:
            pdb_for_stability = pdb_path

        # Step 3: Stability
        if foldx_cfg.get("stability", True):
            stability = run_foldx_stability(
                foldx_bin, pdb_for_stability, work_dir, rotabase,
                temperature, ph, ionic_strength, timeout,
            )
            if stability:
                result.update(stability)
                result["status"] = "success"
            else:
                result["status"] = "failed"
                result["error"] = result.get("error", "") + "|Stability_failed"
        else:
            result["status"] = "skipped"

    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)[:200]

    save_result(result, result_file)
    return result


def save_result(result, path):
    """Save per-protein result as JSON."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(result, f, indent=2)


# =============================================================================
# CHUNK MANAGEMENT (FOR SLURM ARRAY JOBS)
# =============================================================================

def get_protein_list(cfg, organism="both"):
    """
    Build list of (accession, cif_path) for all proteins to process.
    """
    proteins = []
    cif_pattern = re.compile(r"AF-([A-Z0-9]+)-F1-model_v\d+\.cif(\.gz)?$")

    dirs = []
    if organism in ("ecoli", "both"):
        dirs.append(cfg["data"]["alphafold_ecoli_dir"])
    if organism in ("human", "both"):
        dirs.append(cfg["data"]["alphafold_human_dir"])

    for struct_dir in dirs:
        if not os.path.isdir(struct_dir):
            print(f"  WARNING: Structure directory not found: {struct_dir}")
            continue

        for fname in sorted(os.listdir(struct_dir)):
            m = cif_pattern.match(fname)
            if m:
                acc = m.group(1)
                fpath = os.path.join(struct_dir, fname)
                if os.path.getsize(fpath) > 100:
                    proteins.append((acc, fpath))

    print(f"  Found {len(proteins)} proteins to process")
    return proteins


def process_chunk(proteins, chunk_id, chunk_size, output_dir, cfg):
    """
    Process a specific chunk of proteins (for SLURM array jobs).
    chunk_id: 0-based index
    """
    start = chunk_id * chunk_size
    end = min(start + chunk_size, len(proteins))

    if start >= len(proteins):
        print(f"  Chunk {chunk_id} is out of range (only {len(proteins)} proteins)")
        return []

    chunk = proteins[start:end]
    print(f"\n  Processing chunk {chunk_id}: proteins {start}-{end-1} "
          f"({len(chunk)} proteins)")

    results = []
    t0 = time.time()

    for i, (acc, cif_path) in enumerate(chunk, 1):
        result = process_protein(acc, cif_path, output_dir, cfg)
        results.append(result)

        if i % 10 == 0 or i == len(chunk):
            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else 0
            n_ok = sum(1 for r in results if r.get("status") == "success")
            n_fail = sum(1 for r in results if r.get("status") in ("failed", "error"))
            print(f"    [{i}/{len(chunk)}] ok={n_ok} fail={n_fail} "
                  f"rate={rate:.1f}/min elapsed={elapsed:.0f}s")

    # Save chunk summary
    chunk_summary_path = os.path.join(output_dir, "chunks", f"chunk_{chunk_id:04d}.json")
    os.makedirs(os.path.dirname(chunk_summary_path), exist_ok=True)
    with open(chunk_summary_path, "w") as f:
        json.dump({
            "chunk_id": chunk_id,
            "start": start,
            "end": end,
            "n_proteins": len(chunk),
            "n_success": sum(1 for r in results if r.get("status") == "success"),
            "n_failed": sum(1 for r in results if r.get("status") in ("failed", "error")),
            "elapsed_sec": time.time() - t0,
        }, f, indent=2)

    return results


# =============================================================================
# COLLECT AND MERGE RESULTS
# =============================================================================

def collect_results(output_dir, proteins):
    """
    Collect all per-protein JSON results into a single TSV file.
    """
    print(f"\n  Collecting results from {output_dir}...")

    per_protein_dir = os.path.join(output_dir, "per_protein")
    if not os.path.isdir(per_protein_dir):
        print(f"  WARNING: No per-protein results found in {per_protein_dir}")
        return

    fields = [
        "accession", "status", "total_energy",
        "backbone_hbond", "sidechain_hbond",
        "vdw_clashes", "electrostatics",
        "solvation_polar", "solvation_hydrophobic",
        "entropy_mainchain", "entropy_sidechain",
        "error",
    ]

    results = []
    for fname in sorted(os.listdir(per_protein_dir)):
        if not fname.endswith(".json"):
            continue
        fpath = os.path.join(per_protein_dir, fname)
        try:
            with open(fpath) as f:
                data = json.load(f)
            results.append(data)
        except Exception as e:
            print(f"    WARNING: Could not parse {fname}: {e}")

    # Write TSV
    out_tsv = os.path.join(output_dir, "foldx_stability_all.tsv")
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_success = sum(1 for r in results if r.get("status") == "success")
    n_failed = sum(1 for r in results if r.get("status") in ("failed", "error"))
    n_skipped = sum(1 for r in results if r.get("status") == "skipped")

    print(f"  Results collected: {len(results)} proteins")
    print(f"    Success: {n_success}")
    print(f"    Failed:  {n_failed}")
    print(f"    Skipped: {n_skipped}")
    print(f"  Output: {out_tsv}")

    # Summary statistics for successful runs
    energies = [r["total_energy"] for r in results
                if r.get("status") == "success" and r.get("total_energy") is not None]
    if energies:
        import statistics
        print(f"\n  DeltaG summary (n={len(energies)}):")
        print(f"    Mean:   {statistics.mean(energies):.2f} kcal/mol")
        print(f"    Median: {statistics.median(energies):.2f} kcal/mol")
        print(f"    Std:    {statistics.stdev(energies):.2f} kcal/mol")
        print(f"    Min:    {min(energies):.2f} kcal/mol")
        print(f"    Max:    {max(energies):.2f} kcal/mol")

    return out_tsv


# =============================================================================
# SLURM ARRAY JOB GENERATION
# =============================================================================

def generate_slurm_script(cfg, proteins, output_dir, config_path):
    """
    Generate a SLURM array job script for parallel FoldX processing.
    """
    foldx_cfg = cfg.get("foldx", {})
    slurm_cfg = cfg.get("slurm", {})

    chunk_size = foldx_cfg.get("proteins_per_job", 50)
    n_chunks = (len(proteins) + chunk_size - 1) // chunk_size

    mem = slurm_cfg.get("foldx_stability", {}).get("mem", "4G")
    wall_time = slurm_cfg.get("foldx_stability", {}).get("time", "04:00:00")
    cpus = slurm_cfg.get("foldx_stability", {}).get("cpus", 1)
    partition = slurm_cfg.get("partition", "normal")
    account = slurm_cfg.get("account", "my_allocation")

    script_path = os.path.join(output_dir, "foldx_array.slurm")

    script = f"""#!/bin/bash
#SBATCH --job-name=aap_foldx
#SBATCH --partition={partition}
#SBATCH --array=0-{n_chunks - 1}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={wall_time}
#SBATCH --output={output_dir}/logs/foldx_%A_%a.out
#SBATCH --error={output_dir}/logs/foldx_%A_%a.err

# ============================================================================
# FoldX Stability Calculations — SLURM Array Job
# Antah Asti Prarambh Phase 2
#
# Total proteins: {len(proteins)}
# Chunk size: {chunk_size}
# Total chunks: {n_chunks}
# ============================================================================

echo "=== FoldX Array Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib:${{LD_LIBRARY_PATH:-}}"

# Run chunk
python3 -u {os.path.abspath(__file__)} \\
    --config {os.path.abspath(config_path)} \\
    --chunk-id $SLURM_ARRAY_TASK_ID \\
    --chunk-size {chunk_size}

echo ""
echo "End: $(date)"
echo "Exit code: $?"
"""

    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)

    # Also generate a collection script
    collect_script_path = os.path.join(output_dir, "collect_results.slurm")
    collect_script = f"""#!/bin/bash
#SBATCH --job-name=aap_foldx_col
#SBATCH --partition={partition}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output={output_dir}/logs/foldx_collect_%j.out
#SBATCH --dependency=afterok:${{1:-0}}

source /home/vishal.bharti/miniconda3/etc/profile.d/conda.sh
conda activate proteomics
export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib:${{LD_LIBRARY_PATH:-}}"

# Run after all array jobs complete
python3 -u {os.path.abspath(__file__)} \\
    --config {os.path.abspath(config_path)} \\
    --collect

echo "Collection complete: $(date)"
"""

    with open(collect_script_path, "w") as f:
        f.write(collect_script)
    os.chmod(collect_script_path, 0o755)

    # Create logs directory
    os.makedirs(os.path.join(output_dir, "logs"), exist_ok=True)

    print(f"\n  SLURM scripts generated:")
    print(f"    Array job:  {script_path}")
    print(f"    Collection: {collect_script_path}")
    print(f"")
    print(f"  To submit:")
    print(f"    JOB_ID=$(sbatch {script_path} | awk '{{print $4}}')")
    print(f"    sbatch --dependency=afterok:$JOB_ID {collect_script_path}")
    print(f"")
    print(f"  Array tasks: {n_chunks}")
    print(f"  Proteins per task: {chunk_size}")
    print(f"  Total proteins: {len(proteins)}")

    return script_path


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="FoldX stability calculations for full proteomes (Phase 2)"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--organism", choices=["ecoli", "human", "both"], default="both",
        help="Which organism(s) to process (default: both)",
    )

    # Execution modes
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--generate-slurm", action="store_true",
        help="Generate SLURM array job script (don't run FoldX)",
    )
    mode.add_argument(
        "--chunk-id", type=int, default=None,
        help="Process a specific chunk (for SLURM array jobs)",
    )
    mode.add_argument(
        "--collect", action="store_true",
        help="Collect results from completed chunks into single TSV",
    )
    mode.add_argument(
        "--serial", action="store_true",
        help="Run all proteins serially (for small batches / testing)",
    )

    parser.add_argument(
        "--chunk-size", type=int, default=None,
        help="Number of proteins per chunk (default: from config)",
    )

    args = parser.parse_args()

    cfg = load_config(args.config)
    foldx_cfg = cfg.get("foldx", {})
    output_dir = cfg["results"]["foldx_dir"]
    os.makedirs(output_dir, exist_ok=True)

    chunk_size = args.chunk_size or foldx_cfg.get("proteins_per_job", 50)

    print("=" * 70)
    print("  PHASE 2: FoldX Stability Calculations")
    print(f"  Organism: {args.organism}")
    print(f"  Output: {output_dir}")
    print("=" * 70)

    # Get protein list
    proteins = get_protein_list(cfg, args.organism)

    if not proteins and not args.collect:
        print("  ERROR: No proteins found. Check AlphaFold download directories.")
        return 1

    if args.generate_slurm:
        generate_slurm_script(cfg, proteins, output_dir, args.config)

    elif args.chunk_id is not None:
        process_chunk(proteins, args.chunk_id, chunk_size, output_dir, cfg)

    elif args.collect:
        collect_results(output_dir, proteins)

    elif args.serial:
        print(f"\n  Running serially on {len(proteins)} proteins...")
        t0 = time.time()
        results = []

        for i, (acc, cif_path) in enumerate(proteins, 1):
            result = process_protein(acc, cif_path, output_dir, cfg)
            results.append(result)

            if i % 50 == 0 or i == len(proteins):
                elapsed = time.time() - t0
                n_ok = sum(1 for r in results if r.get("status") == "success")
                print(f"  [{i}/{len(proteins)}] ok={n_ok} elapsed={elapsed:.0f}s")

        collect_results(output_dir, proteins)

    else:
        print("\n  No execution mode specified. Use one of:")
        print("    --generate-slurm   Generate SLURM array job scripts")
        print("    --serial           Run all proteins serially")
        print("    --chunk-id N       Process chunk N (for SLURM array)")
        print("    --collect          Collect results from completed chunks")
        return 1

    print(f"\n{'='*70}")
    print("  FoldX pipeline step complete")
    print(f"{'='*70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
