#!/usr/bin/env python3
"""
Phase 2: Download full AlphaFold proteome structures for E. coli K-12 and Human.

Antah Asti Prarambh — Module D (Full-Scale)

Downloads bulk tar archives from AlphaFold EBI FTP, with fallback to individual
download for any missing accessions. Includes resume capability, checkpointing,
integrity verification, and progress tracking.

Estimated resources:
  - Disk: ~2 GB (E. coli) + ~20 GB (Human) = ~22 GB total
  - RAM: <4 GB
  - Time: 1-4 hours depending on network speed
  - SLURM: 1 node, 1 CPU, 4 GB RAM, 12 hour walltime

Usage:
  python download_alphafold_full.py --config config.yaml
  python download_alphafold_full.py --config config.yaml --organism ecoli
  python download_alphafold_full.py --config config.yaml --organism human --workers 4

Resume:
  Simply re-run the same command. Already-downloaded files are skipped.
  A checkpoint file tracks progress for bulk tar extraction.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tarfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml

try:
    import requests
except ImportError:
    requests = None  # Will use wget/curl fallback for bulk downloads


# =============================================================================
# CONFIGURATION
# =============================================================================

def load_config(config_path):
    """Load and resolve config.yaml with variable substitution."""
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
# CHECKPOINT MANAGEMENT
# =============================================================================

class CheckpointManager:
    """Track download progress with a JSON checkpoint file."""

    def __init__(self, checkpoint_path):
        self.path = checkpoint_path
        self.state = self._load()

    def _load(self):
        if os.path.exists(self.path):
            with open(self.path) as f:
                return json.load(f)
        return {
            "bulk_downloads": {},
            "extracted": {},
            "individual_downloads": {},
            "verified": {},
            "failed": [],
            "last_updated": None,
        }

    def save(self):
        self.state["last_updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.state, f, indent=2)
        os.replace(tmp, self.path)

    def mark_bulk_downloaded(self, organism, tar_path):
        self.state["bulk_downloads"][organism] = {
            "path": tar_path,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        self.save()

    def is_bulk_downloaded(self, organism):
        return organism in self.state["bulk_downloads"]

    def mark_extracted(self, organism, count):
        self.state["extracted"][organism] = {
            "count": count,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        self.save()

    def is_extracted(self, organism):
        return organism in self.state["extracted"]

    def mark_individual_downloaded(self, accession):
        self.state["individual_downloads"][accession] = True

    def is_individual_downloaded(self, accession):
        return accession in self.state["individual_downloads"]

    def mark_verified(self, organism, stats):
        self.state["verified"][organism] = {
            **stats,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        self.save()

    def add_failed(self, accession, reason):
        self.state["failed"].append({"accession": accession, "reason": reason})


# =============================================================================
# BULK DOWNLOAD (TAR ARCHIVES FROM EBI FTP)
# =============================================================================

def download_bulk_tar(url, dest_path, max_retries=3, retry_delay=5):
    """
    Download bulk AlphaFold tar archive.
    Uses wget/curl for robustness with large files (20+ GB).
    Supports resume via wget -c / curl -C -.
    """
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)

    for attempt in range(1, max_retries + 1):
        print(f"  Download attempt {attempt}/{max_retries}: {url}")
        print(f"  Destination: {dest_path}")

        # Try wget first (best for large files with resume)
        try:
            cmd = [
                "wget", "-c",  # continue/resume
                "--progress=dot:giga",
                "-O", dest_path,
                "--timeout=120",
                "--tries=3",
                url,
            ]
            result = subprocess.run(cmd, capture_output=False, timeout=43200)  # 12h timeout
            if result.returncode == 0 and os.path.exists(dest_path):
                size_gb = os.path.getsize(dest_path) / (1024 ** 3)
                print(f"  Download complete: {size_gb:.2f} GB")
                return True
        except FileNotFoundError:
            pass  # wget not available, try curl
        except subprocess.TimeoutExpired:
            print(f"  WARNING: wget timed out on attempt {attempt}")
            continue

        # Fallback to curl
        try:
            cmd = [
                "curl", "-C", "-",  # resume
                "-L",  # follow redirects
                "--progress-bar",
                "--retry", "3",
                "--retry-delay", str(retry_delay),
                "-o", dest_path,
                url,
            ]
            result = subprocess.run(cmd, capture_output=False, timeout=43200)
            if result.returncode == 0 and os.path.exists(dest_path):
                size_gb = os.path.getsize(dest_path) / (1024 ** 3)
                print(f"  Download complete: {size_gb:.2f} GB")
                return True
        except FileNotFoundError:
            print("  ERROR: Neither wget nor curl available!")
            return False
        except subprocess.TimeoutExpired:
            print(f"  WARNING: curl timed out on attempt {attempt}")
            continue

        if attempt < max_retries:
            print(f"  Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

    print(f"  FAILED after {max_retries} attempts")
    return False


def extract_tar(tar_path, dest_dir, organism_name):
    """
    Extract CIF/PDB files from AlphaFold bulk tar archive.
    Tracks extraction progress and handles partial extractions on resume.
    """
    os.makedirs(dest_dir, exist_ok=True)

    print(f"\n  Extracting {tar_path} -> {dest_dir}")
    t0 = time.time()
    count = 0
    skipped = 0

    try:
        with tarfile.open(tar_path, "r:*") as tar:
            members = tar.getmembers()
            total = len(members)
            print(f"  Archive contains {total} files")

            for i, member in enumerate(members, 1):
                if not member.isfile():
                    continue

                # Only extract structure files (CIF or PDB)
                basename = os.path.basename(member.name)
                if not (basename.endswith(".cif") or basename.endswith(".cif.gz") or
                        basename.endswith(".pdb") or basename.endswith(".pdb.gz")):
                    continue

                dest_file = os.path.join(dest_dir, basename)
                if os.path.exists(dest_file) and os.path.getsize(dest_file) > 100:
                    skipped += 1
                else:
                    # Extract to flat directory (ignore tar subdirectories)
                    member.name = basename
                    tar.extract(member, dest_dir)
                    count += 1

                if i % 2000 == 0 or i == total:
                    elapsed = time.time() - t0
                    rate = i / elapsed if elapsed > 0 else 0
                    print(f"    [{i}/{total}] extracted={count} skipped={skipped} "
                          f"rate={rate:.0f} files/sec elapsed={elapsed:.0f}s")

    except (tarfile.TarError, EOFError) as e:
        print(f"  WARNING: Tar extraction error: {e}")
        print(f"  Extracted {count} files before error. Consider re-downloading.")
        return count

    elapsed = time.time() - t0
    print(f"  Extraction complete: {count} new + {skipped} existing = "
          f"{count + skipped} total in {elapsed:.0f}s")
    return count + skipped


# =============================================================================
# INDIVIDUAL DOWNLOAD (FALLBACK FOR MISSING ACCESSIONS)
# =============================================================================

def download_single_structure(accession, dest_dir, cfg):
    """Download a single AlphaFold structure (CIF). Try v4, then v6 fallback."""
    af_cfg = cfg.get("alphafold", {})
    urls = [
        af_cfg.get("individual_url_v4", "").format(acc=accession),
        af_cfg.get("individual_url_v6", "").format(acc=accession),
    ]

    for url in urls:
        if not url:
            continue
        fname = os.path.basename(url)
        fpath = os.path.join(dest_dir, fname)

        if os.path.exists(fpath) and os.path.getsize(fpath) > 100:
            return accession, "exists", fpath

        try:
            if requests is not None:
                resp = requests.get(url, timeout=30, stream=True)
                if resp.status_code == 200:
                    with open(fpath, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=65536):
                            f.write(chunk)
                    return accession, "downloaded", fpath
                elif resp.status_code == 404:
                    continue
            else:
                # Fallback to curl
                result = subprocess.run(
                    ["curl", "-sS", "-o", fpath, "-w", "%{http_code}", url],
                    capture_output=True, text=True, timeout=60,
                )
                if result.stdout.strip() == "200" and os.path.exists(fpath):
                    return accession, "downloaded", fpath
                elif os.path.exists(fpath):
                    os.remove(fpath)
        except Exception:
            continue

    return accession, "failed", None


def download_missing_individually(missing_accessions, dest_dir, cfg, checkpoint, workers=4):
    """
    Download missing structures individually using thread pool.
    Used as fallback for accessions not in bulk tar archive.
    """
    if not missing_accessions:
        return 0, 0

    # Filter out already-downloaded
    to_download = [
        acc for acc in missing_accessions
        if not checkpoint.is_individual_downloaded(acc)
    ]

    if not to_download:
        print(f"  All {len(missing_accessions)} accessions already handled")
        return len(missing_accessions), 0

    print(f"\n  Downloading {len(to_download)} missing structures individually "
          f"({workers} workers)...")
    t0 = time.time()
    succeeded = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_single_structure, acc, dest_dir, cfg): acc
            for acc in to_download
        }

        for i, future in enumerate(as_completed(futures), 1):
            acc, status, path = future.result()
            if status in ("downloaded", "exists"):
                succeeded += 1
                checkpoint.mark_individual_downloaded(acc)
            else:
                failed += 1
                checkpoint.add_failed(acc, "404_individual")

            if i % 100 == 0 or i == len(to_download):
                elapsed = time.time() - t0
                print(f"    [{i}/{len(to_download)}] ok={succeeded} fail={failed} "
                      f"elapsed={elapsed:.0f}s")

    checkpoint.save()
    return succeeded, failed


# =============================================================================
# VERIFICATION
# =============================================================================

def verify_downloaded_structures(dest_dir, expected_accessions, organism_name):
    """
    Verify downloaded structures: check file existence and basic integrity.
    Returns dict with statistics.
    """
    print(f"\n  Verifying {organism_name} structures in {dest_dir}...")

    # Scan directory for CIF files
    cif_pattern = re.compile(r"AF-([A-Z0-9]+)-F1-model_v\d+\.cif(\.gz)?$")
    found_accessions = {}
    total_size = 0

    if not os.path.isdir(dest_dir):
        print(f"  WARNING: Directory does not exist: {dest_dir}")
        return {"found": 0, "missing": len(expected_accessions), "total_size_gb": 0}

    for fname in os.listdir(dest_dir):
        m = cif_pattern.match(fname)
        if m:
            acc = m.group(1)
            fpath = os.path.join(dest_dir, fname)
            fsize = os.path.getsize(fpath)
            if fsize > 100:  # Not empty/corrupt
                found_accessions[acc] = fpath
                total_size += fsize

    found_set = set(found_accessions.keys())
    expected_set = set(expected_accessions)
    missing = expected_set - found_set
    extra = found_set - expected_set

    stats = {
        "expected": len(expected_set),
        "found": len(found_set & expected_set),
        "missing": len(missing),
        "extra": len(extra),
        "total_files": len(found_accessions),
        "total_size_gb": round(total_size / (1024 ** 3), 2),
        "coverage_pct": round(100 * len(found_set & expected_set) / max(len(expected_set), 1), 1),
    }

    print(f"    Expected:   {stats['expected']}")
    print(f"    Found:      {stats['found']}")
    print(f"    Missing:    {stats['missing']}")
    print(f"    Coverage:   {stats['coverage_pct']}%")
    print(f"    Disk usage: {stats['total_size_gb']} GB")

    if missing and len(missing) <= 20:
        print(f"    Missing accessions: {', '.join(sorted(missing))}")
    elif missing:
        print(f"    First 20 missing: {', '.join(sorted(missing)[:20])} ...")

    return stats, list(missing)


# =============================================================================
# BUILD STRUCTURE INDEX
# =============================================================================

def build_full_structure_index(cfg, organisms):
    """
    Build a TSV index of all downloaded structures with basic metadata.
    Analogous to Phase 1 structure_index.tsv but for full proteomes.
    """
    index_path = os.path.join(cfg["results"]["structures_dir"], "structure_index_full.tsv")
    os.makedirs(os.path.dirname(index_path), exist_ok=True)

    print(f"\nBuilding full structure index -> {index_path}")

    rows = []
    cif_pattern = re.compile(r"AF-([A-Z0-9]+)-F1-model_v(\d+)\.cif(\.gz)?$")

    for org_name, org_dir in organisms.items():
        if not os.path.isdir(org_dir):
            continue

        for fname in sorted(os.listdir(org_dir)):
            m = cif_pattern.match(fname)
            if not m:
                continue

            acc = m.group(1)
            version = f"v{m.group(2)}"
            fpath = os.path.join(org_dir, fname)
            fsize = os.path.getsize(fpath)

            rows.append({
                "uniprot_accession": acc,
                "organism": org_name,
                "model_version": version,
                "file_name": fname,
                "file_path": fpath,
                "file_size_bytes": fsize,
            })

    if rows:
        import csv
        with open(index_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"  Written {len(rows)} entries to {index_path}")
    else:
        print("  WARNING: No structures found to index")

    return index_path


# =============================================================================
# LOAD EXPECTED ACCESSIONS FROM UNIPROT PROTEOME TSV
# =============================================================================

def load_proteome_accessions(tsv_path):
    """
    Load UniProt accessions from proteome TSV file.
    Handles different column name conventions (Entry, accession, etc.).
    """
    accessions = []

    if not os.path.exists(tsv_path):
        print(f"  WARNING: Proteome TSV not found: {tsv_path}")
        return accessions

    import csv
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Try common column names
        for col in ["Entry", "entry", "accession", "Accession",
                     "uniprot_id", "uniprot_accession", "Entry Name"]:
            if col in reader.fieldnames:
                f.seek(0)
                next(reader)  # skip header
                for row in reader:
                    acc = row.get(col, "").strip()
                    if acc and len(acc) >= 6:
                        accessions.append(acc)
                break

    print(f"  Loaded {len(accessions)} accessions from {os.path.basename(tsv_path)}")
    return accessions


# =============================================================================
# MAIN
# =============================================================================

def process_organism(organism, cfg, checkpoint, workers=4):
    """
    Full download pipeline for one organism:
    1. Download bulk tar archive
    2. Extract CIF files
    3. Verify against expected proteome
    4. Download any missing individually
    5. Final verification
    """
    org_cfg = cfg["proteomes"][organism]
    af_cfg = cfg.get("alphafold", {})

    bulk_url = org_cfg["alphafold_bulk_url"]
    dest_dir = cfg["data"][f"alphafold_{organism}_dir"]
    tar_dir = os.path.join(cfg["data"]["raw_dir"], "alphafold_tars")
    tar_path = os.path.join(tar_dir, os.path.basename(bulk_url))

    print(f"\n{'='*70}")
    print(f"  Processing: {org_cfg['organism']}")
    print(f"  Expected proteins: {org_cfg['expected_proteins']}")
    print(f"  Estimated size: {org_cfg['estimated_size_gb']} GB")
    print(f"{'='*70}")

    # Step 1: Download bulk tar
    if not checkpoint.is_bulk_downloaded(organism):
        os.makedirs(tar_dir, exist_ok=True)
        print(f"\n[1] Downloading bulk tar archive...")
        success = download_bulk_tar(
            bulk_url, tar_path,
            max_retries=af_cfg.get("max_retries", 3),
            retry_delay=af_cfg.get("retry_delay_sec", 5),
        )
        if success:
            checkpoint.mark_bulk_downloaded(organism, tar_path)
        else:
            print(f"  FATAL: Could not download bulk archive for {organism}")
            print(f"  Falling back to individual downloads...")
            # Load expected accessions and download individually
            tsv_key = f"{organism}_proteome_tsv" if organism != "ecoli" else "ecoli_proteome_tsv"
            tsv_path = cfg["phase1_inputs"].get(tsv_key, "")
            expected = load_proteome_accessions(tsv_path)
            if expected:
                download_missing_individually(expected, dest_dir, cfg, checkpoint, workers)
            return
    else:
        print(f"\n[1] Bulk tar already downloaded (checkpoint)")

    # Step 2: Extract
    if not checkpoint.is_extracted(organism):
        print(f"\n[2] Extracting structures from tar archive...")
        count = extract_tar(tar_path, dest_dir, organism)
        checkpoint.mark_extracted(organism, count)
    else:
        print(f"\n[2] Already extracted (checkpoint)")

    # Step 3: Load expected accessions and verify
    print(f"\n[3] Loading expected accessions...")
    tsv_key = "ecoli_proteome_tsv" if organism == "ecoli" else "human_proteome_tsv"
    tsv_path = cfg["phase1_inputs"].get(tsv_key, "")
    expected = load_proteome_accessions(tsv_path)

    if not expected:
        print(f"  WARNING: No expected accessions loaded. Skipping verification.")
        return

    stats, missing = verify_downloaded_structures(dest_dir, expected, organism)

    # Step 4: Download missing individually
    if missing:
        print(f"\n[4] Downloading {len(missing)} missing structures individually...")
        ok, fail = download_missing_individually(missing, dest_dir, cfg, checkpoint, workers)
        print(f"    Individual downloads: {ok} succeeded, {fail} failed")

        # Re-verify
        stats, still_missing = verify_downloaded_structures(dest_dir, expected, organism)
    else:
        still_missing = []

    # Step 5: Save verification results
    checkpoint.mark_verified(organism, stats)

    if still_missing:
        missing_path = os.path.join(
            cfg["results"]["structures_dir"],
            f"missing_accessions_{organism}.txt",
        )
        os.makedirs(os.path.dirname(missing_path), exist_ok=True)
        with open(missing_path, "w") as f:
            for acc in sorted(still_missing):
                f.write(acc + "\n")
        print(f"  Missing accessions written to: {missing_path}")

    print(f"\n  {organism} COMPLETE: {stats.get('coverage_pct', 0)}% coverage")


def main():
    parser = argparse.ArgumentParser(
        description="Download full AlphaFold proteome structures (Phase 2)"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--organism",
        choices=["ecoli", "human", "both"],
        default="both",
        help="Which organism(s) to download (default: both)",
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel download workers for individual fallback (default: 4)",
    )
    parser.add_argument(
        "--skip-bulk", action="store_true",
        help="Skip bulk tar download; only do individual downloads",
    )
    parser.add_argument(
        "--verify-only", action="store_true",
        help="Only verify existing downloads, don't download anything new",
    )

    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    # Set up checkpoint
    checkpoint_path = os.path.join(
        cfg["results"]["structures_dir"], "download_checkpoint.json"
    )
    os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)
    checkpoint = CheckpointManager(checkpoint_path)

    print("=" * 70)
    print("  PHASE 2: AlphaFold Full Proteome Download")
    print(f"  Config: {args.config}")
    print(f"  Organism(s): {args.organism}")
    print(f"  Checkpoint: {checkpoint_path}")
    print("=" * 70)

    # Determine organisms to process
    organisms = []
    if args.organism in ("ecoli", "both"):
        organisms.append("ecoli")
    if args.organism in ("human", "both"):
        organisms.append("human")

    for organism in organisms:
        if args.verify_only:
            dest_dir = cfg["data"][f"alphafold_{organism}_dir"]
            tsv_key = "ecoli_proteome_tsv" if organism == "ecoli" else "human_proteome_tsv"
            tsv_path = cfg["phase1_inputs"].get(tsv_key, "")
            expected = load_proteome_accessions(tsv_path)
            if expected:
                verify_downloaded_structures(dest_dir, expected, organism)
        else:
            process_organism(organism, cfg, checkpoint, args.workers)

    # Build structure index
    if not args.verify_only:
        organisms_dirs = {}
        for org in organisms:
            organisms_dirs[org] = cfg["data"][f"alphafold_{org}_dir"]
        build_full_structure_index(cfg, organisms_dirs)

    # Summary
    print(f"\n{'='*70}")
    print("  DOWNLOAD SUMMARY")
    print(f"{'='*70}")
    for org in organisms:
        verified = checkpoint.state.get("verified", {}).get(org, {})
        print(f"  {org}:")
        print(f"    Coverage: {verified.get('coverage_pct', 'N/A')}%")
        print(f"    Found: {verified.get('found', 'N/A')}/{verified.get('expected', 'N/A')}")
        print(f"    Disk: {verified.get('total_size_gb', 'N/A')} GB")

    if checkpoint.state.get("failed"):
        n_failed = len(checkpoint.state["failed"])
        print(f"\n  Total failed downloads: {n_failed}")
        if n_failed <= 10:
            for entry in checkpoint.state["failed"]:
                print(f"    {entry['accession']}: {entry['reason']}")

    print(f"\n  Checkpoint saved: {checkpoint_path}")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
