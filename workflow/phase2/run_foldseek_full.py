#!/usr/bin/env python3
"""
Phase 2: Full-scale Foldseek all-vs-all structural search and clustering.

Antah Asti Prarambh — Module E (Full-Scale)

Runs Foldseek on complete E. coli + Human proteome structures (~25,000 proteins).
This requires significant RAM (32-64 GB) and is designed for HPC execution.

Pipeline:
  1. Create Foldseek databases from AlphaFold CIF directories
  2. Run all-vs-all structural search (within and between organisms)
  3. Cluster search results
  4. Parse and annotate cluster membership with substrate labels
  5. Generate summary report

Estimated resources:
  - createdb:    16 GB RAM, 4 CPUs, ~30 min
  - search:      64 GB RAM, 16 CPUs, ~6-24 hours
  - cluster:     64 GB RAM, 16 CPUs, ~2-4 hours
  - analysis:    16 GB RAM, 4 CPUs, ~30 min
  - Disk (temp): ~50 GB for Foldseek temporary databases

Usage:
  python run_foldseek_full.py --config config.yaml
  python run_foldseek_full.py --config config.yaml --step search --threads 32
  python run_foldseek_full.py --config config.yaml --step analyze-only

SLURM example:
  sbatch --mem=64G --cpus-per-task=16 --time=24:00:00 \\
    --wrap "python run_foldseek_full.py --config config.yaml --threads 16"
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
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
# UTILITY
# =============================================================================

def run_cmd(cmd, description, timeout=None, env=None):
    """Run a shell command with logging."""
    print(f"\n  [{time.strftime('%H:%M:%S')}] {description}")
    print(f"    Command: {' '.join(cmd)}")
    t0 = time.time()

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )

    elapsed = time.time() - t0
    print(f"    Completed in {elapsed:.0f}s (exit code {result.returncode})")

    if result.returncode != 0:
        print(f"    STDERR: {result.stderr[:2000]}")
        raise RuntimeError(f"Command failed: {' '.join(cmd[:5])}...")

    if result.stdout:
        # Print last few lines of stdout as summary
        lines = result.stdout.strip().split("\n")
        if len(lines) > 5:
            print(f"    Output ({len(lines)} lines, last 3):")
            for line in lines[-3:]:
                print(f"      {line}")
        else:
            for line in lines:
                print(f"      {line}")

    return result


def extract_accession(name):
    """Extract UniProt accession from Foldseek database entry name."""
    m = re.match(r"AF-([A-Z0-9]+)-F1", name)
    if m:
        return m.group(1)
    # Try just the name itself
    m = re.match(r"^([A-Z][A-Z0-9]{4,9})$", name)
    if m:
        return m.group(1)
    return name


def check_foldseek(binary="foldseek"):
    """Verify Foldseek is available and return version."""
    try:
        result = subprocess.run([binary, "version"], capture_output=True, text=True)
        version = result.stdout.strip()
        print(f"  Foldseek version: {version}")
        return version
    except FileNotFoundError:
        print(f"  ERROR: Foldseek not found at '{binary}'")
        print(f"  Install with: conda install -c conda-forge -c bioconda foldseek")
        sys.exit(1)


# =============================================================================
# STEP 1: CREATE FOLDSEEK DATABASES
# =============================================================================

def create_database(foldseek_bin, struct_dir, db_path, threads=4, timeout=7200):
    """
    Create a Foldseek structure database from a directory of CIF/PDB files.
    """
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    # Check if database already exists
    if os.path.exists(db_path) and os.path.exists(db_path + ".index"):
        # Verify it's not empty
        size = os.path.getsize(db_path)
        if size > 1000:
            print(f"  Database already exists: {db_path} ({size / 1e6:.1f} MB)")
            return True

    cmd = [
        foldseek_bin, "createdb",
        struct_dir,
        db_path,
        "--threads", str(threads),
    ]

    run_cmd(cmd, f"Creating Foldseek DB from {struct_dir}", timeout=timeout)
    return True


# =============================================================================
# STEP 2: ALL-VS-ALL STRUCTURAL SEARCH
# =============================================================================

def run_search(foldseek_bin, query_db, target_db, result_db, tmp_dir,
               threads=16, sensitivity=7.5, evalue=0.001,
               split_memory="64G", timeout=86400):
    """
    Run Foldseek all-vs-all structural search.
    This is the most resource-intensive step.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(result_db), exist_ok=True)

    # Check if result already exists
    if os.path.exists(result_db) and os.path.exists(result_db + ".index"):
        size = os.path.getsize(result_db)
        if size > 1000:
            print(f"  Search results already exist: {result_db} ({size / 1e6:.1f} MB)")
            return True

    cmd = [
        foldseek_bin, "search",
        query_db, target_db, result_db, tmp_dir,
        "--threads", str(threads),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--split-memory-limit", split_memory,
        "--alignment-type", "2",  # 3Di+AA (recommended for structural search)
        "-a",  # add backtrace (alignment) for coverage calculation
    ]

    run_cmd(cmd, "Running Foldseek all-vs-all search", timeout=timeout)
    return True


# =============================================================================
# STEP 3: CLUSTERING
# =============================================================================

def run_clustering(foldseek_bin, db_path, result_db, cluster_db, tmp_dir,
                   min_seq_id=0.3, coverage=0.5, coverage_mode=0,
                   cluster_mode=0, threads=16, timeout=28800):
    """
    Cluster Foldseek search results.
    """
    os.makedirs(os.path.dirname(cluster_db), exist_ok=True)

    # Check if cluster already exists
    if os.path.exists(cluster_db) and os.path.exists(cluster_db + ".index"):
        size = os.path.getsize(cluster_db)
        if size > 100:
            print(f"  Cluster results already exist: {cluster_db} ({size / 1e6:.1f} MB)")
            return True

    cmd = [
        foldseek_bin, "cluster",
        db_path, cluster_db, tmp_dir,
        "--threads", str(threads),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(coverage_mode),
        "--cluster-mode", str(cluster_mode),
    ]

    run_cmd(cmd, "Running Foldseek clustering", timeout=timeout)
    return True


# =============================================================================
# STEP 4: EXPORT RESULTS TO TSV
# =============================================================================

def export_results(foldseek_bin, db_path, result_db, output_tsv, tmp_dir,
                   db2_path=None, timeout=3600):
    """
    Export Foldseek results to human-readable TSV format.
    """
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

    if os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 100:
        print(f"  TSV already exists: {output_tsv}")
        return True

    target = db2_path if db2_path else db_path

    cmd = [
        foldseek_bin, "convertalis",
        db_path, target, result_db, output_tsv,
        "--format-output",
        "query,target,evalue,bits,pident,alnlen,qcov,tcov,lddt,prob",
    ]

    run_cmd(cmd, f"Exporting results to {output_tsv}", timeout=timeout)
    return True


def export_cluster_membership(foldseek_bin, db_path, cluster_db, output_tsv,
                               timeout=3600):
    """Export cluster membership table (representative -> member)."""
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

    if os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 100:
        print(f"  Cluster membership TSV already exists: {output_tsv}")
        return True

    cmd = [
        foldseek_bin, "createtsv",
        db_path, db_path, cluster_db, output_tsv,
    ]

    run_cmd(cmd, f"Exporting cluster membership to {output_tsv}", timeout=timeout)
    return True


# =============================================================================
# STEP 5: ANALYZE CLUSTERS
# =============================================================================

def load_substrate_labels(cfg):
    """
    Load substrate/dataset labels from Phase 1 processed files.
    Returns dict: accession -> set of labels.
    """
    labels = defaultdict(set)

    # GroEL substrates
    groel_path = cfg["phase1_inputs"]["groel_substrates"]
    if os.path.exists(groel_path):
        with open(groel_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acc = row.get("current_accession", "").strip()
                if acc:
                    labels[acc].add("groel")
                    klass = row.get("class", "").strip()
                    if klass:
                        labels[acc].add(f"groel_class_{klass}")

    # HSP60 substrates
    hsp60_path = cfg["phase1_inputs"]["hsp60_substrates"]
    if os.path.exists(hsp60_path):
        with open(hsp60_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acc = row.get("uniprot_id", "").strip()
                if acc:
                    labels[acc].add("hsp60")

    # Mito proteome
    mito_path = cfg["phase1_inputs"]["mito_proteome"]
    if os.path.exists(mito_path):
        with open(mito_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acc = row.get("uniprot_id", "").strip()
                if acc:
                    labels[acc].add("mito")
                    if row.get("is_matrix", "0") == "1":
                        labels[acc].add("matrix")

    print(f"  Loaded labels: {sum(1 for v in labels.values() if 'groel' in v)} GroEL, "
          f"{sum(1 for v in labels.values() if 'hsp60' in v)} HSP60, "
          f"{sum(1 for v in labels.values() if 'mito' in v)} mito")

    return labels


def analyze_clusters(cluster_tsv, search_tsv, labels, output_dir):
    """
    Analyze Foldseek clustering results:
    - Cluster size distribution
    - Substrate enrichment per cluster
    - Shared clusters between GroEL and HSP60
    - Summary statistics
    """
    print(f"\n  Analyzing clusters from {cluster_tsv}...")
    os.makedirs(output_dir, exist_ok=True)

    # Parse cluster membership
    member_to_rep = {}
    rep_to_members = defaultdict(list)

    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            rep_acc = extract_accession(parts[0])
            member_acc = extract_accession(parts[1])
            member_to_rep[member_acc] = rep_acc
            rep_to_members[rep_acc].append(member_acc)

    n_proteins = len(member_to_rep)
    n_clusters = len(rep_to_members)

    print(f"    Total proteins: {n_proteins}")
    print(f"    Total clusters: {n_clusters}")

    # Assign numeric IDs (sorted by size descending)
    sorted_reps = sorted(rep_to_members.keys(),
                         key=lambda r: (-len(rep_to_members[r]), r))
    rep_to_id = {rep: i + 1 for i, rep in enumerate(sorted_reps)}

    # Size distribution
    sizes = [len(members) for members in rep_to_members.values()]
    size_dist = Counter(sizes)
    singletons = sum(1 for s in sizes if s == 1)
    small = sum(1 for s in sizes if 2 <= s <= 5)
    medium = sum(1 for s in sizes if 6 <= s <= 20)
    large = sum(1 for s in sizes if 21 <= s <= 100)
    xlarge = sum(1 for s in sizes if s > 100)

    # Per-protein output table
    rows = []
    for member_acc in sorted(member_to_rep.keys()):
        rep_acc = member_to_rep[member_acc]
        cluid = rep_to_id[rep_acc]
        csize = len(rep_to_members[rep_acc])
        member_labels = labels.get(member_acc, set())
        label_str = ",".join(sorted(member_labels)) if member_labels else "background"

        # Determine organism
        organism = "unknown"
        if "groel" in member_labels or any(l.startswith("groel_class") for l in member_labels):
            organism = "ecoli"
        elif "hsp60" in member_labels or "mito" in member_labels or "matrix" in member_labels:
            organism = "human"

        rows.append({
            "uniprot_accession": member_acc,
            "cluster_id": cluid,
            "cluster_representative": rep_acc,
            "cluster_size": csize,
            "organism": organism,
            "dataset_labels": label_str,
        })

    out_table = os.path.join(output_dir, "foldseek_clusters_full.tsv")
    with open(out_table, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"    Cluster table written: {out_table} ({len(rows)} rows)")

    # Shared cluster analysis
    groel_clusters = set()
    hsp60_clusters = set()
    for member_acc, rep_acc in member_to_rep.items():
        member_labels = labels.get(member_acc, set())
        cluid = rep_to_id[rep_acc]
        if "groel" in member_labels:
            groel_clusters.add(cluid)
        if "hsp60" in member_labels:
            hsp60_clusters.add(cluid)

    shared = groel_clusters & hsp60_clusters
    groel_only = groel_clusters - hsp60_clusters
    hsp60_only = hsp60_clusters - groel_clusters

    # Generate summary report
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("  FOLDSEEK FULL-SCALE CLUSTERING — PHASE 2 SUMMARY")
    report_lines.append("=" * 70)
    report_lines.append("")
    report_lines.append("CLUSTERING PARAMETERS:")
    report_lines.append("  Method: Foldseek cluster")
    report_lines.append("  Scope: Full E. coli + Human proteomes")
    report_lines.append("")
    report_lines.append("OVERALL STATISTICS:")
    report_lines.append(f"  Total proteins:    {n_proteins}")
    report_lines.append(f"  Total clusters:    {n_clusters}")
    report_lines.append(f"  Singletons (1):    {singletons} ({100*singletons/max(n_clusters,1):.1f}%)")
    report_lines.append(f"  Small (2-5):       {small}")
    report_lines.append(f"  Medium (6-20):     {medium}")
    report_lines.append(f"  Large (21-100):    {large}")
    report_lines.append(f"  Very large (>100): {xlarge}")
    report_lines.append("")
    report_lines.append("SUBSTRATE CLUSTER OVERLAP:")
    report_lines.append(f"  GroEL substrate clusters:  {len(groel_clusters)}")
    report_lines.append(f"  HSP60 substrate clusters:  {len(hsp60_clusters)}")
    report_lines.append(f"  Shared clusters:           {len(shared)}")
    report_lines.append(f"  GroEL-only clusters:       {len(groel_only)}")
    report_lines.append(f"  HSP60-only clusters:       {len(hsp60_only)}")
    report_lines.append("")

    # Top 30 largest clusters
    report_lines.append("TOP 30 LARGEST CLUSTERS:")
    report_lines.append(f"  {'ID':>6}  {'Size':>6}  {'Rep':>12}  Labels")
    report_lines.append(f"  {'-'*6}  {'-'*6}  {'-'*12}  {'-'*40}")
    for rep in sorted_reps[:30]:
        cluid = rep_to_id[rep]
        members = rep_to_members[rep]
        csize = len(members)
        # Aggregate labels
        all_labels = set()
        for m in members:
            all_labels |= labels.get(m, set())
        label_str = ",".join(sorted(all_labels)) if all_labels else "background"
        report_lines.append(f"  {cluid:>6}  {csize:>6}  {rep:>12}  {label_str[:60]}")

    report_lines.append("")
    report_lines.append("=" * 70)

    out_report = os.path.join(output_dir, "foldseek_full_summary.txt")
    with open(out_report, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"    Report written: {out_report}")

    # Print report
    print("\n" + "\n".join(report_lines))

    return out_table, out_report


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Full-scale Foldseek structural clustering (Phase 2)"
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--step",
        choices=["all", "createdb", "search", "cluster", "export", "analyze-only"],
        default="all",
        help="Which pipeline step to run (default: all)",
    )
    parser.add_argument("--threads", type=int, default=None, help="Override thread count")
    parser.add_argument("--memory", default=None, help="Override split memory limit (e.g. 64G)")
    parser.add_argument(
        "--scope",
        choices=["combined", "ecoli-only", "human-only"],
        default="combined",
        help="Which proteome(s) to cluster (default: combined)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    fs_cfg = cfg.get("foldseek", {})

    foldseek_bin = fs_cfg.get("binary", "foldseek")
    threads = args.threads or fs_cfg.get("threads", 16)
    split_memory = args.memory or fs_cfg.get("split_memory_limit", "64G")
    tmp_dir = fs_cfg.get("tmp_dir", "/tmp/foldseek_tmp")

    # Paths
    foldseek_dir = cfg["results"]["foldseek_dir"]
    db_dir = os.path.join(foldseek_dir, "databases")
    search_dir = os.path.join(foldseek_dir, "search")
    cluster_dir = os.path.join(foldseek_dir, "clusters")
    analysis_dir = os.path.join(foldseek_dir, "analysis")

    for d in [db_dir, search_dir, cluster_dir, analysis_dir]:
        os.makedirs(d, exist_ok=True)

    print("=" * 70)
    print("  PHASE 2: Foldseek Full-Scale Structural Clustering")
    print(f"  Threads: {threads}")
    print(f"  Memory limit: {split_memory}")
    print(f"  Scope: {args.scope}")
    print(f"  Step: {args.step}")
    print("=" * 70)

    version = check_foldseek(foldseek_bin)

    # Define database paths based on scope
    ecoli_struct_dir = cfg["data"]["alphafold_ecoli_dir"]
    human_struct_dir = cfg["data"]["alphafold_human_dir"]

    ecoli_db = os.path.join(db_dir, "ecoli_db")
    human_db = os.path.join(db_dir, "human_db")
    combined_db = os.path.join(db_dir, "combined_db")

    # ── Step 1: Create databases ──
    if args.step in ("all", "createdb"):
        print(f"\n{'='*70}")
        print("  STEP 1: Creating Foldseek databases")
        print(f"{'='*70}")

        if args.scope in ("combined", "ecoli-only"):
            create_database(foldseek_bin, ecoli_struct_dir, ecoli_db, threads)

        if args.scope in ("combined", "human-only"):
            create_database(foldseek_bin, human_struct_dir, human_db, threads)

        if args.scope == "combined":
            # Create combined database by merging (or just create from both dirs)
            # Foldseek doesn't have a merge command; we create a symlink dir
            combined_struct_dir = os.path.join(foldseek_dir, "combined_structures")
            os.makedirs(combined_struct_dir, exist_ok=True)

            # Symlink all structures into combined directory
            for src_dir in [ecoli_struct_dir, human_struct_dir]:
                if os.path.isdir(src_dir):
                    for fname in os.listdir(src_dir):
                        if fname.endswith(".cif") or fname.endswith(".cif.gz"):
                            src = os.path.join(src_dir, fname)
                            dst = os.path.join(combined_struct_dir, fname)
                            if not os.path.exists(dst):
                                os.symlink(src, dst)

            create_database(foldseek_bin, combined_struct_dir, combined_db, threads)

    # ── Step 2: All-vs-all search ──
    if args.step in ("all", "search"):
        print(f"\n{'='*70}")
        print("  STEP 2: All-vs-all structural search")
        print(f"{'='*70}")

        if args.scope == "combined":
            result_db = os.path.join(search_dir, "combined_search")
            run_search(
                foldseek_bin, combined_db, combined_db, result_db, tmp_dir,
                threads=threads, sensitivity=fs_cfg.get("sensitivity", 7.5),
                evalue=fs_cfg.get("search_evalue", 0.001),
                split_memory=split_memory,
            )
        elif args.scope == "ecoli-only":
            result_db = os.path.join(search_dir, "ecoli_search")
            run_search(
                foldseek_bin, ecoli_db, ecoli_db, result_db, tmp_dir,
                threads=threads, sensitivity=fs_cfg.get("sensitivity", 7.5),
                evalue=fs_cfg.get("search_evalue", 0.001),
                split_memory=split_memory,
            )
        elif args.scope == "human-only":
            result_db = os.path.join(search_dir, "human_search")
            run_search(
                foldseek_bin, human_db, human_db, result_db, tmp_dir,
                threads=threads, sensitivity=fs_cfg.get("sensitivity", 7.5),
                evalue=fs_cfg.get("search_evalue", 0.001),
                split_memory=split_memory,
            )

    # ── Step 3: Clustering ──
    if args.step in ("all", "cluster"):
        print(f"\n{'='*70}")
        print("  STEP 3: Clustering")
        print(f"{'='*70}")

        if args.scope == "combined":
            cluster_db_path = os.path.join(cluster_dir, "combined_cluster")
            run_clustering(
                foldseek_bin, combined_db,
                os.path.join(search_dir, "combined_search"),
                cluster_db_path, tmp_dir,
                min_seq_id=fs_cfg.get("min_seq_id", 0.3),
                coverage=fs_cfg.get("coverage", 0.5),
                coverage_mode=fs_cfg.get("coverage_mode", 0),
                cluster_mode=fs_cfg.get("cluster_mode", 0),
                threads=threads,
            )
        elif args.scope == "ecoli-only":
            cluster_db_path = os.path.join(cluster_dir, "ecoli_cluster")
            run_clustering(
                foldseek_bin, ecoli_db,
                os.path.join(search_dir, "ecoli_search"),
                cluster_db_path, tmp_dir,
                min_seq_id=fs_cfg.get("min_seq_id", 0.3),
                coverage=fs_cfg.get("coverage", 0.5),
                threads=threads,
            )
        elif args.scope == "human-only":
            cluster_db_path = os.path.join(cluster_dir, "human_cluster")
            run_clustering(
                foldseek_bin, human_db,
                os.path.join(search_dir, "human_search"),
                cluster_db_path, tmp_dir,
                min_seq_id=fs_cfg.get("min_seq_id", 0.3),
                coverage=fs_cfg.get("coverage", 0.5),
                threads=threads,
            )

    # ── Step 4: Export to TSV ──
    if args.step in ("all", "export"):
        print(f"\n{'='*70}")
        print("  STEP 4: Exporting results to TSV")
        print(f"{'='*70}")

        scope_prefix = args.scope.replace("-only", "")

        # Export search results
        search_tsv = os.path.join(analysis_dir, f"{scope_prefix}_search_results.tsv")
        db_path = combined_db if args.scope == "combined" else (
            ecoli_db if args.scope == "ecoli-only" else human_db
        )
        result_db = os.path.join(search_dir, f"{scope_prefix}_search")

        if os.path.exists(result_db + ".index"):
            export_results(foldseek_bin, db_path, result_db, search_tsv, tmp_dir)

        # Export cluster membership
        cluster_tsv = os.path.join(analysis_dir, f"{scope_prefix}_cluster_membership.tsv")
        cluster_db_path = os.path.join(cluster_dir, f"{scope_prefix}_cluster")

        if os.path.exists(cluster_db_path + ".index"):
            export_cluster_membership(
                foldseek_bin, db_path, cluster_db_path, cluster_tsv
            )

    # ── Step 5: Analyze ──
    if args.step in ("all", "analyze-only"):
        print(f"\n{'='*70}")
        print("  STEP 5: Analyzing clusters")
        print(f"{'='*70}")

        labels = load_substrate_labels(cfg)

        scope_prefix = args.scope.replace("-only", "") if args.step != "analyze-only" else "combined"
        cluster_tsv = os.path.join(analysis_dir, f"{scope_prefix}_cluster_membership.tsv")
        search_tsv = os.path.join(analysis_dir, f"{scope_prefix}_search_results.tsv")

        if os.path.exists(cluster_tsv):
            analyze_clusters(cluster_tsv, search_tsv, labels, analysis_dir)
        else:
            print(f"  WARNING: Cluster membership file not found: {cluster_tsv}")
            print(f"  Run with --step all or --step export first.")

    print(f"\n{'='*70}")
    print("  Foldseek pipeline complete")
    print(f"{'='*70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
