#!/usr/bin/env python3
"""
Module C - Step C4: Orthology group analysis between E. coli and Human proteomes.

Uses MMseqs2 all-vs-all search + connected-component clustering to build orthogroups,
then intersects with GroEL/HSP60 substrate lists.

Fallback approach since OrthoFinder is not installed.
"""

import subprocess
import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from datetime import datetime
import tempfile
import shutil

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────
BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"
MMSEQS = "/Users/vishalbharti/opt/anaconda3/envs/proteomics/bin/mmseqs"

ECOLI_FASTA = f"{BASE}/data/raw/uniprot/ecoli_k12_proteome.fasta"
HUMAN_FASTA = f"{BASE}/data/raw/uniprot/human_proteome.fasta"

GROEL_FILE = f"{BASE}/data/processed/groel_substrates_standardized.tsv"
HSP60_FILE = f"{BASE}/data/processed/hsp60_tier1_substrates.tsv"
RBH_FILE   = f"{BASE}/results/homology/rbh_groel_hsp60_annotated.tsv"

OUT_DIR    = f"{BASE}/results/homology"
WORK_DIR   = f"{BASE}/results/homology/_mmseqs_ortho_work"

# Thresholds
EVALUE_CUTOFF = 1e-5
MIN_COVERAGE  = 0.50   # minimum query or target coverage
MIN_PIDENT    = 25.0    # minimum percent identity

# ──────────────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────────────

def run_cmd(cmd, desc=""):
    """Run a shell command and check for errors."""
    print(f"  >> {desc}: {cmd[:120]}..." if len(cmd) > 120 else f"  >> {desc}: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  !! ERROR (rc={result.returncode}): {result.stderr[:500]}")
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr[:500]}")
    return result.stdout


def parse_fasta_ids(fasta_path):
    """Extract accession IDs from UniProt FASTA headers (sp|ACCESSION|...)."""
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # >sp|P0AFG3|SUCA_ECOLI ... or >tr|...
                parts = line[1:].split("|")
                if len(parts) >= 2:
                    ids.append(parts[1])
                else:
                    ids.append(line[1:].split()[0])
    return set(ids)


def union_find_cluster(edges):
    """
    Simple union-find to compute connected components from edge list.
    edges: list of (nodeA, nodeB)
    Returns: dict {node: component_id}
    """
    parent = {}

    def find(x):
        if x not in parent:
            parent[x] = x
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        union(a, b)

    # Build components
    components = defaultdict(set)
    for node in parent:
        components[find(node)].add(node)

    # Assign numeric IDs
    node_to_comp = {}
    for idx, (_, members) in enumerate(sorted(components.items(), key=lambda x: -len(x[1]))):
        for m in members:
            node_to_comp[m] = idx
    return node_to_comp, components


# ──────────────────────────────────────────────────────────────────────
# STEP 1: MMseqs2 all-vs-all searches
# ──────────────────────────────────────────────────────────────────────

def run_mmseqs_searches():
    """Run bidirectional all-vs-all MMseqs2 searches."""
    os.makedirs(WORK_DIR, exist_ok=True)
    tmpdir = os.path.join(WORK_DIR, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    ecoli_db = os.path.join(WORK_DIR, "ecoli_db")
    human_db = os.path.join(WORK_DIR, "human_db")

    ecoli_vs_human_db = os.path.join(WORK_DIR, "ecoli_vs_human")
    human_vs_ecoli_db = os.path.join(WORK_DIR, "human_vs_ecoli")

    ecoli_vs_human_tsv = os.path.join(WORK_DIR, "ecoli_vs_human.tsv")
    human_vs_ecoli_tsv = os.path.join(WORK_DIR, "human_vs_ecoli.tsv")

    fmt = "query,target,evalue,bits,pident,alnlen,qcov,tcov"

    # Create databases
    print("\n[1/6] Creating MMseqs2 databases...")
    run_cmd(f"{MMSEQS} createdb {ECOLI_FASTA} {ecoli_db}", "Create E. coli DB")
    run_cmd(f"{MMSEQS} createdb {HUMAN_FASTA} {human_db}", "Create Human DB")

    # E. coli -> Human search
    print("\n[2/6] Searching E. coli -> Human...")
    run_cmd(
        f"{MMSEQS} search {ecoli_db} {human_db} {ecoli_vs_human_db} {tmpdir} "
        f"-e {EVALUE_CUTOFF} --max-seqs 300 --threads 4 -s 7.5",
        "E. coli vs Human search"
    )
    run_cmd(
        f"{MMSEQS} convertalis {ecoli_db} {human_db} {ecoli_vs_human_db} {ecoli_vs_human_tsv} "
        f'--format-output "{fmt}"',
        "Convert E->H alignments"
    )

    # Human -> E. coli search
    print("\n[3/6] Searching Human -> E. coli...")
    run_cmd(
        f"{MMSEQS} search {human_db} {ecoli_db} {human_vs_ecoli_db} {tmpdir} "
        f"-e {EVALUE_CUTOFF} --max-seqs 300 --threads 4 -s 7.5",
        "Human vs E. coli search"
    )
    run_cmd(
        f"{MMSEQS} convertalis {human_db} {ecoli_db} {human_vs_ecoli_db} {human_vs_ecoli_tsv} "
        f'--format-output "{fmt}"',
        "Convert H->E alignments"
    )

    return ecoli_vs_human_tsv, human_vs_ecoli_tsv


# ──────────────────────────────────────────────────────────────────────
# STEP 2: Build orthogroups from reciprocal hits
# ──────────────────────────────────────────────────────────────────────

def build_orthogroups(e2h_tsv, h2e_tsv):
    """
    Build orthogroups using reciprocal best hit expansion + connected components.

    Strategy:
    1. Load all hits from both directions
    2. Filter by e-value, coverage, and identity thresholds
    3. Keep only reciprocal hits (A->B exists AND B->A exists)
    4. Cluster reciprocal hits into connected components = orthogroups
    """
    print("\n[4/6] Building orthogroups...")

    cols = ["query", "target", "evalue", "bits", "pident", "alnlen", "qcov", "tcov"]

    # Load hits
    df_e2h = pd.read_csv(e2h_tsv, sep="\t", header=None, names=cols)
    df_h2e = pd.read_csv(h2e_tsv, sep="\t", header=None, names=cols)

    print(f"  Raw hits: E->H = {len(df_e2h):,}, H->E = {len(df_h2e):,}")

    # Filter by thresholds
    def apply_filters(df):
        mask = (
            (df["evalue"] <= EVALUE_CUTOFF) &
            (df["pident"] >= MIN_PIDENT) &
            ((df["qcov"] >= MIN_COVERAGE) | (df["tcov"] >= MIN_COVERAGE))
        )
        return df[mask].copy()

    df_e2h_f = apply_filters(df_e2h)
    df_h2e_f = apply_filters(df_h2e)
    print(f"  After filtering: E->H = {len(df_e2h_f):,}, H->E = {len(df_h2e_f):,}")

    # Build pair sets for reciprocal check
    # E->H: pairs are (ecoli_id, human_id)
    e2h_pairs = set(zip(df_e2h_f["query"], df_e2h_f["target"]))
    # H->E: pairs are (human_id, ecoli_id), flip to (ecoli_id, human_id)
    h2e_pairs = set(zip(df_h2e_f["target"], df_h2e_f["query"]))

    # Reciprocal hits: pair exists in both directions
    reciprocal_pairs = e2h_pairs & h2e_pairs
    print(f"  Reciprocal pairs (both directions): {len(reciprocal_pairs):,}")

    # Also compute best-hit-per-query for comparison
    # For each E. coli protein, its best human hit
    e2h_best = df_e2h_f.sort_values("evalue").groupby("query").first().reset_index()
    # For each human protein, its best E. coli hit
    h2e_best = df_h2e_f.sort_values("evalue").groupby("query").first().reset_index()

    rbh_pairs_from_best = set()
    for _, row in e2h_best.iterrows():
        ec, hu = row["query"], row["target"]
        # Check if human's best hit is this E. coli protein
        hu_best = h2e_best[h2e_best["query"] == hu]
        if len(hu_best) > 0 and hu_best.iloc[0]["target"] == ec:
            rbh_pairs_from_best.add((ec, hu))

    print(f"  Strict RBH pairs (best-hit-only): {len(rbh_pairs_from_best):,}")

    # Build stats for each reciprocal pair (use best stats from e2h direction)
    pair_stats = {}
    for _, row in df_e2h_f.iterrows():
        pair = (row["query"], row["target"])
        if pair in reciprocal_pairs:
            if pair not in pair_stats or row["evalue"] < pair_stats[pair]["evalue"]:
                pair_stats[pair] = {
                    "evalue": row["evalue"],
                    "bits": row["bits"],
                    "pident": row["pident"],
                    "alnlen": row["alnlen"],
                    "qcov": row["qcov"],
                    "tcov": row["tcov"]
                }

    # Cluster into orthogroups using connected components
    # Edges: each reciprocal pair is an edge between ecoli_id and human_id
    # We prefix IDs to distinguish species: "EC:xxx" and "HS:xxx"
    edges = []
    for ec, hu in reciprocal_pairs:
        edges.append((f"EC:{ec}", f"HS:{hu}"))

    node_to_comp, components = union_find_cluster(edges)

    # Build orthogroup table
    orthogroups = []
    for comp_id, members in sorted(components.items(), key=lambda x: -len(x[1])):
        ec_members = sorted([m[3:] for m in members if m.startswith("EC:")])
        hs_members = sorted([m[3:] for m in members if m.startswith("HS:")])

        if not ec_members or not hs_members:
            continue  # skip single-species groups

        # Get best stats for this orthogroup
        best_evalue = float("inf")
        pidents = []
        for ec in ec_members:
            for hu in hs_members:
                pair = (ec, hu)
                if pair in pair_stats:
                    if pair_stats[pair]["evalue"] < best_evalue:
                        best_evalue = pair_stats[pair]["evalue"]
                    pidents.append(pair_stats[pair]["pident"])

        og_id = node_to_comp[f"EC:{ec_members[0]}"]

        orthogroups.append({
            "orthogroup_id": f"OG{og_id:05d}",
            "ecoli_accessions": ",".join(ec_members),
            "human_accessions": ",".join(hs_members),
            "n_ecoli": len(ec_members),
            "n_human": len(hs_members),
            "best_evalue": best_evalue,
            "mean_pident": round(np.mean(pidents), 1) if pidents else 0.0
        })

    df_og = pd.DataFrame(orthogroups)
    # Re-sort and re-number
    df_og = df_og.sort_values("best_evalue").reset_index(drop=True)
    df_og["orthogroup_id"] = [f"OG{i:05d}" for i in range(len(df_og))]

    print(f"  Orthogroups with both species: {len(df_og):,}")
    print(f"  1-to-1 orthogroups: {len(df_og[(df_og.n_ecoli==1) & (df_og.n_human==1)]):,}")
    print(f"  Many-to-many: {len(df_og[(df_og.n_ecoli>1) | (df_og.n_human>1)]):,}")

    return df_og, reciprocal_pairs, rbh_pairs_from_best, pair_stats


# ──────────────────────────────────────────────────────────────────────
# STEP 3: Intersect with substrate lists
# ──────────────────────────────────────────────────────────────────────

def intersect_substrates(df_og):
    """Intersect orthogroups with GroEL and HSP60 substrate lists."""
    print("\n[5/6] Intersecting with substrate lists...")

    # Load substrate lists
    groel_df = pd.read_csv(GROEL_FILE, sep="\t")
    hsp60_df = pd.read_csv(HSP60_FILE, sep="\t")

    groel_ids = set(groel_df["current_accession"].dropna().unique())
    hsp60_ids = set(hsp60_df["uniprot_id"].dropna().unique())

    print(f"  GroEL substrates: {len(groel_ids)}")
    print(f"  HSP60 Tier 1 substrates: {len(hsp60_ids)}")

    # Build lookup: accession -> gene name and class for GroEL
    groel_gene_map = dict(zip(groel_df["current_accession"], groel_df["gene_symbol"]))
    groel_class_map = dict(zip(groel_df["current_accession"], groel_df["groel_class"]))

    # Build lookup for HSP60
    hsp60_gene_map = dict(zip(hsp60_df["uniprot_id"], hsp60_df["gene_name"]))

    # For each orthogroup, check membership
    results = []
    for _, row in df_og.iterrows():
        ec_ids = set(row["ecoli_accessions"].split(","))
        hs_ids = set(row["human_accessions"].split(","))

        groel_in = ec_ids & groel_ids
        hsp60_in = hs_ids & hsp60_ids

        if not groel_in and not hsp60_in:
            continue  # No substrates in this orthogroup

        if groel_in and hsp60_in:
            category = "shared"
        elif groel_in:
            category = "groel_only"
        else:
            category = "hsp60_only"

        groel_genes = ",".join([groel_gene_map.get(a, "?") for a in sorted(groel_in)]) if groel_in else ""
        groel_classes = ",".join([str(groel_class_map.get(a, "?")) for a in sorted(groel_in)]) if groel_in else ""
        hsp60_genes = ",".join([str(hsp60_gene_map.get(a, "?")) for a in sorted(hsp60_in)]) if hsp60_in else ""

        results.append({
            "orthogroup_id": row["orthogroup_id"],
            "groel_accessions": ",".join(sorted(groel_in)) if groel_in else "",
            "groel_genes": groel_genes,
            "groel_classes": groel_classes,
            "hsp60_accessions": ",".join(sorted(hsp60_in)) if hsp60_in else "",
            "hsp60_genes": hsp60_genes,
            "category": category,
            "best_evalue": row["best_evalue"],
            "mean_pident": row["mean_pident"],
            "n_ecoli": row["n_ecoli"],
            "n_human": row["n_human"],
            "ecoli_accessions": row["ecoli_accessions"],
            "human_accessions": row["human_accessions"]
        })

    df_sub = pd.DataFrame(results)

    n_shared = len(df_sub[df_sub.category == "shared"])
    n_groel_only = len(df_sub[df_sub.category == "groel_only"])
    n_hsp60_only = len(df_sub[df_sub.category == "hsp60_only"])

    print(f"\n  Substrate orthogroup summary:")
    print(f"    Shared (GroEL + HSP60):  {n_shared}")
    print(f"    GroEL-only:              {n_groel_only}")
    print(f"    HSP60-only:              {n_hsp60_only}")
    print(f"    Total substrate OGs:     {len(df_sub)}")

    return df_sub


# ──────────────────────────────────────────────────────────────────────
# STEP 4: Compare with RBH results
# ──────────────────────────────────────────────────────────────────────

def compare_with_rbh(df_og, reciprocal_pairs, df_sub):
    """Compare orthogroup results with the original 40 RBH pairs."""
    print("\n[6/6] Comparing with RBH results...")

    rbh_df = pd.read_csv(RBH_FILE, sep="\t")
    rbh_pairs = set(zip(rbh_df["groel_accession"], rbh_df["hsp60_accession"]))
    print(f"  Original RBH pairs: {len(rbh_pairs)}")

    # Check which RBH pairs are in the reciprocal pairs set
    rbh_in_recip = rbh_pairs & reciprocal_pairs
    print(f"  RBH pairs also in orthogroup reciprocal pairs: {len(rbh_in_recip)}")

    # Check which RBH pairs are in the same orthogroup
    # Build protein-to-orthogroup lookup
    ec_to_og = {}
    hs_to_og = {}
    for _, row in df_og.iterrows():
        og = row["orthogroup_id"]
        for ec in row["ecoli_accessions"].split(","):
            ec_to_og[ec] = og
        for hs in row["human_accessions"].split(","):
            hs_to_og[hs] = og

    comparison_rows = []
    for _, row in rbh_df.iterrows():
        ec = row["groel_accession"]
        hs = row["hsp60_accession"]
        ec_og = ec_to_og.get(ec, "")
        hs_og = hs_to_og.get(hs, "")
        in_same_og = ec_og == hs_og and ec_og != ""
        in_recip = (ec, hs) in reciprocal_pairs

        comparison_rows.append({
            "groel_accession": ec,
            "groel_gene": row.get("groel_gene", ""),
            "hsp60_accession": hs,
            "hsp60_gene": row.get("hsp60_gene", ""),
            "rbh_evalue": row.get("evalue", ""),
            "rbh_pident": row.get("pident", ""),
            "in_reciprocal_pairs": in_recip,
            "ecoli_orthogroup": ec_og,
            "human_orthogroup": hs_og,
            "in_same_orthogroup": in_same_og
        })

    df_comp = pd.DataFrame(comparison_rows)

    # Count additional substrate pairs from orthogroups (not in original RBH)
    shared_og = df_sub[df_sub.category == "shared"]
    og_substrate_pairs = set()
    for _, row in shared_og.iterrows():
        groel_accs = row["groel_accessions"].split(",") if row["groel_accessions"] else []
        hsp60_accs = row["hsp60_accessions"].split(",") if row["hsp60_accessions"] else []
        for ec in groel_accs:
            for hs in hsp60_accs:
                og_substrate_pairs.add((ec, hs))

    new_pairs = og_substrate_pairs - rbh_pairs
    print(f"  Substrate pairs in shared orthogroups: {len(og_substrate_pairs)}")
    print(f"  Additional pairs beyond RBH: {len(new_pairs)}")

    # Many-to-many relationships
    m2m_ogs = df_og[(df_og.n_ecoli > 1) | (df_og.n_human > 1)]
    m2m_with_substrates = 0
    for _, row in m2m_ogs.iterrows():
        ec_ids = set(row["ecoli_accessions"].split(","))
        hs_ids = set(row["human_accessions"].split(","))
        groel_df_ids = set(pd.read_csv(GROEL_FILE, sep="\t")["current_accession"].dropna())
        hsp60_df_ids = set(pd.read_csv(HSP60_FILE, sep="\t")["uniprot_id"].dropna())
        if (ec_ids & groel_df_ids) and (hs_ids & hsp60_df_ids):
            m2m_with_substrates += 1

    print(f"  Many-to-many OGs containing substrates from both: {m2m_with_substrates}")

    return df_comp, og_substrate_pairs, new_pairs


# ──────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Module C Step C4: Orthology Group Analysis")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    os.makedirs(OUT_DIR, exist_ok=True)

    # Step 1: Run MMseqs2 searches
    e2h_tsv, h2e_tsv = run_mmseqs_searches()

    # Step 2: Build orthogroups
    df_og, reciprocal_pairs, rbh_strict, pair_stats = build_orthogroups(e2h_tsv, h2e_tsv)

    # Save orthogroups
    og_file = os.path.join(OUT_DIR, "orthogroups_ecoli_human.tsv")
    df_og.to_csv(og_file, sep="\t", index=False)
    print(f"\n  Saved orthogroups: {og_file}")

    # Step 3: Intersect with substrates
    df_sub = intersect_substrates(df_og)

    sub_file = os.path.join(OUT_DIR, "substrate_orthogroups.tsv")
    df_sub.to_csv(sub_file, sep="\t", index=False)
    print(f"  Saved substrate orthogroups: {sub_file}")

    # Step 4: Compare with RBH
    df_comp, og_substrate_pairs, new_pairs = compare_with_rbh(df_og, reciprocal_pairs, df_sub)

    comp_file = os.path.join(OUT_DIR, "orthology_comparison.tsv")
    df_comp.to_csv(comp_file, sep="\t", index=False)
    print(f"  Saved comparison: {comp_file}")

    # Step 5: Generate summary report
    report = generate_report(df_og, df_sub, df_comp, reciprocal_pairs, rbh_strict,
                             og_substrate_pairs, new_pairs)

    report_file = os.path.join(OUT_DIR, "orthology_summary_report.txt")
    with open(report_file, "w") as f:
        f.write(report)
    print(f"  Saved report: {report_file}")

    # Cleanup work directory (keep TSV outputs)
    print(f"\n  (Work directory kept at: {WORK_DIR})")
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)


def generate_report(df_og, df_sub, df_comp, reciprocal_pairs, rbh_strict,
                    og_substrate_pairs, new_pairs):
    """Generate a comprehensive text summary report."""

    n_shared = len(df_sub[df_sub.category == "shared"])
    n_groel_only = len(df_sub[df_sub.category == "groel_only"])
    n_hsp60_only = len(df_sub[df_sub.category == "hsp60_only"])

    rbh_in_same_og = df_comp["in_same_orthogroup"].sum()

    # Class breakdown for shared
    shared = df_sub[df_sub.category == "shared"]
    class_counts = defaultdict(int)
    for _, row in shared.iterrows():
        classes = row["groel_classes"].split(",") if row["groel_classes"] else []
        for c in classes:
            c = c.strip()
            if c:
                class_counts[c] += 1

    lines = []
    lines.append("=" * 70)
    lines.append("MODULE C - STEP C4: ORTHOLOGY GROUP ANALYSIS REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 70)

    lines.append("\n1. METHOD")
    lines.append("-" * 40)
    lines.append("Approach: MMseqs2 all-vs-all reciprocal search + connected-component clustering")
    lines.append("  - Bidirectional MMseqs2 search (E. coli <-> Human)")
    lines.append(f"  - E-value cutoff: {EVALUE_CUTOFF}")
    lines.append(f"  - Min percent identity: {MIN_PIDENT}%")
    lines.append(f"  - Min coverage (query or target): {MIN_COVERAGE*100}%")
    lines.append("  - Reciprocal hit filtering: pair must be found in both search directions")
    lines.append("  - Orthogroup clustering: union-find connected components on reciprocal pairs")
    lines.append("  - This captures many-to-many orthology relationships (paralogs + orthologs)")

    lines.append("\n2. ORTHOGROUP STATISTICS")
    lines.append("-" * 40)
    lines.append(f"Total orthogroups (both species): {len(df_og):,}")
    n_11 = len(df_og[(df_og.n_ecoli == 1) & (df_og.n_human == 1)])
    n_1m = len(df_og[(df_og.n_ecoli == 1) & (df_og.n_human > 1)])
    n_m1 = len(df_og[(df_og.n_ecoli > 1) & (df_og.n_human == 1)])
    n_mm = len(df_og[(df_og.n_ecoli > 1) & (df_og.n_human > 1)])
    lines.append(f"  1-to-1 orthogroups:    {n_11:,}")
    lines.append(f"  1-to-many (EC:HS):     {n_1m:,}")
    lines.append(f"  many-to-1 (EC:HS):     {n_m1:,}")
    lines.append(f"  many-to-many:          {n_mm:,}")
    lines.append(f"Total reciprocal pairs:  {len(reciprocal_pairs):,}")
    lines.append(f"Strict RBH pairs:        {len(rbh_strict):,}")

    # Size distribution
    sizes = df_og["n_ecoli"] + df_og["n_human"]
    lines.append(f"\nOrthogroup size distribution:")
    lines.append(f"  Min:    {sizes.min()}")
    lines.append(f"  Median: {sizes.median():.0f}")
    lines.append(f"  Mean:   {sizes.mean():.1f}")
    lines.append(f"  Max:    {sizes.max()}")

    # Identity distribution
    lines.append(f"\nMean percent identity across orthogroups:")
    lines.append(f"  Min:    {df_og.mean_pident.min():.1f}%")
    lines.append(f"  Median: {df_og.mean_pident.median():.1f}%")
    lines.append(f"  Mean:   {df_og.mean_pident.mean():.1f}%")
    lines.append(f"  Max:    {df_og.mean_pident.max():.1f}%")

    lines.append("\n3. SUBSTRATE ORTHOGROUP ANALYSIS")
    lines.append("-" * 40)
    lines.append(f"Orthogroups containing chaperonin substrates: {len(df_sub):,}")
    lines.append(f"  SHARED (GroEL + HSP60 substrate):  {n_shared}")
    lines.append(f"  GroEL-only (no HSP60 substrate):   {n_groel_only}")
    lines.append(f"  HSP60-only (no GroEL substrate):   {n_hsp60_only}")

    lines.append(f"\nGroEL class breakdown in SHARED orthogroups:")
    for cls in sorted(class_counts.keys()):
        lines.append(f"  Class {cls}: {class_counts[cls]}")

    lines.append(f"\nSubstrate pairs in shared orthogroups: {len(og_substrate_pairs)}")

    lines.append("\n4. COMPARISON WITH SIMPLE RBH (40 PAIRS)")
    lines.append("-" * 40)
    lines.append(f"Original RBH substrate pairs: 40")
    lines.append(f"RBH pairs in same orthogroup: {rbh_in_same_og}")
    lines.append(f"Additional substrate pairs from orthogroups: {len(new_pairs)}")
    lines.append(f"Total substrate pairs (orthogroup method): {len(og_substrate_pairs)}")

    if new_pairs:
        lines.append(f"\nNew substrate pairs found (not in RBH):")
        for ec, hs in sorted(new_pairs):
            lines.append(f"  {ec} -- {hs}")

    lines.append("\n5. TOP SHARED ORTHOGROUPS")
    lines.append("-" * 40)
    shared_sorted = shared.sort_values("best_evalue")
    for i, (_, row) in enumerate(shared_sorted.head(20).iterrows()):
        lines.append(
            f"  {row['orthogroup_id']}: "
            f"GroEL=[{row['groel_genes']}] (class {row['groel_classes']}) <-> "
            f"HSP60=[{row['hsp60_genes']}] | "
            f"e={row['best_evalue']:.1e}, pident={row['mean_pident']:.1f}%"
        )

    lines.append("\n6. INTERPRETATION")
    lines.append("-" * 40)
    lines.append(f"The orthogroup approach identified {n_shared} shared chaperonin-client orthogroups,")
    lines.append(f"containing {len(og_substrate_pairs)} GroEL-HSP60 substrate pairs.")
    if len(og_substrate_pairs) > 40:
        lines.append(f"This is {len(og_substrate_pairs) - 40} more than the 40 pairs from simple RBH,")
        lines.append("due to many-to-many orthology relationships (paralogs with shared function).")
    lines.append(f"\n{n_groel_only} GroEL substrates have human orthologs but the ortholog is NOT")
    lines.append("a known HSP60 substrate - these represent potential cases of:")
    lines.append("  - Lost chaperonin dependence in the human lineage")
    lines.append("  - Undetected HSP60 substrates (experimental coverage gap)")
    lines.append(f"\n{n_hsp60_only} HSP60 substrates have E. coli orthologs but the ortholog is NOT")
    lines.append("a known GroEL substrate - potential cases of:")
    lines.append("  - Gained chaperonin dependence in the eukaryotic/mitochondrial lineage")
    lines.append("  - Undetected GroEL substrates")

    lines.append("\n" + "=" * 70)
    return "\n".join(lines)


if __name__ == "__main__":
    main()
