#!/usr/bin/env python3
"""
Analyze Foldseek structural clustering results for pilot AlphaFold structures.

Module E, Step E4: Parse cluster membership, compute statistics, and generate
per-dataset breakdowns for GroEL substrates, HSP60 interactors,
mitochondrial proteome, and matrix-localised proteins.

Inputs:
    - Foldseek cluster_membership.tsv  (representative  member)
    - structure_index.tsv              (accession -> source_dataset mapping)

Outputs:
    - results/domains/foldseek_clusters.tsv        (per-protein cluster table)
    - results/domains/foldseek_cluster_summary.txt  (human-readable report)
"""

import os
import sys
import re
from collections import defaultdict, Counter

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
PROJECT = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"
CLUSTER_TSV = os.path.join(PROJECT, "results/domains/foldseek/cluster_membership.tsv")
SEARCH_TSV = os.path.join(PROJECT, "results/domains/foldseek/search_results.tsv")
STRUCT_INDEX = os.path.join(PROJECT, "results/structures/structure_index.tsv")

OUT_CLUSTERS = os.path.join(PROJECT, "results/domains/foldseek_clusters.tsv")
OUT_SUMMARY = os.path.join(PROJECT, "results/domains/foldseek_cluster_summary.txt")


def extract_accession(cif_name: str) -> str:
    """Extract UniProt accession from CIF filename stem.

    AF-{ACC}-F1-model_v6 -> ACC
    """
    m = re.match(r"AF-([A-Z0-9]+)-F1-model_v\d+", cif_name)
    if m:
        return m.group(1)
    return cif_name


def load_cluster_membership(path: str):
    """Parse Foldseek createtsv output (representative  member).

    Returns:
        member_to_rep: dict  member_accession -> representative_accession
        rep_to_members: dict  representative_accession -> [member_accessions]
    """
    member_to_rep = {}
    rep_to_members = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            rep_name = parts[0]
            member_name = parts[1]

            rep_acc = extract_accession(rep_name)
            member_acc = extract_accession(member_name)

            member_to_rep[member_acc] = rep_acc
            rep_to_members[rep_acc].append(member_acc)

    return member_to_rep, rep_to_members


def load_source_datasets(path: str):
    """Load structure_index.tsv and return accession -> set of source datasets."""
    df = pd.read_csv(path, sep="\t")
    acc_to_sources = {}
    for _, row in df.iterrows():
        acc = str(row["uniprot_accession"]).strip()
        sources = str(row.get("source_dataset", "")).strip()
        if sources and sources != "nan":
            acc_to_sources[acc] = set(sources.split(","))
        else:
            acc_to_sources[acc] = set()
    return acc_to_sources


def main():
    print("=" * 70)
    print("  MODULE E / STEP E4 — FOLDSEEK STRUCTURAL CLUSTERING ANALYSIS")
    print("=" * 70)

    # ── Load data ─────────────────────────────────────────────────────────
    print("\nLoading cluster membership...")
    member_to_rep, rep_to_members = load_cluster_membership(CLUSTER_TSV)
    print(f"  Total proteins assigned: {len(member_to_rep)}")
    print(f"  Total clusters:          {len(rep_to_members)}")

    print("\nLoading source dataset mapping...")
    acc_to_sources = load_source_datasets(STRUCT_INDEX)
    print(f"  Proteins in structure index: {len(acc_to_sources)}")

    # ── Assign numeric cluster IDs ────────────────────────────────────────
    # Sort representatives by cluster size (largest first) for stable IDs
    sorted_reps = sorted(rep_to_members.keys(),
                         key=lambda r: (-len(rep_to_members[r]), r))
    rep_to_cluid = {rep: i + 1 for i, rep in enumerate(sorted_reps)}

    # ── Compute cluster sizes ─────────────────────────────────────────────
    cluster_sizes = {rep: len(members) for rep, members in rep_to_members.items()}

    # Size distribution
    size_counter = Counter(cluster_sizes.values())

    singletons = sum(1 for s in cluster_sizes.values() if s == 1)
    small = sum(1 for s in cluster_sizes.values() if 2 <= s <= 5)
    medium = sum(1 for s in cluster_sizes.values() if 6 <= s <= 20)
    large = sum(1 for s in cluster_sizes.values() if s > 20)

    # ── Build per-protein output table ────────────────────────────────────
    rows = []
    for member_acc, rep_acc in sorted(member_to_rep.items()):
        cluid = rep_to_cluid[rep_acc]
        csize = cluster_sizes[rep_acc]
        sources = acc_to_sources.get(member_acc, set())
        source_str = ",".join(sorted(sources)) if sources else "unknown"
        rows.append({
            "uniprot_accession": member_acc,
            "cluster_id": cluid,
            "cluster_representative": rep_acc,
            "cluster_size": csize,
            "source_dataset": source_str,
        })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUT_CLUSTERS, sep="\t", index=False)
    print(f"\nCluster table written: {OUT_CLUSTERS}")
    print(f"  Rows: {len(df_out)}")

    # ── Per-dataset analysis ──────────────────────────────────────────────
    datasets = ["groel", "hsp60", "mito", "matrix"]
    dataset_proteins = {}
    dataset_clusters = {}

    for ds in datasets:
        proteins = set()
        for acc, sources in acc_to_sources.items():
            if ds in sources and acc in member_to_rep:
                proteins.add(acc)
        dataset_proteins[ds] = proteins

        # Which clusters are they in?
        cluster_ids = set()
        for acc in proteins:
            rep = member_to_rep[acc]
            cluster_ids.add(rep_to_cluid[rep])
        dataset_clusters[ds] = cluster_ids

    # ── Shared clusters between GroEL and HSP60 ──────────────────────────
    groel_clusters = dataset_clusters.get("groel", set())
    hsp60_clusters = dataset_clusters.get("hsp60", set())
    shared_clusters = groel_clusters & hsp60_clusters

    # For shared clusters, find which proteins belong
    shared_cluster_details = {}
    for cluid in sorted(shared_clusters):
        rep = sorted_reps[cluid - 1]  # cluid is 1-based
        members = rep_to_members[rep]
        groel_members = [m for m in members if "groel" in acc_to_sources.get(m, set())]
        hsp60_members = [m for m in members if "hsp60" in acc_to_sources.get(m, set())]
        if groel_members and hsp60_members:
            shared_cluster_details[cluid] = {
                "representative": rep,
                "size": cluster_sizes[rep],
                "groel_members": groel_members,
                "hsp60_members": hsp60_members,
            }

    # ── Top clusters per dataset ──────────────────────────────────────────
    def top_clusters_for_dataset(ds, n=10):
        """Return top-N clusters by number of proteins from this dataset."""
        cluster_counts = Counter()
        for acc in dataset_proteins[ds]:
            rep = member_to_rep[acc]
            cluid = rep_to_cluid[rep]
            cluster_counts[cluid] += 1
        return cluster_counts.most_common(n)

    # ── Load search results for stats ─────────────────────────────────────
    search_stats = None
    if os.path.exists(SEARCH_TSV):
        try:
            search_df = pd.read_csv(SEARCH_TSV, sep="\t", header=None,
                                     names=["query", "target", "evalue", "bits",
                                            "pident", "alnlen", "qcov", "tcov"])
            # Remove self-hits
            search_df["query_acc"] = search_df["query"].apply(extract_accession)
            search_df["target_acc"] = search_df["target"].apply(extract_accession)
            non_self = search_df[search_df["query_acc"] != search_df["target_acc"]]
            search_stats = {
                "total_hits": len(search_df),
                "non_self_hits": len(non_self),
                "proteins_with_hits": non_self["query_acc"].nunique(),
                "mean_pident": non_self["pident"].mean(),
                "median_pident": non_self["pident"].median(),
                "mean_evalue": non_self["evalue"].mean(),
                "mean_qcov": non_self["qcov"].mean(),
                "mean_tcov": non_self["tcov"].mean(),
            }
        except Exception as e:
            print(f"  Warning: could not parse search results: {e}")

    # ── Generate summary report ───────────────────────────────────────────
    lines = []
    lines.append("=" * 70)
    lines.append("  FOLDSEEK STRUCTURAL CLUSTERING — SUMMARY REPORT")
    lines.append("  Module E, Step E4")
    lines.append("=" * 70)
    lines.append("")
    lines.append("CLUSTERING PARAMETERS:")
    lines.append("  Method:               Foldseek cluster (cascaded, set-cover)")
    lines.append("  Min sequence identity: 0.3 (30%)")
    lines.append("  Coverage threshold:    0.5 (50%)")
    lines.append("  E-value threshold:     0.01 (cluster), 0.001 (search)")
    lines.append("")

    lines.append("OVERALL STATISTICS:")
    lines.append(f"  Total proteins:        {len(member_to_rep)}")
    lines.append(f"  Total clusters:        {len(rep_to_members)}")
    lines.append(f"  Singletons (size=1):   {singletons} ({100*singletons/len(rep_to_members):.1f}%)")
    lines.append(f"  Small (2-5):           {small} ({100*small/len(rep_to_members):.1f}%)")
    lines.append(f"  Medium (6-20):         {medium} ({100*medium/len(rep_to_members):.1f}%)")
    lines.append(f"  Large (>20):           {large} ({100*large/len(rep_to_members):.1f}%)")
    lines.append("")

    # Size distribution table
    lines.append("CLUSTER SIZE DISTRIBUTION:")
    lines.append(f"  {'Size':>6s}  {'Count':>6s}  {'Proteins':>10s}")
    lines.append(f"  {'-'*6}  {'-'*6}  {'-'*10}")
    for size in sorted(size_counter.keys()):
        count = size_counter[size]
        proteins = size * count
        lines.append(f"  {size:>6d}  {count:>6d}  {proteins:>10d}")
    lines.append("")

    # Largest clusters
    lines.append("TOP 20 LARGEST CLUSTERS:")
    lines.append(f"  {'Rank':>4s}  {'ClusterID':>9s}  {'Size':>5s}  {'Representative':>15s}  Members")
    lines.append(f"  {'-'*4}  {'-'*9}  {'-'*5}  {'-'*15}  {'-'*40}")
    for rank, rep in enumerate(sorted_reps[:20], 1):
        cluid = rep_to_cluid[rep]
        csize = cluster_sizes[rep]
        members = rep_to_members[rep]
        member_str = ", ".join(sorted(members)[:8])
        if len(members) > 8:
            member_str += f" ... (+{len(members)-8} more)"
        lines.append(f"  {rank:>4d}  {cluid:>9d}  {csize:>5d}  {rep:>15s}  {member_str}")
    lines.append("")

    # Search results stats
    if search_stats:
        lines.append("ALL-VS-ALL SEARCH STATISTICS:")
        lines.append(f"  Total hits (incl. self): {search_stats['total_hits']}")
        lines.append(f"  Non-self hits:           {search_stats['non_self_hits']}")
        lines.append(f"  Proteins with hits:      {search_stats['proteins_with_hits']}")
        lines.append(f"  Mean % identity:         {search_stats['mean_pident']:.1f}%")
        lines.append(f"  Median % identity:       {search_stats['median_pident']:.1f}%")
        lines.append(f"  Mean query coverage:     {search_stats['mean_qcov']:.3f}")
        lines.append(f"  Mean target coverage:    {search_stats['mean_tcov']:.3f}")
        lines.append("")

    # Per-dataset breakdown
    lines.append("PER-DATASET BREAKDOWN:")
    lines.append(f"  {'Dataset':>10s}  {'Proteins':>10s}  {'Clusters':>10s}  {'Coverage':>10s}")
    lines.append(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
    for ds in datasets:
        n_prot = len(dataset_proteins[ds])
        n_clust = len(dataset_clusters[ds])
        cov = f"{100*n_clust/len(rep_to_members):.1f}%" if len(rep_to_members) > 0 else "N/A"
        lines.append(f"  {ds:>10s}  {n_prot:>10d}  {n_clust:>10d}  {cov:>10s}")
    lines.append("")

    # Top clusters per dataset
    for ds in datasets:
        lines.append(f"  TOP 10 CLUSTERS FOR {ds.upper()}:")
        top = top_clusters_for_dataset(ds, n=10)
        if not top:
            lines.append("    (no proteins in this dataset)")
        else:
            lines.append(f"    {'ClusterID':>9s}  {'#Proteins':>10s}  {'ClusterSize':>11s}  Representative")
            lines.append(f"    {'-'*9}  {'-'*10}  {'-'*11}  {'-'*15}")
            for cluid, count in top:
                rep = sorted_reps[cluid - 1]
                csize = cluster_sizes[rep]
                lines.append(f"    {cluid:>9d}  {count:>10d}  {csize:>11d}  {rep}")
        lines.append("")

    # Shared clusters between GroEL and HSP60
    lines.append("SHARED STRUCTURAL CLUSTERS: GroEL substrates & HSP60 interactors")
    lines.append(f"  GroEL-only clusters:   {len(groel_clusters - hsp60_clusters)}")
    lines.append(f"  HSP60-only clusters:   {len(hsp60_clusters - groel_clusters)}")
    lines.append(f"  Shared clusters:       {len(shared_clusters)}")
    lines.append("")

    if shared_cluster_details:
        lines.append("  SHARED CLUSTER DETAILS (clusters containing both GroEL and HSP60 members):")
        lines.append(f"    {'ClusterID':>9s}  {'Size':>5s}  {'#GroEL':>7s}  {'#HSP60':>7s}  Representative")
        lines.append(f"    {'-'*9}  {'-'*5}  {'-'*7}  {'-'*7}  {'-'*15}")
        for cluid in sorted(shared_cluster_details.keys()):
            info = shared_cluster_details[cluid]
            lines.append(f"    {cluid:>9d}  {info['size']:>5d}  "
                         f"{len(info['groel_members']):>7d}  "
                         f"{len(info['hsp60_members']):>7d}  "
                         f"{info['representative']}")
        lines.append("")

        # Show member proteins for top shared clusters
        sorted_shared = sorted(shared_cluster_details.items(),
                               key=lambda x: -(len(x[1]['groel_members']) + len(x[1]['hsp60_members'])))
        for cluid, info in sorted_shared[:10]:
            lines.append(f"    Cluster {cluid} (rep={info['representative']}, size={info['size']}):")
            lines.append(f"      GroEL members: {', '.join(sorted(info['groel_members']))}")
            lines.append(f"      HSP60 members: {', '.join(sorted(info['hsp60_members']))}")
            lines.append("")
    else:
        lines.append("  No clusters found with members from both GroEL and HSP60 datasets.")
        lines.append("  (This may indicate the two datasets share fold similarities at a")
        lines.append("   broader structural level not captured by strict clustering.)")
        lines.append("")

    lines.append("=" * 70)

    report = "\n".join(lines)
    with open(OUT_SUMMARY, "w") as fh:
        fh.write(report + "\n")
    print(f"\nSummary report written: {OUT_SUMMARY}")
    print("\n" + report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
