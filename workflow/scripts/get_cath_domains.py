#!/usr/bin/env python3
"""
Module E Step E1: Obtain CATH structural domain assignments for all pilot proteins.

Uses the InterPro/Gene3D API to retrieve CATH domain assignments for each
UniProt accession in the pilot set. Falls back gracefully on failures.

Usage:
    python3 get_cath_domains.py
"""

import os
import sys
import time
import json
import requests
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
STRUCTURE_INDEX = BASE / "results" / "structures" / "structure_index.tsv"
OUT_DIR = BASE / "results" / "domains"
OUT_DOMAINS = OUT_DIR / "cath_domain_assignments.tsv"
OUT_SUMMARY = OUT_DIR / "cath_protein_summary.tsv"

# InterPro Gene3D API endpoint
INTERPRO_URL = (
    "https://www.ebi.ac.uk/interpro/api/entry/cathgene3d/"
    "protein/uniprot/{accession}?format=json"
)

# Rate limiting
REQUESTS_PER_SECOND = 1.0  # conservative
MIN_DELAY = 1.0 / REQUESTS_PER_SECOND

# Retry config
MAX_RETRIES = 3
RETRY_BACKOFF = 5  # seconds

# Checkpoint file (for resume on interruption)
CHECKPOINT = OUT_DIR / "_cath_checkpoint.json"

# ── Helpers ──────────────────────────────────────────────────────────────────

def parse_cath_code(accession_str: str) -> str:
    """Extract CATH superfamily code from Gene3D accession like G3DSA:3.40.50.300"""
    if accession_str.startswith("G3DSA:"):
        return accession_str[6:]
    return accession_str


def cath_levels(sf_code: str) -> dict:
    """Parse CATH superfamily code into class/architecture/topology."""
    parts = sf_code.split(".")
    return {
        "cath_class": parts[0] if len(parts) >= 1 else "",
        "cath_architecture": ".".join(parts[:2]) if len(parts) >= 2 else "",
        "cath_topology": ".".join(parts[:3]) if len(parts) >= 3 else "",
    }


def query_interpro(accession: str, session: requests.Session) -> dict:
    """Query InterPro Gene3D API for a single UniProt accession.
    Returns raw JSON dict or None on failure.
    """
    url = INTERPRO_URL.format(accession=accession)
    for attempt in range(MAX_RETRIES):
        try:
            resp = session.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 204:
                # No content — protein exists but no CATH assignments
                return {"count": 0, "results": []}
            elif resp.status_code == 404:
                # Protein not found in InterPro
                return None
            elif resp.status_code == 408 or resp.status_code == 429:
                # Timeout or rate-limited — back off
                wait = RETRY_BACKOFF * (attempt + 1)
                print(f"    Rate-limited/timeout for {accession}, waiting {wait}s...")
                time.sleep(wait)
                continue
            elif resp.status_code >= 500:
                wait = RETRY_BACKOFF * (attempt + 1)
                print(f"    Server error {resp.status_code} for {accession}, retrying in {wait}s...")
                time.sleep(wait)
                continue
            else:
                print(f"    Unexpected status {resp.status_code} for {accession}")
                return None
        except requests.exceptions.RequestException as e:
            wait = RETRY_BACKOFF * (attempt + 1)
            print(f"    Request exception for {accession}: {e}, retrying in {wait}s...")
            time.sleep(wait)
    print(f"    FAILED after {MAX_RETRIES} retries: {accession}")
    return None


def parse_domains(accession: str, data: dict) -> list:
    """Parse InterPro Gene3D JSON into a list of domain dicts."""
    domains = []
    if not data or data.get("count", 0) == 0:
        return domains

    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        raw_sf = meta.get("accession", "")
        sf_code = parse_cath_code(raw_sf)
        sf_name = meta.get("name", "")
        levels = cath_levels(sf_code)

        proteins = entry.get("proteins", [])
        if not proteins:
            continue

        # There should be exactly one protein entry (the one we queried)
        for prot in proteins:
            protein_length = prot.get("protein_length", 0)
            for loc in prot.get("entry_protein_locations", []):
                fragments = loc.get("fragments", [])
                if not fragments:
                    continue

                # A domain can have discontinuous fragments
                # We record the overall span and individual fragments
                starts = [f["start"] for f in fragments]
                ends = [f["end"] for f in fragments]
                overall_start = min(starts)
                overall_end = max(ends)
                domain_length = sum(f["end"] - f["start"] + 1 for f in fragments)

                # Build fragment string for reference
                frag_str = ";".join(f"{f['start']}-{f['end']}" for f in fragments)

                domains.append({
                    "uniprot_accession": accession,
                    "cath_superfamily": sf_code,
                    "cath_superfamily_name": sf_name,
                    "cath_class": levels["cath_class"],
                    "cath_architecture": levels["cath_architecture"],
                    "cath_topology": levels["cath_topology"],
                    "domain_start": overall_start,
                    "domain_end": overall_end,
                    "domain_length": domain_length,
                    "fragments": frag_str,
                    "n_fragments": len(fragments),
                    "protein_length": protein_length,
                    "source": "InterPro_Gene3D",
                })
    return domains


def load_checkpoint() -> dict:
    """Load checkpoint data if it exists."""
    if CHECKPOINT.exists():
        with open(CHECKPOINT) as f:
            return json.load(f)
    return {}


def save_checkpoint(data: dict):
    """Save checkpoint data."""
    with open(CHECKPOINT, "w") as f:
        json.dump(data, f)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load structure index
    print("Loading structure index...")
    si = pd.read_csv(STRUCTURE_INDEX, sep="\t")
    accessions = si["uniprot_accession"].unique().tolist()
    print(f"  Total unique accessions: {len(accessions)}")

    # Load checkpoint (already-queried proteins)
    checkpoint = load_checkpoint()
    cached_domains = checkpoint.get("domains", {})  # accession -> list of domain dicts
    cached_no_assignment = set(checkpoint.get("no_assignment", []))
    cached_failed = set(checkpoint.get("failed", []))

    already_done = set(cached_domains.keys()) | cached_no_assignment | cached_failed
    remaining = [a for a in accessions if a not in already_done]
    print(f"  Already cached: {len(already_done)} proteins")
    print(f"  Remaining to query: {len(remaining)} proteins")

    # Query InterPro for remaining proteins
    session = requests.Session()
    session.headers.update({
        "Accept": "application/json",
        "User-Agent": "AntahAstiPrarambh/1.0 (research project; CATH domain retrieval)"
    })

    new_domains = {}
    new_no_assignment = set()
    new_failed = set()

    total = len(remaining)
    checkpoint_interval = 50  # save checkpoint every N proteins
    last_request_time = 0

    for i, acc in enumerate(remaining, 1):
        # Rate limiting
        elapsed = time.time() - last_request_time
        if elapsed < MIN_DELAY:
            time.sleep(MIN_DELAY - elapsed)

        if i % 50 == 0 or i == 1 or i == total:
            print(f"  Querying {i}/{total}: {acc} ...")

        last_request_time = time.time()
        data = query_interpro(acc, session)

        if data is None:
            new_failed.add(acc)
        elif data.get("count", 0) == 0:
            new_no_assignment.add(acc)
        else:
            domains = parse_domains(acc, data)
            if domains:
                new_domains[acc] = domains
            else:
                new_no_assignment.add(acc)

        # Periodic checkpoint
        if i % checkpoint_interval == 0:
            merged_domains = {**cached_domains, **new_domains}
            merged_no = list(cached_no_assignment | new_no_assignment)
            merged_fail = list(cached_failed | new_failed)
            # Convert domain lists for JSON serialisation
            save_checkpoint({
                "domains": {k: v for k, v in merged_domains.items()},
                "no_assignment": merged_no,
                "failed": merged_fail,
            })
            if i % 200 == 0:
                n_with = len(cached_domains) + len(new_domains)
                n_without = len(cached_no_assignment) + len(new_no_assignment)
                n_fail = len(cached_failed) + len(new_failed)
                print(f"    [Checkpoint] with={n_with}, without={n_without}, failed={n_fail}")

    # Merge all results
    all_domains_dict = {**cached_domains, **new_domains}
    all_no_assignment = cached_no_assignment | new_no_assignment
    all_failed = cached_failed | new_failed

    # Final checkpoint
    save_checkpoint({
        "domains": all_domains_dict,
        "no_assignment": list(all_no_assignment),
        "failed": list(all_failed),
    })

    print(f"\nQuery complete:")
    print(f"  Proteins with CATH assignments: {len(all_domains_dict)}")
    print(f"  Proteins without CATH assignments: {len(all_no_assignment)}")
    print(f"  Failed queries: {len(all_failed)}")

    # ── Build domain table ───────────────────────────────────────────────────
    print("\nBuilding domain assignment table...")

    all_domain_rows = []
    for acc, dom_list in all_domains_dict.items():
        for idx, d in enumerate(dom_list, 1):
            row = dict(d)
            row["domain_index"] = idx
            all_domain_rows.append(row)

    if all_domain_rows:
        df_domains = pd.DataFrame(all_domain_rows)
        # Reorder columns
        col_order = [
            "uniprot_accession", "domain_index", "cath_superfamily",
            "cath_superfamily_name", "cath_class", "cath_architecture",
            "cath_topology", "domain_start", "domain_end", "domain_length",
            "fragments", "n_fragments", "protein_length", "source",
        ]
        df_domains = df_domains[[c for c in col_order if c in df_domains.columns]]
        df_domains.sort_values(["uniprot_accession", "domain_index"], inplace=True)
    else:
        df_domains = pd.DataFrame(columns=[
            "uniprot_accession", "domain_index", "cath_superfamily",
            "cath_superfamily_name", "cath_class", "cath_architecture",
            "cath_topology", "domain_start", "domain_end", "domain_length",
            "fragments", "n_fragments", "protein_length", "source",
        ])

    df_domains.to_csv(OUT_DOMAINS, sep="\t", index=False)
    print(f"  Saved: {OUT_DOMAINS}  ({len(df_domains)} domain rows)")

    # ── Build protein summary table ──────────────────────────────────────────
    print("Building protein summary table...")

    summary_rows = []
    for acc in accessions:
        if acc in all_domains_dict and all_domains_dict[acc]:
            doms = all_domains_dict[acc]
            n_doms = len(doms)
            # Domain architecture: ordered by start position
            sorted_doms = sorted(doms, key=lambda d: d["domain_start"])
            arch = "|".join(d["cath_superfamily"] for d in sorted_doms)
            prot_len = doms[0].get("protein_length", 0)
            has_cath = True
        else:
            n_doms = 0
            arch = ""
            prot_len = 0
            has_cath = False

        # Try to get protein length from structure index if missing
        if prot_len == 0:
            si_row = si[si["uniprot_accession"] == acc]
            if not si_row.empty and "residues_modeled" in si_row.columns:
                val = si_row.iloc[0]["residues_modeled"]
                if pd.notna(val):
                    prot_len = int(val)

        summary_rows.append({
            "uniprot_accession": acc,
            "n_domains": n_doms,
            "domain_architecture": arch,
            "has_cath_assignment": has_cath,
            "protein_length": prot_len,
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"  Saved: {OUT_SUMMARY}  ({len(df_summary)} protein rows)")

    # ── Summary report ───────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("CATH DOMAIN ASSIGNMENT SUMMARY REPORT")
    print("=" * 72)

    n_total = len(accessions)
    n_with = df_summary["has_cath_assignment"].sum()
    n_without = n_total - n_with
    pct = 100 * n_with / n_total if n_total > 0 else 0

    print(f"\nCoverage:")
    print(f"  Total pilot proteins:               {n_total}")
    print(f"  With CATH assignment:                {n_with} ({pct:.1f}%)")
    print(f"  Without CATH assignment:             {n_without} ({100 - pct:.1f}%)")
    print(f"    (of which failed queries:          {len(all_failed)})")

    # Domain count distribution
    print(f"\nDomain count distribution:")
    dom_counts = df_summary[df_summary["has_cath_assignment"]]["n_domains"]
    for nd in sorted(dom_counts.unique()):
        cnt = (dom_counts == nd).sum()
        print(f"  {nd} domain(s): {cnt} proteins ({100 * cnt / n_with:.1f}%)" if n_with > 0 else f"  {nd} domain(s): {cnt}")

    mean_dom = dom_counts.mean() if len(dom_counts) > 0 else 0
    print(f"  Mean domains per protein: {mean_dom:.2f}")

    # Most common superfamilies
    if not df_domains.empty:
        print(f"\nTop 20 most common CATH superfamilies:")
        sf_counts = df_domains.groupby(["cath_superfamily", "cath_superfamily_name"]).size()
        sf_counts = sf_counts.sort_values(ascending=False)
        for (sf, name), cnt in sf_counts.head(20).items():
            print(f"  {sf:20s} {cnt:5d}  {name}")

    # CATH class distribution
    if not df_domains.empty:
        print(f"\nCATH class distribution:")
        class_map = {"1": "Mainly-alpha", "2": "Mainly-beta", "3": "Alpha-Beta", "4": "Few-SS"}
        class_counts = df_domains["cath_class"].value_counts().sort_index()
        for cls, cnt in class_counts.items():
            label = class_map.get(str(cls), "Unknown")
            print(f"  Class {cls} ({label}): {cnt} domains")

    # Breakdown by dataset
    print(f"\nBreakdown by source dataset:")
    # Merge summary with structure index to get dataset info
    si_datasets = si[["uniprot_accession", "source_dataset"]].copy()
    merged = df_summary.merge(si_datasets, on="uniprot_accession", how="left")

    # Expand multi-dataset entries
    dataset_rows = []
    for _, row in merged.iterrows():
        datasets = str(row["source_dataset"]).split(",") if pd.notna(row["source_dataset"]) else ["unknown"]
        for ds in datasets:
            dataset_rows.append({
                "dataset": ds.strip(),
                "has_cath": row["has_cath_assignment"],
                "n_domains": row["n_domains"],
            })
    df_ds = pd.DataFrame(dataset_rows)

    for ds in sorted(df_ds["dataset"].unique()):
        subset = df_ds[df_ds["dataset"] == ds]
        total_ds = len(subset)
        with_cath = subset["has_cath"].sum()
        pct_ds = 100 * with_cath / total_ds if total_ds > 0 else 0
        mean_d = subset[subset["has_cath"]]["n_domains"].mean() if with_cath > 0 else 0
        print(f"  {ds:20s}  total={total_ds:5d}  with_CATH={with_cath:5d} ({pct_ds:.1f}%)  mean_domains={mean_d:.2f}")

    # Multi-domain vs single-domain
    if n_with > 0:
        single = (dom_counts == 1).sum()
        multi = (dom_counts > 1).sum()
        print(f"\nSingle-domain vs multi-domain:")
        print(f"  Single-domain: {single} ({100 * single / n_with:.1f}%)")
        print(f"  Multi-domain:  {multi} ({100 * multi / n_with:.1f}%)")

    print("\n" + "=" * 72)
    print("Done. Proteins without CATH assignments will need Chainsaw/Merizo (Step E2).")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    sys.exit(main())
