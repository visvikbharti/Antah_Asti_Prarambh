#!/usr/bin/env python3
"""
Module C - Steps C1 & C2: Extract FASTA sequences for GroEL and HSP60 substrates.
"""

import pandas as pd
import requests
import time
import sys
from Bio import SeqIO
from collections import OrderedDict

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"

GROEL_TSV   = f"{BASE}/data/processed/groel_substrates_standardized.tsv"
HSP60_TSV   = f"{BASE}/data/processed/hsp60_tier1_substrates.tsv"
ECOLI_FASTA = f"{BASE}/data/raw/uniprot/ecoli_k12_proteome.fasta"
HUMAN_FASTA = f"{BASE}/data/raw/uniprot/human_proteome.fasta"

GROEL_OUT   = f"{BASE}/data/processed/groel_substrates.fasta"
HSP60_OUT   = f"{BASE}/data/processed/hsp60_tier1_substrates.fasta"

# Known missing from K-12 proteome
ECOLI_MISSING_EXPECTED = {"P69408", "Q99390", "P62593", "P29368"}


def parse_accession_from_header(header):
    """Extract accession from UniProt FASTA header: >sp|ACCESSION|ENTRY ..."""
    parts = header.split("|")
    if len(parts) >= 2:
        return parts[1]
    return None


def build_accession_index(fasta_path):
    """Build dict: accession -> SeqRecord from a UniProt FASTA file."""
    idx = OrderedDict()
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = parse_accession_from_header(record.id)
        if acc:
            idx[acc] = record
    return idx


def fetch_from_uniprot(accession, max_retries=3):
    """Fetch a single FASTA from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200 and resp.text.strip().startswith(">"):
                # Parse the FASTA text
                from io import StringIO
                records = list(SeqIO.parse(StringIO(resp.text), "fasta"))
                if records:
                    return records[0]
            print(f"  WARNING: HTTP {resp.status_code} for {accession}", file=sys.stderr)
        except Exception as e:
            print(f"  WARNING: attempt {attempt+1} failed for {accession}: {e}", file=sys.stderr)
        time.sleep(1)
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Step C1: GroEL substrates
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("STEP C1: Extracting GroEL substrate FASTA sequences")
print("=" * 70)

groel_df = pd.read_csv(GROEL_TSV, sep="\t")
groel_accessions = set(groel_df["current_accession"].unique())
print(f"  GroEL substrates to extract: {len(groel_accessions)}")

# Build E. coli proteome index
print("  Indexing E. coli K-12 proteome...")
ecoli_idx = build_accession_index(ECOLI_FASTA)
print(f"  E. coli proteome entries: {len(ecoli_idx)}")

# Find matches and missing
groel_found = {}
groel_missing = []

for acc in groel_accessions:
    if acc in ecoli_idx:
        groel_found[acc] = ecoli_idx[acc]
    else:
        groel_missing.append(acc)

print(f"  Found in proteome: {len(groel_found)}")
print(f"  Missing from proteome: {len(groel_missing)}")
print(f"  Missing accessions: {groel_missing}")

# Fetch missing from UniProt
fetched_count = 0
for acc in groel_missing:
    print(f"  Fetching {acc} from UniProt...")
    record = fetch_from_uniprot(acc)
    if record:
        groel_found[acc] = record
        fetched_count += 1
        print(f"    OK: {record.description[:80]}")
    else:
        print(f"    FAILED: Could not fetch {acc}")

print(f"  Fetched from UniProt: {fetched_count}")
print(f"  Total sequences collected: {len(groel_found)}")

# Write FASTA - preserve order from TSV
written = 0
with open(GROEL_OUT, "w") as fh:
    for _, row in groel_df.iterrows():
        acc = row["current_accession"]
        if acc in groel_found:
            SeqIO.write(groel_found[acc], fh, "fasta")
            written += 1

print(f"  Written to {GROEL_OUT}: {written} sequences")

# Verify
from Bio import SeqIO as SeqIO2
verify_count = sum(1 for _ in SeqIO2.parse(GROEL_OUT, "fasta"))
print(f"  Verification: {verify_count} sequences in output FASTA")


# ═══════════════════════════════════════════════════════════════════════════════
# Step C2: HSP60 Tier 1 substrates
# ═══════════════════════════════════════════════════════════════════════════════
print()
print("=" * 70)
print("STEP C2: Extracting HSP60 Tier 1 substrate FASTA sequences")
print("=" * 70)

hsp60_df = pd.read_csv(HSP60_TSV, sep="\t")
hsp60_accessions = set(hsp60_df["uniprot_id"].unique())
print(f"  HSP60 Tier 1 substrates to extract: {len(hsp60_accessions)}")

# Build human proteome index
print("  Indexing human proteome...")
human_idx = build_accession_index(HUMAN_FASTA)
print(f"  Human proteome entries: {len(human_idx)}")

# Find matches and missing
hsp60_found = {}
hsp60_missing = []

for acc in hsp60_accessions:
    if acc in human_idx:
        hsp60_found[acc] = human_idx[acc]
    else:
        hsp60_missing.append(acc)

print(f"  Found in proteome: {len(hsp60_found)}")
print(f"  Missing from proteome: {len(hsp60_missing)}")
if hsp60_missing:
    print(f"  Missing accessions: {hsp60_missing}")

# Fetch missing from UniProt
fetched_count = 0
for acc in hsp60_missing:
    print(f"  Fetching {acc} from UniProt...")
    record = fetch_from_uniprot(acc)
    if record:
        hsp60_found[acc] = record
        fetched_count += 1
        print(f"    OK: {record.description[:80]}")
    else:
        print(f"    FAILED: Could not fetch {acc}")

if fetched_count > 0:
    print(f"  Fetched from UniProt: {fetched_count}")
print(f"  Total sequences collected: {len(hsp60_found)}")

# Write FASTA - preserve order from TSV
written = 0
with open(HSP60_OUT, "w") as fh:
    for _, row in hsp60_df.iterrows():
        acc = row["uniprot_id"]
        if acc in hsp60_found:
            SeqIO.write(hsp60_found[acc], fh, "fasta")
            written += 1

print(f"  Written to {HSP60_OUT}: {written} sequences")

# Verify
verify_count = sum(1 for _ in SeqIO2.parse(HSP60_OUT, "fasta"))
print(f"  Verification: {verify_count} sequences in output FASTA")

print()
print("Steps C1 & C2 complete.")
