#!/usr/bin/env python3
"""
Validate and remap all 252 UniProt accessions from the GroEL substrate dataset
against the current UniProt database. Write a cleaned, standardized output file.

Many old E. coli accessions have been DEMERGED into strain-specific entries.
This script follows the demerge chain and picks the E. coli K12 entry.
"""

import csv
import re
import time
import json
import sys
import requests
import pandas as pd
from pathlib import Path

# ─── Paths ───────────────────────────────────────────────────────────────────
BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
CLEAN_CSV = BASE / "data/raw/custom/kerner_2005_groel_interactors_clean.csv"
TABLE_S3  = BASE / "data/raw/custom/kerner_2005_groel_interactors_table_s3.csv"
OUTPUT    = BASE / "data/processed/groel_substrates_standardized.tsv"

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

# ─── 1. Read input files ────────────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Reading input files")
print("=" * 70)

df_clean = pd.read_csv(CLEAN_CSV)
df_s3 = pd.read_csv(TABLE_S3)

print(f"  Clean file: {len(df_clean)} rows")
print(f"  Table S3:   {len(df_s3)} rows")

# ─── 2. Parse SCOP folds from table_s3 raw_text_extracted ───────────────────
print("\n" + "=" * 70)
print("STEP 2: Parsing SCOP fold codes from table_s3")
print("=" * 70)

def parse_scop_folds(raw_text):
    """Extract SCOP fold codes like c.47, b.84, d.38, etc. from raw text."""
    if pd.isna(raw_text):
        return ""
    folds = re.findall(r'\b([a-g]\.\d+)\b', raw_text)
    seen = set()
    unique_folds = []
    for f in folds:
        if f not in seen:
            seen.add(f)
            unique_folds.append(f)
    return "; ".join(unique_folds)

def parse_s3_location(raw_text):
    """Extract location info from table S3 raw text as fallback.
    Patterns: 'Cytoplasmic. Cytoplasmic', 'Periplasmic. Unknown',
    'Integral Outer Membrane', 'Membrane-associated', etc.
    """
    if pd.isna(raw_text):
        return ""
    text = str(raw_text)
    # Look for location keywords in the text
    text_lower = text.lower()
    if "outer membrane" in text_lower:
        return "Outer membrane"
    if "inner membrane" in text_lower:
        return "Cell inner membrane"
    if "periplasmic" in text_lower:
        return "Periplasm"
    if "membrane-associated" in text_lower or "membrane associated" in text_lower:
        return "Membrane"
    if "integral" in text_lower and "membrane" in text_lower:
        return "Cell membrane"
    if "cytoplasmic" in text_lower or "cytoplasm" in text_lower:
        return "Cytoplasm"
    return ""

scop_map = {}
s3_loc_map = {}
for _, row in df_s3.iterrows():
    acc = row['accession']
    scop_map[acc] = parse_scop_folds(row.get('raw_text_extracted', ''))
    s3_loc_map[acc] = parse_s3_location(row.get('raw_text_extracted', ''))

sample = list(scop_map.items())[:5]
for acc, folds in sample:
    print(f"  {acc}: scop={folds}, s3_loc={s3_loc_map.get(acc, '')}")
print(f"  ... ({len(scop_map)} total)")

# ─── 3. Query UniProt for each accession ────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Querying UniProt REST API for all 252 accessions")
print("=" * 70)

def parse_entry(data):
    """Parse a full UniProt entry JSON into our fields."""
    entry = {
        "current_accession": data.get("primaryAccession", ""),
        "entry_name": data.get("uniProtkbId", ""),
        "gene_names": "",
        "protein_name": "",
        "organism": "",
        "length": "",
        "reviewed": "",
        "subcellular_location": "",
    }

    # Entry type / reviewed
    etype = data.get("entryType", "")
    if "Swiss-Prot" in str(etype) or etype == "UniProtKB reviewed (Swiss-Prot)":
        entry["reviewed"] = "reviewed"
    elif "TrEMBL" in str(etype) or "unreviewed" in str(etype).lower():
        entry["reviewed"] = "unreviewed"
    else:
        entry["reviewed"] = str(etype)

    # Gene names
    genes = data.get("genes", [])
    if genes:
        gn = []
        for g in genes:
            if "geneName" in g:
                gn.append(g["geneName"].get("value", ""))
        entry["gene_names"] = "; ".join(gn) if gn else ""

    # Protein name
    prot_desc = data.get("proteinDescription", {})
    rec_name = prot_desc.get("recommendedName", {})
    if rec_name:
        entry["protein_name"] = rec_name.get("fullName", {}).get("value", "")
    elif prot_desc.get("submissionNames"):
        entry["protein_name"] = prot_desc["submissionNames"][0].get("fullName", {}).get("value", "")

    # Organism
    entry["organism"] = data.get("organism", {}).get("scientificName", "")

    # Length
    entry["length"] = data.get("sequence", {}).get("length", "")

    # Subcellular location
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locs = []
            for subloc in comment.get("subcellularLocations", []):
                loc = subloc.get("location", {}).get("value", "")
                if loc:
                    locs.append(loc)
            entry["subcellular_location"] = "; ".join(locs)
            break

    return entry


def fetch_entry(acc, retries=2):
    """Fetch a single UniProt entry by accession. Returns (entry_data, status)."""
    for attempt in range(retries + 1):
        try:
            resp = requests.get(f"{UNIPROT_BASE}/{acc}", timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                etype = data.get("entryType", "")

                if etype == "Inactive":
                    reason = data.get("inactiveReason", {})
                    reason_type = reason.get("inactiveReasonType", "")
                    targets = reason.get("mergeDemergeTo", [])
                    return data, reason_type, targets
                else:
                    entry = parse_entry(data)
                    return entry, "ACTIVE", []

            elif resp.status_code == 400:
                return None, "NOT_FOUND", []
            elif resp.status_code == 429:
                time.sleep(3)
                continue
            else:
                if attempt < retries:
                    time.sleep(1)
                    continue
                return None, f"HTTP_{resp.status_code}", []
        except Exception as e:
            if attempt < retries:
                time.sleep(1)
                continue
            return None, f"ERROR:{e}", []
    return None, "FAILED", []


def resolve_accession(orig_acc):
    """
    Resolve an accession, following demerge/merge chains.
    Returns (entry_data_dict, status_string).
    """
    result, reason, targets = fetch_entry(orig_acc)

    if reason == "ACTIVE":
        return result, "current"

    if reason in ("DEMERGED", "MERGED"):
        # For demerged: try each target, prefer K12 (taxon 83333),
        # then prefer reviewed Swiss-Prot
        if not targets:
            return None, "obsolete"

        # Try to find the E. coli K12 entry among targets
        best_entry = None
        for target_acc in targets:
            time.sleep(0.3)
            t_result, t_reason, _ = fetch_entry(target_acc)
            if t_reason == "ACTIVE" and t_result:
                organism = t_result.get("organism", "")
                # Prefer K12 strain
                if "K12" in organism or "k12" in organism.lower():
                    if t_result.get("reviewed") == "reviewed":
                        return t_result, "demerged" if reason == "DEMERGED" else "merged"
                    if best_entry is None or best_entry.get("reviewed") != "reviewed":
                        best_entry = t_result
                elif best_entry is None:
                    best_entry = t_result

        if best_entry:
            return best_entry, "demerged" if reason == "DEMERGED" else "merged"

        return None, "obsolete"

    if reason == "DELETED":
        return None, "obsolete"

    if reason == "NOT_FOUND":
        return None, "obsolete"

    return None, "obsolete"


# Process all accessions
accessions = df_clean['accession'].tolist()
unique_accs = list(dict.fromkeys(accessions))  # dedupe preserving order
print(f"  Total unique accessions: {len(unique_accs)}")

acc_results = {}
status_summary = {}

for i, acc in enumerate(unique_accs):
    if (i + 1) % 10 == 0 or i == 0:
        print(f"  Processing {i+1}/{len(unique_accs)}: {acc}...")

    entry, status = resolve_accession(acc)
    acc_results[acc] = (entry, status)
    status_summary[status] = status_summary.get(status, 0) + 1

    time.sleep(0.25)  # Rate limiting

print(f"\n  Status summary:")
for s, c in sorted(status_summary.items()):
    print(f"    {s:12s}: {c}")

# ─── 4. Parse location category ─────────────────────────────────────────────
def categorize_location(subloc_text):
    """Categorize subcellular location into standard categories."""
    if not subloc_text or subloc_text == "":
        return "unknown"
    text_lower = subloc_text.lower()
    if "periplasm" in text_lower:
        return "periplasmic"
    if "outer membrane" in text_lower:
        return "outer_membrane"
    if "inner membrane" in text_lower or "cell inner membrane" in text_lower:
        return "inner_membrane"
    if "plasma membrane" in text_lower or "cell membrane" in text_lower:
        return "inner_membrane"
    if "membrane" in text_lower:
        return "inner_membrane"
    if "secret" in text_lower or "extracellular" in text_lower:
        return "secreted"
    if "cytoplasm" in text_lower or "cytosol" in text_lower:
        return "cytoplasmic"
    return "unknown"

# ─── 5. Build final merged dataframe ────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Building standardized output")
print("=" * 70)

rows = []
failed = []

for _, row in df_clean.iterrows():
    orig_acc = row['accession']
    entry, status = acc_results.get(orig_acc, (None, "obsolete"))

    if entry is None:
        failed.append(orig_acc)
        current_acc = orig_acc
        entry_name = row.get("entry_name", "")
        gene_sym = row.get("gene_symbol", "")
        protein_name = ""
        organism = "Escherichia coli"
        length = ""
        reviewed = ""
        subloc = ""
    else:
        current_acc = entry.get("current_accession", orig_acc)
        entry_name = entry.get("entry_name", row.get("entry_name", ""))
        gene_sym = entry.get("gene_names", "") or row.get("gene_symbol", "")
        protein_name = entry.get("protein_name", "")
        organism = entry.get("organism", "")
        length = entry.get("length", "")
        reviewed = entry.get("reviewed", "")
        subloc = entry.get("subcellular_location", "")

    # Use S3 location as fallback if UniProt has no annotation
    if not subloc or subloc == "" or str(subloc) == "nan":
        subloc_fallback = s3_loc_map.get(orig_acc, "")
        if subloc_fallback:
            subloc = subloc_fallback + " [from Kerner 2005]"

    loc_cat = categorize_location(subloc)
    scop = scop_map.get(orig_acc, "")

    out_row = {
        "original_accession": orig_acc,
        "current_accession": current_acc,
        "accession_status": status,
        "entry_name": entry_name,
        "gene_symbol": gene_sym,
        "protein_name": protein_name,
        "organism": organism,
        "length": length,
        "reviewed": reviewed,
        "groel_class": row.get("predicted_groel_class", ""),
        "mass_kDa": row.get("mass_kDa", ""),
        "subcellular_location": subloc,
        "location_category": loc_cat,
        "scop_folds": scop,
        "description_clean": row.get("description_clean", ""),
    }
    rows.append(out_row)

df_out = pd.DataFrame(rows)

# ─── 6. Write output ────────────────────────────────────────────────────────
print(f"\n  Writing {len(df_out)} rows to {OUTPUT}")
df_out.to_csv(OUTPUT, sep="\t", index=False)
print(f"  Done! File size: {OUTPUT.stat().st_size:,} bytes")

# ─── 7. Summary report ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY REPORT")
print("=" * 70)

print(f"\n  Total proteins: {len(df_out)}")

print(f"\n  Accession status:")
status_counts = df_out['accession_status'].value_counts()
for s in sorted(status_counts.index):
    print(f"    {s:12s}: {status_counts[s]}")

print(f"\n  GroEL class distribution:")
class_counts = df_out['groel_class'].value_counts()
for cls in sorted(class_counts.index):
    print(f"    Class {cls:8s}: {class_counts[cls]}")

print(f"\n  Location category distribution:")
loc_counts = df_out['location_category'].value_counts()
for loc in sorted(loc_counts.index):
    print(f"    {loc:20s}: {loc_counts[loc]}")

cyto = (df_out['location_category'] == 'cytoplasmic').sum()
non_cyto = len(df_out) - cyto
print(f"\n  Cytoplasmic vs non-cytoplasmic:")
print(f"    Cytoplasmic:     {cyto}")
print(f"    Non-cytoplasmic: {non_cyto}")

print(f"\n  Reviewed status:")
rev_counts = df_out['reviewed'].value_counts()
for rev in sorted(rev_counts.index):
    print(f"    {rev:12s}: {rev_counts[rev]}")

if failed:
    print(f"\n  Failed to resolve ({len(failed)}):")
    for acc in failed:
        print(f"    {acc}")
else:
    print(f"\n  All accessions resolved successfully!")

print(f"\n  Sample rows:")
print(df_out[['original_accession', 'current_accession', 'accession_status', 'gene_symbol', 'groel_class', 'location_category', 'scop_folds']].head(10).to_string(index=False))

print(f"\n  Output written to: {OUTPUT}")
print("=" * 70)
