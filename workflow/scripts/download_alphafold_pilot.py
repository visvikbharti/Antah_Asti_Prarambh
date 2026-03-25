#!/usr/bin/env python3
"""
Download AlphaFold structures for all pilot proteins in Antah Asti Prarambh.
Then build a structure index with pLDDT statistics.

Module D (Steps D1-D2).

Note: AlphaFold DB has migrated to v6 as of 2025. We download v6 CIF files
and fall back to v4 if v6 is unavailable for a given accession.
"""

import os
import sys
import time
import re
import statistics
from collections import defaultdict

import pandas as pd
import requests

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT = "/Users/vishalbharti/Downloads/Antah_Asti_Prarambh"
GROEL_TSV   = os.path.join(PROJECT, "data/processed/groel_substrates_standardized.tsv")
HSP60_TSV   = os.path.join(PROJECT, "data/processed/hsp60_tier1_substrates.tsv")
MITO_TSV    = os.path.join(PROJECT, "data/processed/human_mito_proteome.tsv")
OUT_DIR     = os.path.join(PROJECT, "data/raw/alphafold/pilot")
INDEX_PATH  = os.path.join(PROJECT, "results/structures/structure_index.tsv")

AF_URL_V6 = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif"
AF_URL_V4 = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.cif"

BATCH_SIZE = 50
DELAY      = 0.5   # seconds between batches
TIMEOUT    = 30     # per-request timeout

# Regex for valid UniProt accession
UNIPROT_RE = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')


def is_valid_accession(acc):
    """Check if a string looks like a valid UniProt accession."""
    if not acc or len(acc) < 6:
        return False
    # Allow standard 6 or 10-char UniProt accessions
    return bool(re.match(r'^[A-Z][0-9A-Z]{4,9}$', acc))


# ── Step 1: Collect unique accessions with source dataset labels ───────────
def collect_accessions():
    """Return dict  accession -> set of source labels."""
    acc_sources = defaultdict(set)

    # GroEL substrates
    df = pd.read_csv(GROEL_TSV, sep="\t")
    for acc in df["current_accession"].dropna().unique():
        acc = str(acc).strip()
        if is_valid_accession(acc):
            acc_sources[acc].add("groel")
    print(f"  GroEL unique accessions:  {sum(1 for s in acc_sources.values() if 'groel' in s)}")

    # HSP60 Tier 1 substrates
    df = pd.read_csv(HSP60_TSV, sep="\t")
    for acc in df["uniprot_id"].dropna().unique():
        acc = str(acc).strip()
        if is_valid_accession(acc):
            acc_sources[acc].add("hsp60")
    print(f"  HSP60 unique accessions:  {sum(1 for s in acc_sources.values() if 'hsp60' in s)}")

    # Human mito proteome (includes matrix subset)
    df = pd.read_csv(MITO_TSV, sep="\t")
    for _, row in df.iterrows():
        acc = str(row["uniprot_id"]).strip()
        if not is_valid_accession(acc):
            continue
        acc_sources[acc].add("mito")
        # Check if it is matrix-localised
        if row.get("is_matrix", 0) == 1:
            acc_sources[acc].add("matrix")
    mito_count = sum(1 for s in acc_sources.values() if "mito" in s)
    matrix_count = sum(1 for s in acc_sources.values() if "matrix" in s)
    print(f"  Mito unique accessions:   {mito_count}")
    print(f"  Matrix subset:            {matrix_count}")
    print(f"  Total unique accessions:  {len(acc_sources)}")
    return acc_sources


# ── Step 2: Download CIF files ─────────────────────────────────────────────
def download_structures(acc_sources):
    """Download AlphaFold CIF files. Try v6 first, fall back to v4.
    Returns (succeeded, failed, skipped) lists.
    succeeded: list of (accession, version_str, filename)
    failed: list of (accession, reason)
    """
    os.makedirs(OUT_DIR, exist_ok=True)

    accessions = sorted(acc_sources.keys())
    total = len(accessions)
    succeeded = []
    failed = []
    skipped = []

    session = requests.Session()
    session.headers.update({
        "User-Agent": "AntahAstiPrarambh/1.0 (research project; AlphaFold bulk download)"
    })

    print(f"\nDownloading {total} AlphaFold structures …")
    t0 = time.time()

    for i, acc in enumerate(accessions, 1):
        # Check if any version already exists
        fname_v6 = f"AF-{acc}-F1-model_v6.cif"
        fname_v4 = f"AF-{acc}-F1-model_v4.cif"
        fpath_v6 = os.path.join(OUT_DIR, fname_v6)
        fpath_v4 = os.path.join(OUT_DIR, fname_v4)

        # Skip if already downloaded (either version)
        for fpath, ver in [(fpath_v6, "v6"), (fpath_v4, "v4")]:
            if os.path.exists(fpath) and os.path.getsize(fpath) > 100:
                skipped.append(acc)
                succeeded.append((acc, ver, os.path.basename(fpath)))
                if i % 100 == 0 or i == total:
                    elapsed = time.time() - t0
                    print(f"  [{i}/{total}]  {acc}  SKIPPED (exists)  "
                          f"elapsed={elapsed:.0f}s  ok={len(succeeded)} fail={len(failed)}")
                break
        else:
            # Try v6 first
            downloaded = False
            for url_tmpl, fname, ver in [
                (AF_URL_V6, fname_v6, "v6"),
                (AF_URL_V4, fname_v4, "v4"),
            ]:
                url = url_tmpl.format(acc=acc)
                try:
                    resp = session.get(url, timeout=TIMEOUT, stream=True)
                    if resp.status_code == 200:
                        fpath = os.path.join(OUT_DIR, fname)
                        with open(fpath, "wb") as f:
                            for chunk in resp.iter_content(chunk_size=65536):
                                f.write(chunk)
                        succeeded.append((acc, ver, fname))
                        downloaded = True
                        break
                    elif resp.status_code == 404:
                        continue  # try next version
                    else:
                        # Non-404 error, still try next version
                        continue
                except requests.exceptions.RequestException:
                    continue

            if not downloaded:
                failed.append((acc, "404 (v6 and v4)"))

            # Polite delay between requests
            if i % BATCH_SIZE == 0:
                time.sleep(DELAY)

        # Progress
        if i % 100 == 0 or i == total:
            elapsed = time.time() - t0
            print(f"  [{i}/{total}]  elapsed={elapsed:.0f}s  "
                  f"ok={len(succeeded)} fail={len(failed)} skip={len(skipped)}")

    print(f"\nDownload complete: {len(succeeded)} succeeded, "
          f"{len(failed)} failed, {len(skipped)} skipped (already existed)")
    return succeeded, failed, skipped


# ── Step 3: Parse a CIF file for pLDDT and residue count (CA atoms) ───────
def parse_cif_fast(filepath):
    """
    Fast text-based parser for AlphaFold mmCIF files.
    Extracts CA atom records to get residue count and pLDDT values.

    Returns (residues_modeled, mean_plddt, min_plddt, frac_gt70) or Nones.
    """
    try:
        in_atom_site = False
        col_names = []
        b_col = None      # index of B_iso_or_equiv
        atom_col = None    # index of label_atom_id (to filter CA)
        seq_col = None     # index of label_seq_id (residue number)
        columns_resolved = False

        plddt_values = []
        residue_ids = set()

        with open(filepath, "r") as f:
            for line in f:
                stripped = line.strip()

                # Detect start of _atom_site loop
                if stripped.startswith("_atom_site."):
                    in_atom_site = True
                    col_name = stripped.split(".")[1].strip()
                    col_names.append(col_name)
                    continue

                if in_atom_site and not stripped.startswith("_atom_site."):
                    # We've finished reading column headers
                    if not columns_resolved:
                        columns_resolved = True
                        for idx, cn in enumerate(col_names):
                            if cn == "B_iso_or_equiv":
                                b_col = idx
                            elif cn == "label_atom_id":
                                atom_col = idx
                            elif cn == "label_seq_id":
                                seq_col = idx
                        if b_col is None:
                            return None, None, None, None

                    # Now parse data lines
                    if stripped.startswith("ATOM") or stripped.startswith("HETATM"):
                        parts = stripped.split()
                        max_needed = max(b_col, atom_col or 0, seq_col or 0)
                        if len(parts) > max_needed:
                            atom_name = parts[atom_col] if atom_col is not None else ""
                            if atom_name == "CA":
                                try:
                                    plddt = float(parts[b_col])
                                    plddt_values.append(plddt)
                                except ValueError:
                                    pass
                                if seq_col is not None:
                                    residue_ids.add(parts[seq_col])
                    elif stripped.startswith("#") or stripped.startswith("loop_") or stripped.startswith("_"):
                        # End of atom_site block
                        in_atom_site = False

        if not plddt_values:
            return None, None, None, None

        n_res = len(residue_ids) if residue_ids else len(plddt_values)
        mean_p = statistics.mean(plddt_values)
        min_p  = min(plddt_values)
        frac70 = sum(1 for v in plddt_values if v > 70) / len(plddt_values)
        return n_res, round(mean_p, 2), round(min_p, 2), round(frac70, 4)

    except Exception as e:
        print(f"    WARNING: parse error for {filepath}: {e}")
        return None, None, None, None


# ── Step 4: Build the structure index ──────────────────────────────────────
def build_index(acc_sources, succeeded, failed):
    """Build structure_index.tsv with per-protein metadata."""
    os.makedirs(os.path.dirname(INDEX_PATH), exist_ok=True)

    # Build lookup: accession -> (version, filename)
    acc_info = {}
    for acc, ver, fname in succeeded:
        acc_info[acc] = (ver, fname)

    all_accs = sorted(acc_sources.keys())

    rows = []
    n_parsed = 0
    t0 = time.time()

    print(f"\nBuilding structure index for {len(all_accs)} proteins …")

    for i, acc in enumerate(all_accs, 1):
        sources = ",".join(sorted(acc_sources[acc]))

        if acc in acc_info:
            ver, fname = acc_info[acc]
            fpath = os.path.join(OUT_DIR, fname)
            rel_path = os.path.relpath(fpath, PROJECT)
            fsize = os.path.getsize(fpath) if os.path.exists(fpath) else 0

            n_res, mean_p, min_p, frac70 = parse_cif_fast(fpath)
            n_parsed += 1
            rows.append({
                "uniprot_accession": acc,
                "model_path": rel_path,
                "model_version": ver,
                "file_size_bytes": fsize,
                "has_structure": True,
                "source_dataset": sources,
                "residues_modeled": n_res if n_res else "",
                "mean_plddt": mean_p if mean_p is not None else "",
                "min_plddt": min_p if min_p is not None else "",
                "fraction_plddt_gt70": frac70 if frac70 is not None else "",
            })
        else:
            rows.append({
                "uniprot_accession": acc,
                "model_path": "",
                "model_version": "",
                "file_size_bytes": 0,
                "has_structure": False,
                "source_dataset": sources,
                "residues_modeled": "",
                "mean_plddt": "",
                "min_plddt": "",
                "fraction_plddt_gt70": "",
            })

        if i % 200 == 0 or i == len(all_accs):
            elapsed = time.time() - t0
            print(f"  [{i}/{len(all_accs)}]  parsed={n_parsed}  elapsed={elapsed:.0f}s")

    df = pd.DataFrame(rows)
    df.to_csv(INDEX_PATH, sep="\t", index=False)
    print(f"\nStructure index written to: {INDEX_PATH}")
    return df


# ── Step 5: Summary statistics ─────────────────────────────────────────────
def print_summary(df, failed):
    """Print a summary report."""
    n_struct = df["has_structure"].sum()
    n_no     = (~df["has_structure"]).sum()
    total_bytes = df["file_size_bytes"].sum()
    total_mb = total_bytes / (1024 * 1024)

    print("\n" + "=" * 65)
    print("  MODULE D — STRUCTURE DOWNLOAD & INDEX SUMMARY")
    print("=" * 65)
    print(f"  Total proteins in index:          {len(df)}")
    print(f"  Structures downloaded:            {n_struct}")
    print(f"  Structures missing (no AF model): {n_no}")
    print(f"  Disk usage:                       {total_mb:.1f} MB")

    # Version breakdown
    if n_struct > 0:
        print(f"\n  Version breakdown:")
        for ver in sorted(df.loc[df["has_structure"], "model_version"].unique()):
            n = (df["model_version"] == ver).sum()
            print(f"    {ver}: {n}")

    # pLDDT stats for proteins with structures
    plddt = df.loc[df["has_structure"], "mean_plddt"]
    plddt = pd.to_numeric(plddt, errors="coerce").dropna()
    if len(plddt) > 0:
        print(f"\n  pLDDT distribution (mean per protein):")
        print(f"    count:   {len(plddt)}")
        print(f"    mean:    {plddt.mean():.1f}")
        print(f"    median:  {plddt.median():.1f}")
        print(f"    std:     {plddt.std():.1f}")
        print(f"    min:     {plddt.min():.1f}")
        print(f"    max:     {plddt.max():.1f}")
        print(f"    Q25:     {plddt.quantile(0.25):.1f}")
        print(f"    Q75:     {plddt.quantile(0.75):.1f}")

    frac = df.loc[df["has_structure"], "fraction_plddt_gt70"]
    frac = pd.to_numeric(frac, errors="coerce").dropna()
    if len(frac) > 0:
        print(f"\n  Fraction of residues with pLDDT > 70:")
        print(f"    mean:    {frac.mean():.3f}")
        print(f"    median:  {frac.median():.3f}")

    if failed:
        print(f"\n  Failed downloads ({len(failed)}):")
        for acc, reason in failed[:20]:
            print(f"    {acc}: {reason}")
        if len(failed) > 20:
            print(f"    … and {len(failed) - 20} more")

    # Source breakdown
    print(f"\n  Source dataset breakdown:")
    for src in ["groel", "hsp60", "mito", "matrix"]:
        mask = df["source_dataset"].str.contains(src, na=False)
        n = mask.sum()
        n_ok = (mask & df["has_structure"]).sum()
        print(f"    {src:10s}: {n:5d} total,  {n_ok:5d} with structure")

    print("=" * 65)


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("Module D: AlphaFold Structure Download & Index\n")

    # D1a: Collect accessions
    print("Step 1: Collecting unique accessions …")
    acc_sources = collect_accessions()

    # D1b: Download
    print("\nStep 2: Downloading AlphaFold structures …")
    succeeded, failed, skipped = download_structures(acc_sources)

    # D2: Build index
    print("\nStep 3: Building structure index …")
    df = build_index(acc_sources, succeeded, failed)

    # Summary
    print_summary(df, failed)

    return 0


if __name__ == "__main__":
    sys.exit(main())
