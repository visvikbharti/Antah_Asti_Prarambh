#!/usr/bin/env python3
"""
parse_mitocarta.py
==================
Parse MitoCarta 3.0 Excel file and produce processed datasets for the
Antah Asti Prarambh project.

MitoCarta 3.0 Reference:
    Rath S et al. "MitoCarta3.0: an updated mitochondrial proteome now with
    sub-organelle localization and pathway annotations."
    Nucleic Acids Research, 2021; 49(D1):D1541-D1547.
    PMID: 33174596  DOI: 10.1093/nar/gkaa1011

Download URL (if file not already present):
    https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls
    Landing page: https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta3.0.html

Usage:
    python parse_mitocarta.py [--download]

    --download   Attempt to download the MitoCarta 3.0 file before parsing.
                 Without this flag, the script expects the file to already be
                 at data/raw/mitocarta/Human.MitoCarta3.0.xls

Outputs:
    data/processed/human_mito_proteome.tsv       All 1,136 mitochondrial proteins
    data/processed/human_matrix_proteome.tsv      Matrix-localized subset
    data/processed/mitocarta_summary_report.txt   Summary with HSP60 cross-reference

Dependencies:
    pip install pandas openpyxl requests xlrd
"""

import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Paths (relative to project root)
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

MITOCARTA_RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "mitocarta")
MITOCARTA_FILE = os.path.join(MITOCARTA_RAW_DIR, "Human.MitoCarta3.0.xls")

PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FULL_PROTEOME_OUT = os.path.join(PROCESSED_DIR, "human_mito_proteome.tsv")
MATRIX_PROTEOME_OUT = os.path.join(PROCESSED_DIR, "human_matrix_proteome.tsv")
REPORT_OUT = os.path.join(PROCESSED_DIR, "mitocarta_summary_report.txt")

HSP60_FILE = os.path.join(PROJECT_ROOT, "data", "raw", "custom", "hsp60_interactome_clean.tsv")

DOWNLOAD_URL = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_mitocarta():
    """Download MitoCarta 3.0 from the Broad Institute."""
    import requests

    os.makedirs(MITOCARTA_RAW_DIR, exist_ok=True)
    if os.path.exists(MITOCARTA_FILE):
        print(f"[INFO] File already exists: {MITOCARTA_FILE}")
        return

    print(f"[INFO] Downloading MitoCarta 3.0 from {DOWNLOAD_URL} ...")
    resp = requests.get(DOWNLOAD_URL, timeout=120)
    resp.raise_for_status()
    with open(MITOCARTA_FILE, "wb") as fh:
        fh.write(resp.content)
    print(f"[INFO] Saved {len(resp.content)} bytes to {MITOCARTA_FILE}")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_mitocarta():
    """Parse MitoCarta 3.0 and return a clean DataFrame."""
    if not os.path.exists(MITOCARTA_FILE):
        sys.exit(
            f"[ERROR] MitoCarta file not found at {MITOCARTA_FILE}\n"
            f"        Download it from: {DOWNLOAD_URL}\n"
            f"        Or re-run with --download flag."
        )

    print(f"[INFO] Reading {MITOCARTA_FILE} ...")
    df = pd.read_excel(MITOCARTA_FILE, sheet_name="A Human MitoCarta3.0")
    print(f"[INFO] Loaded {len(df)} proteins from MitoCarta 3.0")

    result = pd.DataFrame()
    result["uniprot_id"] = df["UniProt"]
    result["gene_symbol"] = df["Symbol"]
    result["protein_name"] = df["Description"]
    result["mitocarta_score"] = df["MitoCarta2.0_Score"]
    result["sub_mito_localization"] = df["MitoCarta3.0_SubMitoLocalization"]
    result["pathways"] = df["MitoCarta3.0_MitoPathways"].fillna("")

    # Binary localization flags
    result["is_matrix"] = (
        result["sub_mito_localization"]
        .str.contains("Matrix", case=False, na=False)
        .astype(int)
    )
    result["is_im"] = (
        result["sub_mito_localization"]
        .str.contains("MIM|Membrane", case=False, na=False)
        .astype(int)
    )
    result["is_ims"] = (
        result["sub_mito_localization"]
        .str.contains("IMS", case=False, na=False)
        .astype(int)
    )
    result["is_om"] = (
        result["sub_mito_localization"]
        .str.contains("MOM", case=False, na=False)
        .astype(int)
    )

    return df, result


# ---------------------------------------------------------------------------
# Cross-reference with HSP60 interactome
# ---------------------------------------------------------------------------

def cross_reference_hsp60(result):
    """Compare HSP60 interactome annotations against MitoCarta 3.0."""
    if not os.path.exists(HSP60_FILE):
        print(f"[WARN] HSP60 file not found at {HSP60_FILE}; skipping cross-reference.")
        return None

    hsp = pd.read_csv(HSP60_FILE, sep="\t")
    print(f"[INFO] Loaded {len(hsp)} HSP60 interactors")

    # Build lookups (first occurrence wins to avoid Series ambiguity)
    mc3_by_uniprot = {}
    mc3_by_gene = {}
    for _, row in result.iterrows():
        uid = str(row["uniprot_id"]).strip()
        gene = str(row["gene_symbol"]).strip().upper()
        mc3_by_uniprot.setdefault(uid, row)
        mc3_by_gene.setdefault(gene, row)

    hsp_in_mc2 = hsp[hsp["mitocarta"] == "X"]
    hsp_matrix_mc2 = (hsp["matrix"] == "X").sum()

    gained, lost, matrix_changes = [], [], []
    hsp_matrix_mc3_count = 0
    total_in_mc3 = 0
    not_in_mc3_list = []

    for _, hrow in hsp.iterrows():
        uid = str(hrow["uniprot_id"]).strip()
        gene = str(hrow["gene_name"]).strip().upper()
        in_mc2 = hrow["mitocarta"] == "X"
        was_matrix_mc2 = hrow["matrix"] == "X"

        mc3_row = mc3_by_uniprot.get(uid) or mc3_by_gene.get(gene)
        in_mc3 = mc3_row is not None

        if in_mc3:
            total_in_mc3 += 1
            is_matrix_mc3 = int(mc3_row["is_matrix"]) == 1
            sub_loc = str(mc3_row["sub_mito_localization"])
            if is_matrix_mc3:
                hsp_matrix_mc3_count += 1
            if not in_mc2:
                gained.append((uid, gene, sub_loc))
            if in_mc2:
                if was_matrix_mc2 and not is_matrix_mc3:
                    matrix_changes.append(
                        (uid, gene, "Matrix in MC2 -> NOT Matrix in MC3", sub_loc)
                    )
                elif not was_matrix_mc2 and is_matrix_mc3:
                    matrix_changes.append(
                        (uid, gene, "NOT Matrix in MC2 -> Matrix in MC3", sub_loc)
                    )
        else:
            not_in_mc3_list.append((uid, gene, "Yes" if in_mc2 else "No"))
            if in_mc2:
                lost.append((uid, gene))

    return {
        "hsp": hsp,
        "hsp_in_mc2": len(hsp_in_mc2),
        "hsp_matrix_mc2": hsp_matrix_mc2,
        "total_in_mc3": total_in_mc3,
        "hsp_matrix_mc3": hsp_matrix_mc3_count,
        "not_in_mc3": len(hsp) - total_in_mc3,
        "gained": gained,
        "lost": lost,
        "matrix_changes": matrix_changes,
        "not_in_mc3_list": not_in_mc3_list,
    }


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_report(df_raw, result, matrix_df, xref):
    """Write a human-readable summary report."""
    os.makedirs(os.path.dirname(REPORT_OUT), exist_ok=True)
    with open(REPORT_OUT, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("MitoCarta 3.0 Analysis Summary Report\n")
        f.write("=" * 70 + "\n\n")

        f.write("SOURCE\n" + "-" * 40 + "\n")
        f.write("Database: MitoCarta 3.0 (Rath et al., Nucleic Acids Research 2021)\n")
        f.write("PMID: 33174596\n")
        f.write("DOI: 10.1093/nar/gkaa1011\n")
        f.write(f"Downloaded from: {DOWNLOAD_URL}\n\n")

        f.write("MITOCARTA 3.0 OVERVIEW\n" + "-" * 40 + "\n")
        f.write(f"Total human mitochondrial proteins: {len(result)}\n\n")

        f.write("Sub-mitochondrial compartment breakdown (raw annotations):\n")
        for loc, count in df_raw["MitoCarta3.0_SubMitoLocalization"].value_counts().items():
            f.write(f"  {loc:20s}: {count:5d} ({count / len(result) * 100:5.1f}%)\n")
        f.write("\n")

        f.write("Binary localization flags:\n")
        f.write(f"  Matrix:                {result['is_matrix'].sum():5d}\n")
        f.write(f"  Inner membrane (MIM):  {result['is_im'].sum():5d}\n")
        f.write(f"  Intermembrane space:   {result['is_ims'].sum():5d}\n")
        f.write(f"  Outer membrane:        {result['is_om'].sum():5d}\n")
        unknown = (df_raw["MitoCarta3.0_SubMitoLocalization"] == "unknown").sum()
        f.write(f"  Unknown:               {unknown:5d}\n\n")

        f.write("OUTPUT FILES\n" + "-" * 40 + "\n")
        f.write(f"1. human_mito_proteome.tsv   ({len(result)} proteins)\n")
        f.write(f"2. human_matrix_proteome.tsv ({len(matrix_df)} proteins)\n\n")

        if xref is not None:
            hsp = xref["hsp"]
            f.write("=" * 70 + "\n")
            f.write("CROSS-REFERENCE: HSP60 Interactome vs MitoCarta 3.0\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Total HSP60 interactors: {len(hsp)}\n\n")
            f.write("Membership comparison:\n")
            f.write(f"  In MC2: {xref['hsp_in_mc2']:4d} / {len(hsp)}\n")
            f.write(f"  In MC3: {xref['total_in_mc3']:4d} / {len(hsp)}\n")
            f.write(f"  Not in MC3: {xref['not_in_mc3']:4d} / {len(hsp)}\n\n")

            f.write("Matrix annotation comparison:\n")
            f.write(f"  Matrix in MC2: {xref['hsp_matrix_mc2']}\n")
            f.write(f"  Matrix in MC3: {xref['hsp_matrix_mc3']}\n\n")

            f.write("=" * 70 + "\n")
            f.write("DISCREPANCIES\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Gained in MC3: {len(xref['gained'])}\n")
            for uid, gene, loc in xref["gained"]:
                f.write(f"  {uid:10s} {gene:12s} MC3: {loc}\n")

            f.write(f"\nLost from MC3: {len(xref['lost'])}\n")
            for uid, gene in xref["lost"]:
                f.write(f"  {uid:10s} {gene}\n")

            f.write(f"\nMatrix localization changes: {len(xref['matrix_changes'])}\n")
            for uid, gene, change, loc in xref["matrix_changes"]:
                f.write(f"  {uid:10s} {gene:12s} {change}\n")

            f.write(f"\nNot in MC3 at all: {len(xref['not_in_mc3_list'])}\n")
            f.write(f"  {'UniProt':12s} {'Gene':12s} {'Was MC2?'}\n")
            f.write(f"  {'-'*12} {'-'*12} {'-'*8}\n")
            for uid, gene, was in xref["not_in_mc3_list"]:
                f.write(f"  {uid:12s} {gene:12s} {was}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 70 + "\n")

    print(f"[INFO] Report saved to {REPORT_OUT}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Parse MitoCarta 3.0")
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download MitoCarta 3.0 file before parsing",
    )
    args = parser.parse_args()

    if args.download:
        download_mitocarta()

    # Parse
    df_raw, result = parse_mitocarta()

    # Save full proteome
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    result.to_csv(FULL_PROTEOME_OUT, sep="\t", index=False)
    print(f"[INFO] Saved {FULL_PROTEOME_OUT} ({len(result)} proteins)")

    # Save matrix proteome
    matrix_df = result[result["is_matrix"] == 1].copy()
    matrix_df.to_csv(MATRIX_PROTEOME_OUT, sep="\t", index=False)
    print(f"[INFO] Saved {MATRIX_PROTEOME_OUT} ({len(matrix_df)} proteins)")

    # Cross-reference with HSP60
    xref = cross_reference_hsp60(result)

    # Write report
    write_report(df_raw, result, matrix_df, xref)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()
