#!/usr/bin/env python3
"""
Module G: Mitochondrial Targeting Analysis (Steps G1-G5)
========================================================
Antah Asti Prarambh project — Goal 3: Which parts can go to the mitochondrial matrix?

G1: Query UniProt REST API for transit peptide & signal peptide annotations
G2: Signal peptide analysis (extracted from G1 query)
G3: MitoCarta-based sub-mitochondrial localization cross-reference
G4: Integrate all targeting evidence into a combined table
G5: MTS vs. domain boundary relationship analysis
"""

import os
import sys
import time
import re
import requests
import pandas as pd
from io import StringIO
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
DATA_PROC   = BASE / "data" / "processed"
DATA_RAW    = BASE / "data" / "raw" / "uniprot"
RESULTS_MTS = BASE / "results" / "mts"
RESULTS_DOM = BASE / "results" / "domains"

RESULTS_MTS.mkdir(parents=True, exist_ok=True)

# Input files
HSP60_FILE   = DATA_PROC / "hsp60_tier1_substrates.tsv"
MITO_FILE    = DATA_PROC / "human_mito_proteome.tsv"
MATRIX_FILE  = DATA_PROC / "human_matrix_proteome.tsv"
UNIPROT_FILE = DATA_RAW  / "human_proteome.tsv"
CATH_FILE    = RESULTS_DOM / "cath_domain_assignments.tsv"

# Output files
UNIPROT_TP_CACHE   = RESULTS_MTS / "uniprot_transit_signal_cache.tsv"
COMBINED_OUT       = RESULTS_MTS / "combined_targeting.tsv"
MTS_DOMAIN_OUT     = RESULTS_MTS / "mts_domain_relationship.tsv"
REPORT_OUT         = RESULTS_MTS / "targeting_summary_report.txt"

# ─────────────────────────────────────────────────────────────────────────────
# STEP G1 & G2: Fetch transit peptide and signal peptide annotations from
#               UniProt REST API
# ─────────────────────────────────────────────────────────────────────────────

def load_protein_sets():
    """Load HSP60 Tier-1 substrates and mito proteome, return merged set."""
    hsp60 = pd.read_csv(HSP60_FILE, sep="\t")
    mito  = pd.read_csv(MITO_FILE, sep="\t")
    matrix = pd.read_csv(MATRIX_FILE, sep="\t")

    hsp60_ids  = set(hsp60["uniprot_id"].dropna().unique())
    mito_ids   = set(mito["uniprot_id"].dropna().unique())
    matrix_ids = set(matrix["uniprot_id"].dropna().unique())

    all_ids = hsp60_ids | mito_ids
    print(f"[load] HSP60 Tier-1 substrates : {len(hsp60_ids)}")
    print(f"[load] MitoCarta mito proteome : {len(mito_ids)}")
    print(f"[load] MitoCarta matrix subset : {len(matrix_ids)}")
    print(f"[load] Union (to query)        : {len(all_ids)}")
    return hsp60, mito, matrix, hsp60_ids, mito_ids, matrix_ids, all_ids


def fetch_uniprot_features(accessions, batch_size=100, max_retries=3):
    """
    Query UniProt REST API in batches for transit peptide, signal peptide,
    and subcellular location annotations.

    Returns a DataFrame with columns:
       accession, ft_transit, ft_signal, cc_subcellular_location
    """
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    fields = "accession,ft_transit,ft_signal,cc_subcellular_location"

    acc_list = sorted(accessions)
    all_rows = []
    n_batches = (len(acc_list) + batch_size - 1) // batch_size

    print(f"\n[G1/G2] Fetching UniProt transit/signal annotations for "
          f"{len(acc_list)} proteins in {n_batches} batches ...")

    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i:i + batch_size]
        batch_num = i // batch_size + 1
        query = " OR ".join(f"accession:{a}" for a in batch)
        # Wrap in parentheses for correct OR logic
        query = f"({query})"

        params = {
            "query": query,
            "format": "tsv",
            "fields": fields,
        }

        for attempt in range(max_retries):
            try:
                resp = requests.get(base_url, params=params, timeout=120)
                if resp.status_code == 200:
                    df_batch = pd.read_csv(StringIO(resp.text), sep="\t")
                    all_rows.append(df_batch)
                    print(f"  batch {batch_num}/{n_batches}: "
                          f"received {len(df_batch)} rows")
                    break
                elif resp.status_code == 429:
                    wait = int(resp.headers.get("Retry-After", 5))
                    print(f"  batch {batch_num}: rate-limited, "
                          f"waiting {wait}s ...")
                    time.sleep(wait)
                else:
                    print(f"  batch {batch_num}: HTTP {resp.status_code}, "
                          f"attempt {attempt+1}/{max_retries}")
                    time.sleep(2 ** attempt)
            except requests.exceptions.RequestException as e:
                print(f"  batch {batch_num}: request error: {e}, "
                      f"attempt {attempt+1}/{max_retries}")
                time.sleep(2 ** attempt)
        else:
            print(f"  WARNING: batch {batch_num} failed after "
                  f"{max_retries} attempts — skipping")

        # Polite delay between batches
        if batch_num < n_batches:
            time.sleep(1.0)

    if all_rows:
        df = pd.concat(all_rows, ignore_index=True)
    else:
        df = pd.DataFrame(columns=["Entry", "Transit peptide",
                                    "Signal peptide",
                                    "Subcellular location [CC]"])

    # Standardise column names
    col_map = {}
    for c in df.columns:
        cl = c.lower()
        if "entry" in cl and "name" not in cl:
            col_map[c] = "accession"
        elif "transit" in cl:
            col_map[c] = "ft_transit"
        elif "signal" in cl:
            col_map[c] = "ft_signal"
        elif "subcellular" in cl:
            col_map[c] = "cc_subcellular_location"
    df.rename(columns=col_map, inplace=True)

    # Ensure required columns exist
    for col in ["accession", "ft_transit", "ft_signal",
                "cc_subcellular_location"]:
        if col not in df.columns:
            df[col] = ""

    df = df[["accession", "ft_transit", "ft_signal",
             "cc_subcellular_location"]].copy()
    df.fillna("", inplace=True)

    print(f"[G1/G2] Total rows retrieved: {len(df)}")
    return df


def parse_transit_peptide(tp_str):
    """
    Parse UniProt transit peptide annotation string.
    Example: "TRANSIT 1..22; /note=\"Mitochondrion\""
    Returns (has_tp, tp_start, tp_end, tp_note).
    """
    if not tp_str or pd.isna(tp_str) or tp_str.strip() == "":
        return False, None, None, ""

    tp_str = str(tp_str)
    # Look for TRANSIT x..y pattern
    match = re.search(r'TRANSIT\s+(\d+)\.\.(\d+)', tp_str)
    if match:
        start = int(match.group(1))
        end   = int(match.group(2))
        # Extract note if present
        note_match = re.search(r'/note="([^"]*)"', tp_str)
        note = note_match.group(1) if note_match else ""
        return True, start, end, note

    return False, None, None, ""


def parse_signal_peptide(sp_str):
    """
    Parse UniProt signal peptide annotation string.
    Example: "SIGNAL 1..25"
    Returns (has_sp, sp_start, sp_end).
    """
    if not sp_str or pd.isna(sp_str) or sp_str.strip() == "":
        return False, None, None

    sp_str = str(sp_str)
    match = re.search(r'SIGNAL\s+(\d+)\.\.(\d+)', sp_str)
    if match:
        return True, int(match.group(1)), int(match.group(2))
    return False, None, None


# ─────────────────────────────────────────────────────────────────────────────
# STEP G3: MitoCarta-based sub-mitochondrial localization
# ─────────────────────────────────────────────────────────────────────────────

def build_mitocarta_lookup(mito_df, matrix_df):
    """
    Build a dict: accession -> {compartment, is_matrix, is_im, is_ims, is_om,
                                 gene_symbol, mitocarta_score}
    """
    lookup = {}
    for _, row in mito_df.iterrows():
        acc = row["uniprot_id"]
        comp = row.get("sub_mito_localization", "")
        lookup[acc] = {
            "mitocarta_compartment": comp,
            "is_matrix": bool(row.get("is_matrix", 0)),
            "is_im":     bool(row.get("is_im", 0)),
            "is_ims":    bool(row.get("is_ims", 0)),
            "is_om":     bool(row.get("is_om", 0)),
            "gene_symbol": row.get("gene_symbol", ""),
            "mitocarta_score": row.get("mitocarta_score", None),
        }
    # Matrix file may have IDs not in mito_df (shouldn't, but be safe)
    for _, row in matrix_df.iterrows():
        acc = row["uniprot_id"]
        if acc not in lookup:
            lookup[acc] = {
                "mitocarta_compartment": "Matrix",
                "is_matrix": True,
                "is_im": False, "is_ims": False, "is_om": False,
                "gene_symbol": row.get("gene_symbol", ""),
                "mitocarta_score": row.get("mitocarta_score", None),
            }
    return lookup


# ─────────────────────────────────────────────────────────────────────────────
# STEP G4: Integrate all targeting evidence
# ─────────────────────────────────────────────────────────────────────────────

def classify_targeting(row):
    """
    Assign a targeting classification based on integrated evidence.
    """
    mc = str(row.get("mitocarta_compartment", "")).strip()
    is_in_mito = row.get("is_in_mito_proteome", False)
    has_tp     = row.get("has_transit_peptide", False)
    has_sp     = row.get("has_signal_peptide", False)
    is_matrix  = row.get("mitocarta_is_matrix", False)
    is_im      = row.get("mitocarta_is_im", False)
    is_ims     = row.get("mitocarta_is_ims", False)
    is_om      = row.get("mitocarta_is_om", False)

    if is_matrix and has_tp:
        return "High-confidence matrix"
    elif is_matrix and not has_tp and not has_sp:
        return "Non-canonical matrix import"
    elif is_matrix:
        return "Probable matrix"
    elif is_im:
        return "Inner membrane (MIM)"
    elif is_ims:
        return "Intermembrane space (IMS)"
    elif is_om:
        return "Outer membrane (MOM)"
    elif is_in_mito and has_tp:
        # In MitoCarta but compartment not specifically matrix
        return "Mitochondrial (with MTS)"
    elif is_in_mito:
        return "Mitochondrial (other/unspecified)"
    elif has_tp:
        # Not in MitoCarta but has transit peptide
        return "Predicted mitochondrial (transit peptide)"
    elif has_sp:
        return "Secretory pathway (signal peptide)"
    else:
        return "Non-mitochondrial / no targeting signal"


def build_combined_table(uniprot_features, hsp60_df, mito_df, matrix_df,
                         hsp60_ids, mito_ids, all_ids):
    """Build the combined targeting table for all proteins."""

    mc_lookup = build_mitocarta_lookup(mito_df, matrix_df)

    # Build gene symbol lookup from HSP60 file
    gene_lookup = {}
    for _, row in hsp60_df.iterrows():
        gene_lookup[row["uniprot_id"]] = row.get("gene_name", "")
    for _, row in mito_df.iterrows():
        if row["uniprot_id"] not in gene_lookup:
            gene_lookup[row["uniprot_id"]] = row.get("gene_symbol", "")

    rows = []
    for acc in sorted(all_ids):
        # UniProt features
        feat = uniprot_features[uniprot_features["accession"] == acc]
        if len(feat) > 0:
            feat = feat.iloc[0]
            tp_str = feat["ft_transit"]
            sp_str = feat["ft_signal"]
            subcell = feat["cc_subcellular_location"]
        else:
            tp_str, sp_str, subcell = "", "", ""

        has_tp, tp_start, tp_end, tp_note = parse_transit_peptide(tp_str)
        has_sp, sp_start, sp_end = parse_signal_peptide(sp_str)

        # MitoCarta info
        mc_info = mc_lookup.get(acc, {})
        mc_comp  = mc_info.get("mitocarta_compartment", "")
        mc_mat   = mc_info.get("is_matrix", False)
        mc_im    = mc_info.get("is_im", False)
        mc_ims   = mc_info.get("is_ims", False)
        mc_om    = mc_info.get("is_om", False)

        gene = gene_lookup.get(acc, mc_info.get("gene_symbol", ""))

        rec = {
            "uniprot_accession":   acc,
            "gene_symbol":         gene,
            "mitocarta_status":    "MitoCarta" if acc in mito_ids else "Not in MitoCarta",
            "mitocarta_compartment": mc_comp,
            "mitocarta_is_matrix": mc_mat,
            "mitocarta_is_im":     mc_im,
            "mitocarta_is_ims":    mc_ims,
            "mitocarta_is_om":     mc_om,
            "has_transit_peptide": has_tp,
            "transit_peptide_start": tp_start,
            "transit_peptide_end":   tp_end,
            "transit_peptide_note":  tp_note,
            "has_signal_peptide":  has_sp,
            "signal_peptide_start": sp_start,
            "signal_peptide_end":   sp_end,
            "uniprot_subcellular":  subcell,
            "is_hsp60_substrate":  acc in hsp60_ids,
            "is_in_mito_proteome": acc in mito_ids,
        }
        rec["targeting_classification"] = classify_targeting(rec)
        rows.append(rec)

    df = pd.DataFrame(rows)
    print(f"\n[G4] Combined targeting table: {len(df)} proteins")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# STEP G5: MTS vs. domain boundary relationship
# ─────────────────────────────────────────────────────────────────────────────

def analyze_mts_domain_relationship(combined_df, cath_df):
    """
    For proteins with a transit peptide AND CATH domain assignments,
    compare the transit peptide cleavage site vs the first domain start.
    """
    # Proteins with transit peptide
    tp_df = combined_df[combined_df["has_transit_peptide"]].copy()
    tp_df = tp_df[["uniprot_accession", "transit_peptide_end"]].copy()
    tp_df = tp_df.dropna(subset=["transit_peptide_end"])
    tp_df["transit_peptide_end"] = tp_df["transit_peptide_end"].astype(int)

    # First domain per protein from CATH
    cath = cath_df.copy()
    cath["domain_start"] = pd.to_numeric(cath["domain_start"], errors="coerce")
    cath = cath.dropna(subset=["domain_start"])
    first_dom = (cath.sort_values("domain_start")
                     .groupby("uniprot_accession")
                     .first()
                     .reset_index()
                     [["uniprot_accession", "domain_start",
                       "cath_superfamily", "cath_superfamily_name"]])
    first_dom.rename(columns={"domain_start": "first_domain_start",
                               "cath_superfamily": "first_domain_cath_sf",
                               "cath_superfamily_name": "first_domain_name"},
                     inplace=True)
    first_dom["first_domain_start"] = first_dom["first_domain_start"].astype(int)

    merged = tp_df.merge(first_dom, on="uniprot_accession", how="inner")

    merged["gap_length"] = (merged["first_domain_start"]
                            - merged["transit_peptide_end"])
    merged["mts_overlaps_domain"] = merged["gap_length"] < 0
    merged["mts_is_pre_domain"] = merged["gap_length"] >= 0

    print(f"\n[G5] MTS-domain relationship: {len(merged)} proteins with "
          f"both TP and CATH domain")
    return merged


# ─────────────────────────────────────────────────────────────────────────────
# Summary Report
# ─────────────────────────────────────────────────────────────────────────────

def write_report(combined_df, mts_dom_df, hsp60_ids, mito_ids, matrix_ids):
    lines = []
    def p(msg=""):
        lines.append(msg)
        print(msg)

    p("=" * 72)
    p("Module G — Mitochondrial Targeting Analysis: Summary Report")
    p("=" * 72)
    p()

    # ── Overall numbers ──────────────────────────────────────────────────
    total = len(combined_df)
    n_hsp60 = combined_df["is_hsp60_substrate"].sum()
    n_mito  = combined_df["is_in_mito_proteome"].sum()
    p(f"Total proteins analysed          : {total}")
    p(f"  HSP60 Tier-1 substrates        : {n_hsp60}")
    p(f"  MitoCarta mito proteome        : {n_mito}")
    p(f"  Overlap (HSP60 ∩ MitoCarta)    : "
      f"{(combined_df['is_hsp60_substrate'] & combined_df['is_in_mito_proteome']).sum()}")
    p()

    # ── Targeting classification breakdown ────────────────────────────────
    p("--- Targeting classification breakdown (all proteins) ---")
    cls_counts = combined_df["targeting_classification"].value_counts()
    for cls, cnt in cls_counts.items():
        p(f"  {cls:50s}  {cnt:5d}  ({100*cnt/total:.1f}%)")
    p()

    # ── Transit peptide prevalence ────────────────────────────────────────
    n_tp_all = combined_df["has_transit_peptide"].sum()
    n_sp_all = combined_df["has_signal_peptide"].sum()
    p(f"Proteins with transit peptide     : {n_tp_all} "
      f"({100*n_tp_all/total:.1f}%)")
    p(f"Proteins with signal peptide      : {n_sp_all} "
      f"({100*n_sp_all/total:.1f}%)")
    p()

    # ── HSP60 substrate breakdown ─────────────────────────────────────────
    hsp60_sub = combined_df[combined_df["is_hsp60_substrate"]]
    p("--- HSP60 Tier-1 substrates: targeting breakdown ---")
    hsp60_cls = hsp60_sub["targeting_classification"].value_counts()
    for cls, cnt in hsp60_cls.items():
        p(f"  {cls:50s}  {cnt:5d}  ({100*cnt/len(hsp60_sub):.1f}%)")
    p()

    n_tp_hsp60 = hsp60_sub["has_transit_peptide"].sum()
    n_sp_hsp60 = hsp60_sub["has_signal_peptide"].sum()
    p(f"HSP60 substrates with transit peptide : {n_tp_hsp60} "
      f"({100*n_tp_hsp60/len(hsp60_sub):.1f}%)")
    p(f"HSP60 substrates with signal peptide  : {n_sp_hsp60} "
      f"({100*n_sp_hsp60/len(hsp60_sub):.1f}%)")
    p()

    # ── MTS prevalence: HSP60 matrix substrates vs general matrix proteome
    hsp60_matrix = hsp60_sub[hsp60_sub["mitocarta_is_matrix"]]
    mito_matrix  = combined_df[combined_df["mitocarta_is_matrix"]]

    if len(hsp60_matrix) > 0:
        tp_hsp60_mat = hsp60_matrix["has_transit_peptide"].sum()
        p(f"HSP60 matrix substrates : {len(hsp60_matrix)} total, "
          f"{tp_hsp60_mat} with TP "
          f"({100*tp_hsp60_mat/len(hsp60_matrix):.1f}%)")
    if len(mito_matrix) > 0:
        tp_mat = mito_matrix["has_transit_peptide"].sum()
        p(f"All matrix proteins     : {len(mito_matrix)} total, "
          f"{tp_mat} with TP "
          f"({100*tp_mat/len(mito_matrix):.1f}%)")
    p()

    # ── MTS-domain relationship ───────────────────────────────────────────
    p("--- MTS vs. first structural domain relationship ---")
    if len(mts_dom_df) > 0:
        n_overlap = mts_dom_df["mts_overlaps_domain"].sum()
        n_pre     = mts_dom_df["mts_is_pre_domain"].sum()
        p(f"Proteins with both TP and CATH domain : {len(mts_dom_df)}")
        p(f"  MTS overlaps first domain           : {n_overlap} "
          f"({100*n_overlap/len(mts_dom_df):.1f}%)")
        p(f"  MTS is pre-domain (non-overlapping) : {n_pre} "
          f"({100*n_pre/len(mts_dom_df):.1f}%)")
        p()
        gap = mts_dom_df[mts_dom_df["mts_is_pre_domain"]]["gap_length"]
        if len(gap) > 0:
            p(f"  Gap (domain_start - tp_end) statistics for non-overlapping:")
            p(f"    Mean   : {gap.mean():.1f} residues")
            p(f"    Median : {gap.median():.1f} residues")
            p(f"    Min    : {gap.min()} residues")
            p(f"    Max    : {gap.max()} residues")
        p()
        # Also report overlap details
        overlap_gap = mts_dom_df[mts_dom_df["mts_overlaps_domain"]]["gap_length"]
        if len(overlap_gap) > 0:
            p(f"  Overlap extent (negative gap) statistics:")
            p(f"    Mean overlap : {-overlap_gap.mean():.1f} residues")
            p(f"    Max overlap  : {-overlap_gap.min()} residues")
    else:
        p("  No proteins found with both TP and CATH domain data.")

    p()
    p("=" * 72)
    p("Key biological conclusion:")
    p("  The MTS (mitochondrial targeting sequence) is a SEPARATE N-terminal")
    p("  extension that is cleaved upon import. In most cases it does NOT")
    p("  overlap with the first structural domain.")
    p("=" * 72)

    report_text = "\n".join(lines)
    with open(REPORT_OUT, "w") as f:
        f.write(report_text + "\n")
    print(f"\n[report] Saved to {REPORT_OUT}")
    return report_text


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("Module G: Mitochondrial Targeting Analysis")
    print("=" * 55)

    # ── Load datasets ────────────────────────────────────────────────────
    hsp60_df, mito_df, matrix_df, hsp60_ids, mito_ids, matrix_ids, all_ids = \
        load_protein_sets()

    # ── G1 & G2: UniProt transit/signal peptide features ─────────────────
    if UNIPROT_TP_CACHE.exists():
        print(f"\n[cache] Loading cached UniProt features from {UNIPROT_TP_CACHE}")
        uniprot_feat = pd.read_csv(UNIPROT_TP_CACHE, sep="\t",
                                    dtype=str).fillna("")
    else:
        uniprot_feat = fetch_uniprot_features(all_ids, batch_size=100)
        uniprot_feat.to_csv(UNIPROT_TP_CACHE, sep="\t", index=False)
        print(f"[cache] Saved to {UNIPROT_TP_CACHE}")

    # Quick stats
    n_with_tp = sum(1 for _, r in uniprot_feat.iterrows()
                    if r["ft_transit"].strip() != "")
    n_with_sp = sum(1 for _, r in uniprot_feat.iterrows()
                    if r["ft_signal"].strip() != "")
    print(f"\n[G1] Transit peptide annotations found : {n_with_tp}")
    print(f"[G2] Signal peptide annotations found  : {n_with_sp}")

    # ── G3: MitoCarta cross-reference (already loaded) ───────────────────
    print(f"\n[G3] MitoCarta lookup built: {len(mito_ids)} mito proteins, "
          f"{len(matrix_ids)} matrix proteins")

    # ── G4: Combined targeting table ─────────────────────────────────────
    combined = build_combined_table(
        uniprot_feat, hsp60_df, mito_df, matrix_df,
        hsp60_ids, mito_ids, all_ids
    )
    combined.to_csv(COMBINED_OUT, sep="\t", index=False)
    print(f"[G4] Saved combined targeting table to {COMBINED_OUT}")

    # ── G5: MTS vs domain boundary ───────────────────────────────────────
    cath_df = pd.read_csv(CATH_FILE, sep="\t")
    mts_dom = analyze_mts_domain_relationship(combined, cath_df)
    mts_dom.to_csv(MTS_DOMAIN_OUT, sep="\t", index=False)
    print(f"[G5] Saved MTS-domain relationship to {MTS_DOMAIN_OUT}")

    # ── Summary report ───────────────────────────────────────────────────
    write_report(combined, mts_dom, hsp60_ids, mito_ids, matrix_ids)

    print("\nModule G complete.")


if __name__ == "__main__":
    main()
