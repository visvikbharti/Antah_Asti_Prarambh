#!/usr/bin/env python3
"""
Module D, Step D3: Run DSSP on all pilot AlphaFold structures.

Assigns secondary structure to every residue using mkdssp, then produces:
  1. Per-protein summary  -> results/structures/dssp_summary.tsv
  2. Per-residue table    -> results/structures/dssp_per_residue.tsv
     (skipped if estimated size > 500 MB)
  3. Individual DSSP files -> results/structures/dssp/{ACCESSION}.dssp
"""

import os, sys, glob, subprocess, re, time
from pathlib import Path
from collections import Counter

# ── paths ─────────────────────────────────────────────────────────────
BASE      = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
PILOT_DIR = BASE / "data/raw/alphafold/pilot"
OUT_DIR   = BASE / "results/structures"
DSSP_DIR  = OUT_DIR / "dssp"
MKDSSP    = "/Users/vishalbharti/opt/anaconda3/bin/mkdssp"
TIMEOUT   = 30  # seconds per structure

DSSP_DIR.mkdir(parents=True, exist_ok=True)

# ── DSSP code grouping ───────────────────────────────────────────────
HELIX_CODES  = set("HGI")
STRAND_CODES = set("EB")
COIL_CODES   = set("TS -")   # T, S, '-', and ' ' (space = coil in DSSP)

def classify(code):
    """Return 'helix', 'strand', or 'coil' for a one-character DSSP code."""
    if code in HELIX_CODES:
        return "helix"
    if code in STRAND_CODES:
        return "strand"
    return "coil"

# ── parse DSSP output file ───────────────────────────────────────────
def parse_dssp(dssp_path):
    """
    Parse a DSSP file and return a list of (residue_number, residue_name, ss_code)
    tuples. ss_code is one of H,G,I,E,B,T,S,- (space mapped to '-').
    """
    residues = []
    in_residue_section = False
    with open(dssp_path) as fh:
        for line in fh:
            # The residue section starts after the header line beginning with
            # "  #  RESIDUE AA STRUCTURE"
            if line.startswith("  #  RESIDUE AA"):
                in_residue_section = True
                continue
            if not in_residue_section:
                continue
            # Skip chain-break lines (marked with '!' in column 13)
            if len(line) > 13 and line[13] == '!':
                continue
            # Fixed-width fields:
            #   columns  0-4  : sequential residue number
            #   columns  5-10 : PDB residue number
            #   column  11    : insertion code
            #   column  13    : amino acid (one-letter)
            #   column  16    : secondary structure code
            try:
                res_num  = int(line[5:10].strip())
                res_name = line[13].strip()
                ss_code  = line[16]
                if ss_code == ' ':
                    ss_code = '-'
                residues.append((res_num, res_name, ss_code))
            except (ValueError, IndexError):
                continue
    return residues

# ── extract accession from filename ──────────────────────────────────
def accession_from_path(p):
    """AF-XXXXX-F1-model_v6.cif -> XXXXX"""
    name = Path(p).stem                       # AF-XXXXX-F1-model_v6
    parts = name.split("-")
    # accession is everything between the first 'AF' and 'F1'
    # e.g. AF-A0A087WU62-F1-model_v6 -> A0A087WU62
    f1_idx = parts.index("F1")
    return "-".join(parts[1:f1_idx])

# ── main ─────────────────────────────────────────────────────────────
def main():
    cif_files = sorted(glob.glob(str(PILOT_DIR / "AF-*-F1-model_v*.cif")))
    n_total = len(cif_files)
    print(f"Found {n_total} CIF files in {PILOT_DIR}")

    # Accumulators
    summaries = []        # list of dicts for per-protein summary
    all_residues = []     # list of (acc, res_num, res_name, ss_code) – may be skipped
    failed = []
    collect_residues = True  # we'll estimate size after first 100

    t0 = time.time()

    for idx, cif_path in enumerate(cif_files):
        acc = accession_from_path(cif_path)
        dssp_out = DSSP_DIR / f"{acc}.dssp"

        # ── run mkdssp ────────────────────────────────────────────
        try:
            result = subprocess.run(
                [MKDSSP, "-i", cif_path, "-o", str(dssp_out)],
                capture_output=True, text=True, timeout=TIMEOUT
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr[:300])
        except Exception as exc:
            failed.append((acc, str(exc)[:200]))
            if (idx + 1) % 100 == 0 or idx == 0:
                print(f"  [{idx+1}/{n_total}] FAILED  {acc}: {str(exc)[:80]}")
            continue

        # ── parse DSSP ────────────────────────────────────────────
        try:
            residues = parse_dssp(dssp_out)
        except Exception as exc:
            failed.append((acc, f"parse error: {exc}"))
            continue

        if not residues:
            failed.append((acc, "no residues parsed"))
            continue

        total   = len(residues)
        counts  = Counter(classify(r[2]) for r in residues)
        n_helix  = counts.get("helix", 0)
        n_strand = counts.get("strand", 0)
        n_coil   = counts.get("coil", 0)

        ss_string = "".join(r[2] for r in residues)

        summaries.append({
            "uniprot_accession": acc,
            "total_residues": total,
            "n_helix": n_helix,
            "n_strand": n_strand,
            "n_coil": n_coil,
            "frac_helix":  round(n_helix  / total, 4),
            "frac_strand": round(n_strand / total, 4),
            "frac_coil":   round(n_coil   / total, 4),
            "ss_string": ss_string,
        })

        if collect_residues:
            for res_num, res_name, ss_code in residues:
                all_residues.append((acc, res_num, res_name, ss_code))

        # After 100 proteins, estimate total per-residue file size
        if collect_residues and len(summaries) == 100:
            avg_res = sum(s["total_residues"] for s in summaries) / 100
            # ~40 bytes per line (acc + res_num + res_name + ss_code + tabs + newline)
            est_bytes = avg_res * n_total * 45
            est_mb = est_bytes / 1e6
            print(f"  Estimated per-residue TSV size: {est_mb:.0f} MB")
            if est_mb > 500:
                print("  >> Skipping per-residue TSV (too large). Individual DSSP files retained.")
                collect_residues = False
                all_residues.clear()

        # ── progress ──────────────────────────────────────────────
        if (idx + 1) % 100 == 0 or idx + 1 == n_total:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"  [{idx+1}/{n_total}]  ok={len(summaries)}  fail={len(failed)}  "
                  f"rate={rate:.1f} prot/s  elapsed={elapsed:.0f}s")

    # ── write per-protein summary ─────────────────────────────────
    summary_path = OUT_DIR / "dssp_summary.tsv"
    cols = ["uniprot_accession", "total_residues", "n_helix", "n_strand", "n_coil",
            "frac_helix", "frac_strand", "frac_coil", "ss_string"]
    with open(summary_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in summaries:
            fh.write("\t".join(str(s[c]) for c in cols) + "\n")
    print(f"\nWrote per-protein summary: {summary_path}  ({len(summaries)} rows)")

    # ── write per-residue TSV (if collected) ──────────────────────
    if collect_residues and all_residues:
        per_res_path = OUT_DIR / "dssp_per_residue.tsv"
        with open(per_res_path, "w") as fh:
            fh.write("uniprot_accession\tresidue_number\tresidue_name\tss_code\n")
            for acc, rn, rname, ss in all_residues:
                fh.write(f"{acc}\t{rn}\t{rname}\t{ss}\n")
        size_mb = os.path.getsize(per_res_path) / 1e6
        print(f"Wrote per-residue TSV: {per_res_path}  ({len(all_residues)} rows, {size_mb:.1f} MB)")
    elif not collect_residues:
        print("Per-residue TSV skipped (estimated >500 MB). Individual DSSP files available.")

    # ── failures ──────────────────────────────────────────────────
    if failed:
        fail_path = OUT_DIR / "dssp_failures.tsv"
        with open(fail_path, "w") as fh:
            fh.write("uniprot_accession\terror\n")
            for acc, err in failed:
                fh.write(f"{acc}\t{err}\n")
        print(f"Wrote failure log: {fail_path}  ({len(failed)} failures)")

    # ── summary report ────────────────────────────────────────────
    total_elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("DSSP SUMMARY REPORT")
    print("=" * 70)
    print(f"Total CIF files found : {n_total}")
    print(f"Successfully processed: {len(summaries)}")
    print(f"Failed                : {len(failed)}")
    print(f"Elapsed time          : {total_elapsed:.1f} s")
    print()

    if summaries:
        import statistics
        frac_h = [s["frac_helix"]  for s in summaries]
        frac_e = [s["frac_strand"] for s in summaries]
        frac_c = [s["frac_coil"]   for s in summaries]
        tot_r  = [s["total_residues"] for s in summaries]

        print(f"Mean residues/protein : {statistics.mean(tot_r):.1f}  (median {statistics.median(tot_r):.0f})")
        print()
        print("Secondary structure fractions across all proteins:")
        for label, vals in [("Helix  (H+G+I)", frac_h),
                            ("Strand (E+B)",    frac_e),
                            ("Coil   (T+S+-)",  frac_c)]:
            mn = statistics.mean(vals)
            md = statistics.median(vals)
            sd = statistics.stdev(vals) if len(vals) > 1 else 0
            lo = min(vals)
            hi = max(vals)
            print(f"  {label}:  mean={mn:.3f}  median={md:.3f}  "
                  f"sd={sd:.3f}  min={lo:.3f}  max={hi:.3f}")

        # Distribution bins
        print("\nDistribution of helix fraction:")
        for lo_bin in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            hi_bin = lo_bin + 0.1
            n = sum(1 for v in frac_h if lo_bin <= v < hi_bin)
            bar = "#" * (n // 5) if n > 0 else ""
            print(f"  [{lo_bin:.1f}-{hi_bin:.1f}): {n:5d}  {bar}")

        print("\nDistribution of strand fraction:")
        for lo_bin in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            hi_bin = lo_bin + 0.1
            n = sum(1 for v in frac_e if lo_bin <= v < hi_bin)
            bar = "#" * (n // 5) if n > 0 else ""
            print(f"  [{lo_bin:.1f}-{hi_bin:.1f}): {n:5d}  {bar}")

        print("\nDistribution of coil fraction:")
        for lo_bin in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            hi_bin = lo_bin + 0.1
            n = sum(1 for v in frac_c if lo_bin <= v < hi_bin)
            bar = "#" * (n // 5) if n > 0 else ""
            print(f"  [{lo_bin:.1f}-{hi_bin:.1f}): {n:5d}  {bar}")

    print("\n" + "=" * 70)
    print("Done.")

if __name__ == "__main__":
    main()
