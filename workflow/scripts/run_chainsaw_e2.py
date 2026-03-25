#!/usr/bin/env python3
"""
E2: Run Chainsaw domain segmentation on 239 proteins lacking CATH assignments.

Uses Chainsaw (Wells et al., 2024) — a deep learning model for protein domain
boundary prediction from AlphaFold structures.

Reads: results/domains/cath_protein_summary.tsv (to identify no-CATH proteins)
       results/structures/structure_index.tsv (to find structure files)
Writes: results/domains/chainsaw_domain_predictions.tsv
        results/domains/ml_domain_assignments.tsv (merged CATH + Chainsaw)
        results/domains/chainsaw_report.txt
"""

import subprocess
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
import re

BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
CHAINSAW_DIR = Path("/tmp/chainsaw")
ALPHAFOLD_DIR = BASE / "data/raw/alphafold/pilot"

# Load data
cath_summary = pd.read_csv(BASE / "results/domains/cath_protein_summary.tsv", sep="\t")
struct_idx = pd.read_csv(BASE / "results/structures/structure_index.tsv", sep="\t")

# Identify proteins without CATH
no_cath = cath_summary[cath_summary["has_cath_assignment"] == False]["uniprot_accession"].tolist()
print(f"Proteins without CATH assignment: {len(no_cath)}")

# Check which have structures
struct_map = dict(zip(struct_idx["uniprot_accession"], struct_idx["has_structure"]))
has_struct = [acc for acc in no_cath if struct_map.get(acc, False)]
no_struct = [acc for acc in no_cath if not struct_map.get(acc, False)]
print(f"  With AlphaFold structure: {len(has_struct)}")
print(f"  Without structure: {len(no_struct)}")

# Create temp directory for Chainsaw input — symlink needed structures
import tempfile
tmpdir = tempfile.mkdtemp(prefix="chainsaw_input_")
print(f"Working directory: {tmpdir}")

linked = 0
for acc in has_struct:
    # Find CIF file
    cif_path = ALPHAFOLD_DIR / f"AF-{acc}-F1-model_v4.cif"
    if not cif_path.exists():
        cif_path = ALPHAFOLD_DIR / f"AF-{acc}-F1-model_v6.cif"
    if cif_path.exists():
        os.symlink(cif_path, os.path.join(tmpdir, cif_path.name))
        linked += 1
    else:
        print(f"  WARNING: No CIF found for {acc}")

print(f"Linked {linked} CIF files")

# Run Chainsaw in batches to avoid memory issues
output_file = BASE / "results/domains/chainsaw_raw_output.tsv"

print(f"\nRunning Chainsaw on {linked} structures...")
start_time = time.time()

cmd = [
    sys.executable, str(CHAINSAW_DIR / "get_predictions.py"),
    "--structure_directory", tmpdir,
    "--output", str(output_file),
    "--use_first_chain",
]

result = subprocess.run(
    cmd,
    capture_output=True,
    text=True,
    cwd=str(CHAINSAW_DIR),
    timeout=3600,  # 1 hour max
)

elapsed = time.time() - start_time
print(f"Chainsaw completed in {elapsed:.1f} seconds")

if result.returncode != 0:
    print(f"STDERR (last 2000 chars):\n{result.stderr[-2000:]}")

# Check output
if output_file.exists():
    raw = pd.read_csv(output_file, sep="\t")
    print(f"\nChainsaw output: {len(raw)} rows")
else:
    print("ERROR: No output file produced")
    sys.exit(1)

# Parse results
# chain_id format: AF-{ACC}-F1-model_v{version}
raw["uniprot_accession"] = raw["chain_id"].apply(
    lambda x: x.replace("AF-", "").split("-F1-")[0] if "AF-" in str(x) else x
)

# Replace NULL chopping with empty string
raw["chopping"] = raw["chopping"].replace("NULL", np.nan)

# Parse chopping into domain boundaries
def parse_chopping(chopping_str):
    """Parse Chainsaw chopping string like '1-100,150-300' into domain list."""
    if pd.isna(chopping_str) or chopping_str == "NULL":
        return []
    domains = []
    for domain_str in str(chopping_str).split(","):
        segments = []
        for seg in domain_str.split("_"):
            match = re.match(r"(\d+)-(\d+)", seg.strip())
            if match:
                segments.append((int(match.group(1)), int(match.group(2))))
        if segments:
            start = min(s[0] for s in segments)
            end = max(s[1] for s in segments)
            total_res = sum(s[1] - s[0] + 1 for s in segments)
            domains.append({
                "segments": segments,
                "start": start,
                "end": end,
                "n_residues": total_res,
                "n_segments": len(segments),
            })
    return domains

# Build domain assignment table
domain_rows = []
for _, row in raw.iterrows():
    acc = row["uniprot_accession"]
    domains = parse_chopping(row["chopping"])
    if len(domains) == 0:
        # No domains predicted — single-domain or disordered
        domain_rows.append({
            "uniprot_accession": acc,
            "domain_index": 0,
            "domain_start": 1,
            "domain_end": row["nres"],
            "domain_length": row["nres"],
            "n_segments": 1,
            "chopping": f"1-{row['nres']}",
            "method": "chainsaw_v3",
            "n_domains_total": max(1, row["ndom"]),
            "confidence": row["confidence"],
            "is_single_domain": True,
        })
    else:
        for i, dom in enumerate(domains):
            seg_str = "_".join(f"{s[0]}-{s[1]}" for s in dom["segments"])
            domain_rows.append({
                "uniprot_accession": acc,
                "domain_index": i + 1,
                "domain_start": dom["start"],
                "domain_end": dom["end"],
                "domain_length": dom["n_residues"],
                "n_segments": dom["n_segments"],
                "chopping": seg_str,
                "method": "chainsaw_v3",
                "n_domains_total": row["ndom"],
                "confidence": row["confidence"],
                "is_single_domain": False,
            })

chainsaw_df = pd.DataFrame(domain_rows)

# Save Chainsaw predictions
chainsaw_df.to_csv(BASE / "results/domains/chainsaw_domain_predictions.tsv", sep="\t", index=False)

# Summary stats
n_proteins = chainsaw_df["uniprot_accession"].nunique()
n_with_domains = chainsaw_df[~chainsaw_df["is_single_domain"]]["uniprot_accession"].nunique()
n_single = chainsaw_df[chainsaw_df["is_single_domain"]]["uniprot_accession"].nunique()

print(f"\n=== Chainsaw Results Summary ===")
print(f"Proteins processed: {n_proteins}")
print(f"  Multi-domain: {n_with_domains} ({100*n_with_domains/n_proteins:.1f}%)")
print(f"  Single-domain/no domain: {n_single} ({100*n_single/n_proteins:.1f}%)")

ndom_dist = raw.groupby("ndom").size()
print(f"\nDomain count distribution:")
for ndom, count in sorted(ndom_dist.items()):
    print(f"  {ndom} domains: {count} proteins")

# Mean confidence
print(f"\nMean confidence: {raw['confidence'].mean():.3f}")
print(f"Median confidence: {raw['confidence'].median():.3f}")

# Now merge with existing CATH assignments to create unified domain table
cath_domains = pd.read_csv(BASE / "results/domains/cath_domain_assignments.tsv", sep="\t")

# Create unified assignment table
# CATH proteins keep their assignments; Chainsaw fills the gaps
cath_accs = set(cath_summary[cath_summary["has_cath_assignment"] == True]["uniprot_accession"])
chainsaw_accs = set(chainsaw_df["uniprot_accession"])

print(f"\n=== Unified Domain Assignments ===")
print(f"CATH-assigned: {len(cath_accs)} proteins")
print(f"Chainsaw-assigned: {len(chainsaw_accs)} proteins")
print(f"Overlap: {len(cath_accs & chainsaw_accs)} (should be 0)")
print(f"Total coverage: {len(cath_accs | chainsaw_accs)}/{len(cath_summary)} ({100*len(cath_accs | chainsaw_accs)/len(cath_summary):.1f}%)")

# Build merged summary
merged_rows = []
for _, row in cath_summary.iterrows():
    acc = row["uniprot_accession"]
    if row["has_cath_assignment"]:
        merged_rows.append({
            "uniprot_accession": acc,
            "n_domains": row["n_domains"],
            "domain_source": "CATH/Gene3D",
            "has_domain_assignment": True,
            "domain_architecture": row.get("domain_architecture", ""),
        })
    elif acc in chainsaw_accs:
        cs_row = raw[raw["uniprot_accession"] == acc].iloc[0]
        merged_rows.append({
            "uniprot_accession": acc,
            "n_domains": max(1, cs_row["ndom"]),
            "domain_source": "Chainsaw_v3",
            "has_domain_assignment": True,
            "domain_architecture": cs_row.get("chopping", ""),
        })
    else:
        merged_rows.append({
            "uniprot_accession": acc,
            "n_domains": 0,
            "domain_source": "none",
            "has_domain_assignment": False,
            "domain_architecture": "",
        })

merged_df = pd.DataFrame(merged_rows)
merged_df.to_csv(BASE / "results/domains/ml_domain_assignments.tsv", sep="\t", index=False)

n_assigned = merged_df["has_domain_assignment"].sum()
print(f"Final: {n_assigned}/{len(merged_df)} proteins have domain assignments ({100*n_assigned/len(merged_df):.1f}%)")

# Write report
report = []
report.append("=" * 70)
report.append("E2: Chainsaw Domain Segmentation Report")
report.append("=" * 70)
report.append("")
report.append(f"Tool: Chainsaw v3 (Wells et al., 2024)")
report.append(f"Input: {len(has_struct)} proteins without CATH/Gene3D assignments")
report.append(f"  ({len(no_struct)} additional proteins had no AlphaFold structure)")
report.append(f"Runtime: {elapsed:.1f} seconds ({elapsed/60:.1f} min)")
report.append("")
report.append(f"Results:")
report.append(f"  Proteins processed: {n_proteins}")
report.append(f"  Multi-domain: {n_with_domains} ({100*n_with_domains/n_proteins:.1f}%)")
report.append(f"  Single-domain/no domain: {n_single} ({100*n_single/n_proteins:.1f}%)")
report.append(f"  Mean confidence: {raw['confidence'].mean():.3f}")
report.append("")
report.append("Domain count distribution:")
for ndom, count in sorted(ndom_dist.items()):
    report.append(f"  {ndom} domains: {count} proteins")
report.append("")
report.append("Unified coverage (CATH + Chainsaw):")
report.append(f"  CATH/Gene3D: {len(cath_accs)} proteins")
report.append(f"  Chainsaw: {len(chainsaw_accs)} proteins")
report.append(f"  Total: {n_assigned}/{len(merged_df)} ({100*n_assigned/len(merged_df):.1f}%)")
report.append("")
report.append("Output files:")
report.append(f"  results/domains/chainsaw_domain_predictions.tsv")
report.append(f"  results/domains/ml_domain_assignments.tsv")
report.append("=" * 70)

report_text = "\n".join(report)
(BASE / "results/domains/chainsaw_report.txt").write_text(report_text)
print(f"\n{report_text}")

# Clean up
import shutil
shutil.rmtree(tmpdir, ignore_errors=True)
print("\nDone.")
