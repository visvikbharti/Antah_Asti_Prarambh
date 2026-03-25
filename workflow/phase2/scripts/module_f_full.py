#!/usr/bin/env python3
"""
Module F: Full-scale N-vs-C terminus structural analysis.

Computes structural metrics (contact order, pLDDT, secondary structure,
sequence composition) for N-terminal domain vs C-terminal region across
all proteins with domain assignments.

Scientific basis:
  - Contact order: Plaxco et al. (1998) J Mol Biol 277:985-994
    RCO = (1/NL) * sum(|i-j|) for contacts within 8A, seq_sep >= 6
  - Three-region model: pre-domain tail / N-domain / C-region
  - pLDDT is model confidence (NOT thermodynamic stability)

Requires: Module E output (unified domain assignments + Chainsaw predictions)
"""

import os
import sys
import subprocess
import tempfile
import re
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

# Try gemmi for robust CIF parsing
try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False
    print("WARNING: gemmi not available, using fallback CIF parser")

# ===========================================================================
# Configuration
# ===========================================================================
PROJECT_DIR = os.environ.get("PROJECT_DIR",
    "/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc")
RESULTS = f"{PROJECT_DIR}/results/phase2"
STRIDE_BIN = "/lustre/vishal.bharti/software/chainsaw/stride/stride"

CIF_DIRS = [
    f"{PROJECT_DIR}/data/raw/alphafold/full/ecoli",
    f"{PROJECT_DIR}/data/raw/alphafold/full/human",
    f"{PROJECT_DIR}/data/raw/alphafold/pilot",
]

# Contact order parameters (Plaxco et al. 1998)
CONTACT_DISTANCE_CUTOFF = 8.0  # Angstroms, CA-CA
MIN_SEQUENCE_SEPARATION = 6    # residues

# Amino acid property sets
HYDROPHOBIC = set("AVILFWM")
CHARGED = set("KRDE")
POSITIVE = set("KR")
NEGATIVE = set("DE")
POLAR = set("STNQYHC")
AROMATIC = set("FWY")
SMALL = set("GAS")

# SS code classification (DSSP/Stride convention)
HELIX_CODES = set("HGI")
STRAND_CODES = set("EB")

print("=" * 70)
print("Module F: Full-Scale N-vs-C Terminus Structural Analysis")
print("=" * 70)


# ===========================================================================
# Helper functions
# ===========================================================================

def parse_chainsaw_chopping(chopping_str):
    """Parse Chainsaw chopping notation into domain boundaries.

    Format: comma separates DOMAINS, underscore separates SEGMENTS within domain.
    Example: "12-42_308-396,56-302" -> 2 domains:
        Domain 0: [(12,42), (308,396)]
        Domain 1: [(56,302)]

    Returns list of domains, each domain is list of (start, end) tuples.
    """
    if not chopping_str or chopping_str == "NULL" or pd.isna(chopping_str):
        return []
    domains = []
    for domain_str in chopping_str.split(","):
        segments = []
        for seg_str in domain_str.strip().split("_"):
            match = re.match(r"(\d+)-(\d+)", seg_str.strip())
            if match:
                segments.append((int(match.group(1)), int(match.group(2))))
        if segments:
            domains.append(segments)
    return domains


def find_cif_file(accession):
    """Find AlphaFold CIF file for a protein accession."""
    for d in CIF_DIRS:
        for version in ["v6", "v4"]:
            path = os.path.join(d, f"AF-{accession}-F1-model_{version}.cif")
            if os.path.exists(path):
                return path
    return None


def extract_structure_data_gemmi(cif_path):
    """Read CIF file with gemmi. Returns (ca_coords, plddt, sequence, resnum_map)."""
    st = gemmi.read_structure(cif_path)
    model = st[0]
    chain = model[0]

    ca_coords = []
    plddt = []
    sequence = []
    resnum_map = {}  # resnum -> index in arrays

    idx = 0
    for residue in chain:
        if residue.het_flag == "H":
            continue  # skip heteroatoms
        info = gemmi.find_tabulated_residue(residue.name)
        if not info.is_amino_acid():
            continue
        ca = residue.find_atom("CA", "*")
        if ca is None:
            continue
        ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        plddt.append(ca.b_iso)
        one_letter = info.one_letter_code
        sequence.append(one_letter if one_letter != "?" else "X")
        resnum_map[residue.seqid.num] = idx
        idx += 1

    ca_coords = np.array(ca_coords) if ca_coords else np.zeros((0, 3))
    plddt = np.array(plddt)
    return ca_coords, plddt, sequence, resnum_map


def extract_structure_data_fallback(cif_path):
    """Fallback CIF parser using text processing."""
    ca_coords = []
    plddt = []
    sequence = []
    resnum_map = {}

    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }

    idx = 0
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            atom_name = parts[3]
            if atom_name != "CA":
                continue
            resname = parts[5]
            if resname not in three_to_one:
                continue
            try:
                x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                b = float(parts[14])
                resnum = int(parts[8])
            except (ValueError, IndexError):
                continue
            ca_coords.append([x, y, z])
            plddt.append(b)
            sequence.append(three_to_one.get(resname, "X"))
            resnum_map[resnum] = idx
            idx += 1

    ca_coords = np.array(ca_coords) if ca_coords else np.zeros((0, 3))
    plddt = np.array(plddt)
    return ca_coords, plddt, sequence, resnum_map


def extract_structure_data(cif_path):
    """Read CIF file. Returns (ca_coords, plddt, sequence, resnum_map)."""
    if HAS_GEMMI:
        return extract_structure_data_gemmi(cif_path)
    return extract_structure_data_fallback(cif_path)


def run_stride(cif_path, tmpdir):
    """Run Stride for secondary structure. Returns dict: resnum -> SS code."""
    if not os.path.exists(STRIDE_BIN):
        return {}

    ss_map = {}
    try:
        # Convert CIF to PDB using gemmi
        if HAS_GEMMI:
            st = gemmi.read_structure(cif_path)
            pdb_path = os.path.join(tmpdir, "temp.pdb")
            st.write_pdb(pdb_path)
        else:
            return {}

        result = subprocess.run(
            [STRIDE_BIN, pdb_path],
            capture_output=True, text=True, timeout=30
        )

        for line in result.stdout.split("\n"):
            if line.startswith("ASG"):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        resnum = int(parts[3])
                        ss_code = parts[5]
                        ss_map[resnum] = ss_code
                    except (ValueError, IndexError):
                        pass
    except Exception:
        pass

    return ss_map


def compute_contact_order(ca_coords, residue_indices, min_sep=MIN_SEQUENCE_SEPARATION,
                          cutoff=CONTACT_DISTANCE_CUTOFF):
    """Compute absolute and relative contact order for a set of residues.

    Uses Plaxco et al. (1998) formula:
      RCO = (1/(N*L)) * sum(|i-j|) for contacts where d(CA_i, CA_j) < cutoff
      and |i-j| >= min_sep

    Args:
        ca_coords: Nx3 array of all CA coordinates
        residue_indices: list of 0-based indices into ca_coords for this region
        min_sep: minimum sequence separation (default 6)
        cutoff: distance cutoff in Angstroms (default 8.0)

    Returns: (absolute_co, relative_co, n_contacts)
    """
    if len(residue_indices) < min_sep + 1:
        return np.nan, np.nan, 0

    indices = np.array(sorted(residue_indices))
    coords = ca_coords[indices]
    n = len(indices)

    # Compute pairwise distances (upper triangle only)
    total_sep = 0.0
    n_contacts = 0
    for i in range(n):
        for j in range(i + 1, n):
            seq_sep = abs(indices[j] - indices[i])
            if seq_sep < min_sep:
                continue
            dx = coords[j] - coords[i]
            dist = np.sqrt(np.sum(dx * dx))
            if dist < cutoff:
                total_sep += seq_sep
                n_contacts += 1

    if n_contacts == 0:
        return np.nan, np.nan, 0

    aco = total_sep / n_contacts
    rco = aco / n
    return aco, rco, n_contacts


def compute_region_metrics(sequence_chars, ss_codes, plddt_values,
                           ca_coords, residue_indices):
    """Compute all metrics for a protein region.

    Returns dict with sequence, structure, and contact order metrics.
    """
    n = len(residue_indices)
    if n == 0:
        return None

    # Sequence metrics
    seq = [sequence_chars[i] for i in residue_indices if i < len(sequence_chars)]
    n_seq = len(seq)
    if n_seq == 0:
        return None

    n_pos = sum(1 for aa in seq if aa in POSITIVE)
    n_neg = sum(1 for aa in seq if aa in NEGATIVE)

    metrics = {
        "length": n_seq,
        "net_charge": n_pos - n_neg,
        "frac_hydrophobic": sum(1 for aa in seq if aa in HYDROPHOBIC) / n_seq,
        "frac_charged": sum(1 for aa in seq if aa in CHARGED) / n_seq,
        "frac_polar": sum(1 for aa in seq if aa in POLAR) / n_seq,
        "frac_aromatic": sum(1 for aa in seq if aa in AROMATIC) / n_seq,
        "frac_small": sum(1 for aa in seq if aa in SMALL) / n_seq,
    }

    # pLDDT metrics
    plddt_region = plddt_values[residue_indices] if len(plddt_values) > 0 else np.array([])
    if len(plddt_region) > 0:
        metrics["mean_plddt"] = float(np.mean(plddt_region))
        metrics["frac_plddt_gt70"] = float(np.mean(plddt_region > 70))
        metrics["frac_plddt_gt90"] = float(np.mean(plddt_region > 90))
    else:
        metrics["mean_plddt"] = np.nan
        metrics["frac_plddt_gt70"] = np.nan
        metrics["frac_plddt_gt90"] = np.nan

    # Secondary structure metrics
    if ss_codes:
        ss_region = [ss_codes.get(i, "C") for i in residue_indices]
        n_ss = len(ss_region)
        metrics["frac_helix"] = sum(1 for s in ss_region if s in HELIX_CODES) / n_ss
        metrics["frac_strand"] = sum(1 for s in ss_region if s in STRAND_CODES) / n_ss
        metrics["frac_coil"] = 1.0 - metrics["frac_helix"] - metrics["frac_strand"]
    else:
        metrics["frac_helix"] = np.nan
        metrics["frac_strand"] = np.nan
        metrics["frac_coil"] = np.nan

    # Contact order
    aco, rco, n_contacts = compute_contact_order(ca_coords, residue_indices)
    metrics["abs_contact_order"] = aco
    metrics["relative_contact_order"] = rco
    metrics["n_contacts"] = n_contacts

    return metrics


# ===========================================================================
# Load data
# ===========================================================================
print("\n--- Loading data ---")

# Unified domain assignments from Module E
unified = pd.read_csv(f"{RESULTS}/domains/unified_domain_assignments_full.tsv", sep="\t")
print(f"Unified domain assignments: {len(unified)} proteins")

# Chainsaw predictions (for boundary parsing)
chainsaw = pd.read_csv(f"{RESULTS}/domains/chainsaw_full_predictions.tsv", sep="\t")
# Extract accession from chain_id
if "accession" not in chainsaw.columns:
    chainsaw["accession"] = chainsaw["chain_id"].str.extract(
        r"AF-([A-Z0-9]+)-", expand=False)
print(f"Chainsaw predictions: {len(chainsaw)} proteins")

# Phase 1 CATH domain assignments (for boundary info)
cath_path = f"{PROJECT_DIR}/results/domains/cath_domain_assignments.tsv"
if os.path.exists(cath_path):
    cath = pd.read_csv(cath_path, sep="\t")
    # Use correct column name
    cath_acc_col = "uniprot_accession" if "uniprot_accession" in cath.columns else "accession"
    cath_start_col = "domain_start" if "domain_start" in cath.columns else "start"
    cath_end_col = "domain_end" if "domain_end" in cath.columns else "end"
    print(f"Phase 1 CATH domains: {len(cath)} records (acc col: {cath_acc_col})")
else:
    cath = pd.DataFrame()
    print("No Phase 1 CATH data found")

# Substrate lists (correct column names)
groel = pd.read_csv(f"{PROJECT_DIR}/data/processed/groel_substrates_standardized.tsv", sep="\t")
hsp60 = pd.read_csv(f"{PROJECT_DIR}/data/processed/hsp60_tier1_substrates.tsv", sep="\t")
matrix = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_matrix_proteome.tsv", sep="\t")
mito = pd.read_csv(f"{PROJECT_DIR}/data/processed/human_mito_proteome.tsv", sep="\t")

groel_acc_col = "current_accession" if "current_accession" in groel.columns else "accession"
hsp60_acc_col = "uniprot_id" if "uniprot_id" in hsp60.columns else "accession"
matrix_acc_col = "uniprot_id" if "uniprot_id" in matrix.columns else "accession"
mito_acc_col = "uniprot_id" if "uniprot_id" in mito.columns else "accession"

groel_acc = set(groel[groel_acc_col].dropna().values)
hsp60_acc = set(hsp60[hsp60_acc_col].dropna().values)
matrix_acc = set(matrix[matrix_acc_col].dropna().values)
mito_acc = set(mito[mito_acc_col].dropna().values)

print(f"Substrates: GroEL={len(groel_acc)}, HSP60={len(hsp60_acc)}")
print(f"Backgrounds: Matrix={len(matrix_acc)}, Mito={len(mito_acc)}")

# GroEL class lookup
groel_class = {}
class_col = "groel_class" if "groel_class" in groel.columns else None
if class_col:
    for _, row in groel.iterrows():
        groel_class[row[groel_acc_col]] = row[class_col]

# Homolog pairs
homolog_path = f"{PROJECT_DIR}/data/processed/groel_hsp60_homologs.tsv"
if os.path.exists(homolog_path):
    homologs = pd.read_csv(homolog_path, sep="\t")
    homolog_accs = set()
    for col in ["groel_accession", "hsp60_accession"]:
        if col in homologs.columns:
            homolog_accs.update(homologs[col].dropna().values)
    print(f"Homolog pairs: {len(homologs)}, unique proteins: {len(homolog_accs)}")
else:
    homolog_accs = set()

# E. coli proteome accessions (for background)
ecoli_path = f"{PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.tsv"
if os.path.exists(ecoli_path):
    ecoli = pd.read_csv(ecoli_path, sep="\t")
    ecoli_acc_col = "Entry" if "Entry" in ecoli.columns else ecoli.columns[0]
    ecoli_acc = set(ecoli[ecoli_acc_col].dropna().values)
    print(f"E. coli proteome: {len(ecoli_acc)}")
else:
    ecoli_acc = set()


# ===========================================================================
# Parse domain boundaries
# ===========================================================================
print("\n--- Parsing domain boundaries ---")

# Build boundary lookup: accession -> {source, n_domains, domains: [[(start,end),...], ...]}
boundary_lookup = {}

# 1. From CATH (preferred)
if not cath.empty:
    for acc, group in cath.groupby(cath_acc_col):
        domains = []
        for _, row in group.iterrows():
            try:
                s = int(row[cath_start_col])
                e = int(row[cath_end_col])
                domains.append([(s, e)])  # each CATH domain is contiguous
            except (ValueError, KeyError):
                pass
        if domains:
            domains.sort(key=lambda d: d[0][0])  # sort by start of first segment
            boundary_lookup[acc] = {
                "source": "CATH",
                "n_domains": len(domains),
                "domains": domains,
            }
    print(f"  CATH boundaries: {len(boundary_lookup)} proteins")

# 2. From Chainsaw (for proteins not in CATH)
n_chainsaw = 0
for _, row in chainsaw.iterrows():
    acc = row.get("accession")
    if pd.isna(acc) or acc in boundary_lookup:
        continue
    chopping = str(row.get("chopping", ""))
    ndom = int(row.get("ndom", 0))
    if ndom < 1 or chopping == "NULL" or not chopping:
        continue
    domains = parse_chainsaw_chopping(chopping)
    if domains:
        domains.sort(key=lambda d: d[0][0])
        boundary_lookup[acc] = {
            "source": "Chainsaw",
            "n_domains": len(domains),
            "domains": domains,
        }
        n_chainsaw += 1
print(f"  Chainsaw boundaries: {n_chainsaw} proteins")
print(f"  Total with boundaries: {len(boundary_lookup)} proteins")

multi_domain_accs = {acc for acc, info in boundary_lookup.items()
                     if info["n_domains"] >= 2}
print(f"  Multi-domain: {len(multi_domain_accs)}")


# ===========================================================================
# Identify proteins to process
# ===========================================================================
print("\n--- Identifying proteins to process ---")

# Priority: substrates, backgrounds, then homolog partners
target_accs = set()
target_accs.update(groel_acc)
target_accs.update(hsp60_acc)
target_accs.update(mito_acc)
target_accs.update(ecoli_acc)
target_accs.update(homolog_accs)

# Filter to proteins that have domain boundaries
processable = target_accs & set(boundary_lookup.keys())
print(f"  Target proteins: {len(target_accs)}")
print(f"  With domain boundaries: {len(processable)}")

# Focus on multi-domain for N-vs-C analysis (but compute all for background)
multi_to_process = processable & multi_domain_accs
print(f"  Multi-domain to process: {len(multi_to_process)}")


# ===========================================================================
# Process structures
# ===========================================================================
print("\n--- Processing structures ---")

region_records = []
contact_order_records = []
structure_metric_records = []
sequence_metric_records = []
paired_records = []

processed = 0
skipped_no_cif = 0
skipped_error = 0

# Create temp directory for Stride
tmpdir = tempfile.mkdtemp(prefix="module_f_", dir=f"{PROJECT_DIR}/tmp")

for acc in sorted(processable):
    info = boundary_lookup[acc]
    domains = info["domains"]
    n_dom = info["n_domains"]
    source = info["source"]

    # Find CIF file
    cif_path = find_cif_file(acc)
    if cif_path is None:
        skipped_no_cif += 1
        continue

    try:
        # Extract structural data
        ca_coords, plddt_arr, seq_chars, resnum_map = extract_structure_data(cif_path)
        if len(ca_coords) < 10:
            skipped_error += 1
            continue

        protein_length = len(ca_coords)

        # Run Stride for secondary structure
        ss_map_resnum = run_stride(cif_path, tmpdir)
        # Convert resnum-based SS to index-based
        ss_map_idx = {}
        for resnum, ss_code in ss_map_resnum.items():
            if resnum in resnum_map:
                ss_map_idx[resnum_map[resnum]] = ss_code

        # Define first domain indices
        first_domain = domains[0]
        first_domain_indices = []
        for seg_start, seg_end in first_domain:
            for resnum in range(seg_start, seg_end + 1):
                if resnum in resnum_map:
                    first_domain_indices.append(resnum_map[resnum])

        if not first_domain_indices:
            skipped_error += 1
            continue

        first_domain_start = first_domain[0][0]
        first_domain_end = max(seg_end for _, seg_end in first_domain)

        # Pre-domain tail: indices before first domain
        pre_domain_indices = [resnum_map[rn] for rn in resnum_map
                              if rn < first_domain_start]

        # C-region: indices after last segment of first domain
        c_region_indices = [resnum_map[rn] for rn in resnum_map
                            if rn > first_domain_end]

        # Dataset assignment
        datasets = []
        if acc in groel_acc: datasets.append("groel")
        if acc in hsp60_acc: datasets.append("hsp60")
        if acc in matrix_acc: datasets.append("matrix_bg")
        if acc in mito_acc: datasets.append("mito_bg")
        if not datasets: datasets.append("proteome_bg")

        # Region boundary record
        region_records.append({
            "accession": acc,
            "source": source,
            "protein_length": protein_length,
            "n_domains": n_dom,
            "is_multi_domain": n_dom >= 2,
            "first_domain_start": first_domain_start,
            "first_domain_end": first_domain_end,
            "pre_domain_length": len(pre_domain_indices),
            "n_domain_length": len(first_domain_indices),
            "c_region_length": len(c_region_indices),
            "datasets": ",".join(datasets),
            "groel_class": groel_class.get(acc, ""),
        })

        # Compute metrics per region
        for region_name, region_indices in [
            ("pre_domain", pre_domain_indices),
            ("n_domain", first_domain_indices),
            ("c_region", c_region_indices),
        ]:
            if len(region_indices) < 5:
                continue

            m = compute_region_metrics(
                seq_chars, ss_map_idx, plddt_arr, ca_coords,
                np.array(sorted(region_indices))
            )
            if m is None:
                continue

            # Contact order record
            contact_order_records.append({
                "accession": acc,
                "region": region_name,
                "n_residues": m["length"],
                "n_contacts": m["n_contacts"],
                "abs_contact_order": m["abs_contact_order"],
                "relative_contact_order": m["relative_contact_order"],
            })

            # Structure metric record
            structure_metric_records.append({
                "accession": acc,
                "region": region_name,
                "frac_helix": m["frac_helix"],
                "frac_strand": m["frac_strand"],
                "frac_coil": m["frac_coil"],
                "mean_plddt": m["mean_plddt"],
                "frac_plddt_gt70": m["frac_plddt_gt70"],
                "frac_plddt_gt90": m["frac_plddt_gt90"],
            })

            # Sequence metric record
            sequence_metric_records.append({
                "accession": acc,
                "region": region_name,
                "length": m["length"],
                "net_charge": m["net_charge"],
                "frac_hydrophobic": m["frac_hydrophobic"],
                "frac_charged": m["frac_charged"],
                "frac_polar": m["frac_polar"],
                "frac_aromatic": m["frac_aromatic"],
                "frac_small": m["frac_small"],
            })

        # Paired N-vs-C comparison (multi-domain only)
        if n_dom >= 2 and len(first_domain_indices) >= 5 and len(c_region_indices) >= 5:
            n_metrics = compute_region_metrics(
                seq_chars, ss_map_idx, plddt_arr, ca_coords,
                np.array(sorted(first_domain_indices))
            )
            c_metrics = compute_region_metrics(
                seq_chars, ss_map_idx, plddt_arr, ca_coords,
                np.array(sorted(c_region_indices))
            )

            if n_metrics and c_metrics:
                pair = {
                    "accession": acc,
                    "source": source,
                    "n_domains": n_dom,
                    "datasets": ",".join(datasets),
                    "groel_class": groel_class.get(acc, ""),
                    "protein_length": protein_length,
                }
                for metric in ["length", "net_charge", "frac_hydrophobic",
                               "frac_helix", "frac_strand", "mean_plddt",
                               "frac_plddt_gt70", "relative_contact_order"]:
                    n_val = n_metrics.get(metric, np.nan)
                    c_val = c_metrics.get(metric, np.nan)
                    pair[f"{metric}_n_domain"] = n_val
                    pair[f"{metric}_c_region"] = c_val
                    if not np.isnan(n_val) and not np.isnan(c_val):
                        pair[f"{metric}_diff"] = n_val - c_val
                    else:
                        pair[f"{metric}_diff"] = np.nan

                paired_records.append(pair)

        processed += 1
        if processed % 500 == 0:
            print(f"  Processed {processed} proteins... ({skipped_no_cif} no CIF, {skipped_error} errors)")

    except Exception as e:
        skipped_error += 1
        if skipped_error <= 5:
            print(f"  Error processing {acc}: {e}")

# Cleanup
try:
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
except Exception:
    pass

print(f"\n  Total processed: {processed}")
print(f"  Skipped (no CIF): {skipped_no_cif}")
print(f"  Skipped (error): {skipped_error}")
print(f"  Paired records: {len(paired_records)}")


# ===========================================================================
# Save outputs
# ===========================================================================
print("\n--- Saving outputs ---")
os.makedirs(f"{RESULTS}/stability", exist_ok=True)

# Region boundaries
if region_records:
    regions_df = pd.DataFrame(region_records)
    path = f"{RESULTS}/stability/region_boundaries_full.tsv"
    regions_df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path} ({len(regions_df)} proteins)")

# Contact order
if contact_order_records:
    co_df = pd.DataFrame(contact_order_records)
    path = f"{RESULTS}/stability/contact_order_full.tsv"
    co_df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path} ({len(co_df)} records)")

# Structure metrics
if structure_metric_records:
    sm_df = pd.DataFrame(structure_metric_records)
    path = f"{RESULTS}/stability/structure_metrics_full.tsv"
    sm_df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path} ({len(sm_df)} records)")

# Sequence metrics
if sequence_metric_records:
    sq_df = pd.DataFrame(sequence_metric_records)
    path = f"{RESULTS}/stability/sequence_metrics_full.tsv"
    sq_df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path} ({len(sq_df)} records)")

# Paired N-vs-C
if paired_records:
    paired_df = pd.DataFrame(paired_records)
    path = f"{RESULTS}/stability/n_vs_c_paired_full.tsv"
    paired_df.to_csv(path, sep="\t", index=False)
    print(f"Saved: {path} ({len(paired_df)} paired comparisons)")

    # Also save as stability comparisons for Module H
    os.makedirs(f"{RESULTS}/stats", exist_ok=True)
    paired_df.to_csv(f"{RESULTS}/stats/stability_comparisons_full.tsv",
                     sep="\t", index=False)

# Quick summary statistics
if paired_records:
    paired_df = pd.DataFrame(paired_records)
    print("\n--- Quick N-vs-C Summary ---")
    for ds in ["groel", "hsp60", "matrix_bg", "mito_bg"]:
        sub = paired_df[paired_df["datasets"].str.contains(ds)]
        if len(sub) < 5:
            continue
        rco_diff = sub["relative_contact_order_diff"].dropna()
        plddt_diff = sub["mean_plddt_diff"].dropna()
        print(f"\n{ds} (n={len(sub)}):")
        if len(rco_diff) > 0:
            print(f"  RCO diff (N-C): mean={rco_diff.mean():.4f}, "
                  f"median={rco_diff.median():.4f}")
        if len(plddt_diff) > 0:
            print(f"  pLDDT diff (N-C): mean={plddt_diff.mean():.2f}, "
                  f"median={plddt_diff.median():.2f}")

print("\n" + "=" * 70)
print("Module F complete.")
print("=" * 70)
