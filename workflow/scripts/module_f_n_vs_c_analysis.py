#!/usr/bin/env python3
"""
Module F: N-domain vs C-region analysis for the Antah Asti Prarambh project.

Steps F1-F5: Define three structural regions per protein (pre-domain tail,
first CATH domain = N-domain, remainder = C-region), compute sequence-derived
and structure-derived metrics per region, compute contact order, and perform
within-protein paired comparisons.

Scientific definitions:
  - Pre-domain tail: residues before the first CATH domain starts
  - N-domain: the first CATH domain by residue position
  - C-region: everything after the first domain ends
  - Single-domain proteins: compare N-half vs C-half of the domain
  - pLDDT is a CONFIDENCE metric, not a stability metric

Only the 1,151 proteins WITH CATH domain assignments are analyzed.
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path("/Users/vishalbharti/Downloads/Antah_Asti_Prarambh")
DOMAIN_ASSIGN  = BASE / "results/domains/cath_domain_assignments.tsv"
PROTEIN_SUMM   = BASE / "results/domains/cath_protein_summary.tsv"
DSSP_FILE      = BASE / "results/structures/dssp_per_residue.tsv"
STRUCT_INDEX   = BASE / "results/structures/structure_index.tsv"
GROEL_FILE     = BASE / "data/processed/groel_substrates_standardized.tsv"
CIF_DIR        = BASE / "data/raw/alphafold/pilot"
OUT_DIR        = BASE / "results/termini"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Three-letter to one-letter amino acid mapping
AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Amino acid property sets (one-letter codes)
HYDROPHOBIC = set('AVILFWM')
CHARGED     = set('KRDE')
POLAR       = set('STNQYHC')
AROMATIC    = set('FWY')
SMALL       = set('GAS')
POSITIVE    = set('KR')
NEGATIVE    = set('DE')


def parse_cif_ca_atoms(cif_path):
    """
    Parse CA atom records from an AlphaFold CIF file.
    Returns a dict: residue_number -> {
        'comp_id': str (3-letter), 'x': float, 'y': float, 'z': float,
        'bfactor': float (pLDDT)
    }
    Text-based parser, one protein at a time for memory efficiency.
    """
    ca_atoms = {}
    in_atom_section = False
    with open(cif_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                # parts[3] = label_atom_id, parts[5] = label_comp_id
                # parts[8] = label_seq_id, parts[10,11,12] = x,y,z
                # parts[14] = B_iso_or_equiv (pLDDT)
                if parts[3] == 'CA':
                    try:
                        resnum = int(parts[8])
                        ca_atoms[resnum] = {
                            'comp_id': parts[5],
                            'x': float(parts[10]),
                            'y': float(parts[11]),
                            'z': float(parts[12]),
                            'bfactor': float(parts[14])
                        }
                    except (ValueError, IndexError):
                        continue
    return ca_atoms


def compute_net_charge(residues_1letter):
    """Net charge at pH 7: count(K+R) - count(D+E)."""
    pos = sum(1 for r in residues_1letter if r in POSITIVE)
    neg = sum(1 for r in residues_1letter if r in NEGATIVE)
    return pos - neg


def compute_fractions(residues_1letter):
    """Compute amino acid composition fractions."""
    n = len(residues_1letter)
    if n == 0:
        return {'frac_hydrophobic': np.nan, 'frac_charged': np.nan,
                'frac_polar': np.nan, 'frac_aromatic': np.nan, 'frac_small': np.nan}
    return {
        'frac_hydrophobic': sum(1 for r in residues_1letter if r in HYDROPHOBIC) / n,
        'frac_charged':     sum(1 for r in residues_1letter if r in CHARGED) / n,
        'frac_polar':       sum(1 for r in residues_1letter if r in POLAR) / n,
        'frac_aromatic':    sum(1 for r in residues_1letter if r in AROMATIC) / n,
        'frac_small':       sum(1 for r in residues_1letter if r in SMALL) / n,
    }


def compute_ss_fractions(ss_codes):
    """Compute secondary structure fractions from DSSP codes."""
    n = len(ss_codes)
    if n == 0:
        return {'frac_helix': np.nan, 'frac_strand': np.nan, 'frac_coil': np.nan}
    helix_codes = {'H', 'G', 'I'}   # alpha, 3-10, pi helix
    strand_codes = {'E', 'B'}        # extended strand, isolated bridge
    n_helix  = sum(1 for s in ss_codes if s in helix_codes)
    n_strand = sum(1 for s in ss_codes if s in strand_codes)
    n_coil   = n - n_helix - n_strand
    return {
        'frac_helix':  n_helix / n,
        'frac_strand': n_strand / n,
        'frac_coil':   n_coil / n,
    }


def compute_plddt_metrics(plddt_values):
    """Compute pLDDT confidence metrics from B-factor values."""
    if len(plddt_values) == 0:
        return {'mean_plddt': np.nan, 'frac_plddt_gt70': np.nan, 'frac_plddt_gt90': np.nan}
    arr = np.array(plddt_values)
    return {
        'mean_plddt':      float(np.mean(arr)),
        'frac_plddt_gt70': float(np.mean(arr > 70)),
        'frac_plddt_gt90': float(np.mean(arr > 90)),
    }


def compute_contact_order(coords, residue_numbers, min_seq_sep=6, distance_cutoff=8.0):
    """
    Compute relative contact order (Plaxco-Simons style).
    coords: Nx3 array of CA coordinates
    residue_numbers: array of residue numbers
    Returns: (n_residues, n_contacts, abs_contact_order, relative_contact_order)
    """
    n = len(coords)
    if n < 2:
        return n, 0, np.nan, np.nan

    coords = np.array(coords, dtype=np.float64)
    resnums = np.array(residue_numbers, dtype=np.int64)

    # Compute pairwise distances (vectorized for up to ~2000 residues)
    # For memory safety, process in chunks if needed
    if n > 3000:
        # Chunk-based approach for very large proteins
        total_sep = 0.0
        n_contacts = 0
        chunk_size = 1000
        for i_start in range(0, n, chunk_size):
            i_end = min(i_start + chunk_size, n)
            for j_start in range(i_start, n, chunk_size):
                j_end = min(j_start + chunk_size, n)
                diff = coords[i_start:i_end, np.newaxis, :] - coords[np.newaxis, j_start:j_end, :]
                dists = np.sqrt(np.sum(diff**2, axis=2))
                seq_sep = np.abs(resnums[i_start:i_end, np.newaxis] - resnums[np.newaxis, j_start:j_end])
                mask = (dists < distance_cutoff) & (seq_sep >= min_seq_sep)
                # Avoid double counting: only count i < j
                for ii in range(i_end - i_start):
                    for jj in range(j_end - j_start):
                        gi = i_start + ii
                        gj = j_start + jj
                        if gi < gj and mask[ii, jj]:
                            total_sep += seq_sep[ii, jj]
                            n_contacts += 1
        if n_contacts == 0:
            return n, 0, np.nan, np.nan
        aco = total_sep / n_contacts
        rco = aco / n
        return n, n_contacts, float(aco), float(rco)
    else:
        # Vectorized approach
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dists = np.sqrt(np.sum(diff**2, axis=2))
        seq_sep = np.abs(resnums[:, np.newaxis] - resnums[np.newaxis, :])

        # Upper triangle only (i < j)
        mask = np.triu(np.ones((n, n), dtype=bool), k=1)
        contact_mask = mask & (dists < distance_cutoff) & (seq_sep >= min_seq_sep)

        n_contacts = int(np.sum(contact_mask))
        if n_contacts == 0:
            return n, 0, np.nan, np.nan

        total_sep = float(np.sum(seq_sep[contact_mask]))
        aco = total_sep / n_contacts
        rco = aco / n
        return n, n_contacts, float(aco), float(rco)


# ═══════════════════════════════════════════════════════════════════════════
# STEP F1: Define three regions per protein
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("STEP F1: Define three regions per protein")
print("=" * 72)

# Load data
domain_df = pd.read_csv(DOMAIN_ASSIGN, sep='\t')
summary_df = pd.read_csv(PROTEIN_SUMM, sep='\t')
struct_df = pd.read_csv(STRUCT_INDEX, sep='\t')

# Only proteins WITH CATH assignments
cath_proteins = summary_df[summary_df['has_cath_assignment'] == True].copy()
print(f"Proteins with CATH assignments: {len(cath_proteins)}")

# Source dataset lookup
source_lookup = dict(zip(struct_df['uniprot_accession'], struct_df['source_dataset']))

# Build region boundaries
region_rows = []
for _, prow in cath_proteins.iterrows():
    acc = prow['uniprot_accession']
    prot_len = prow['protein_length']
    n_dom = prow['n_domains']

    # Get domains for this protein, sorted by start position
    prot_domains = domain_df[domain_df['uniprot_accession'] == acc].copy()
    prot_domains = prot_domains.sort_values('domain_start')

    if len(prot_domains) == 0:
        continue

    first_dom = prot_domains.iloc[0]
    first_start = int(first_dom['domain_start'])
    first_end = int(first_dom['domain_end'])

    # Pre-domain tail
    pre_start = 1
    pre_end = first_start - 1
    pre_len = max(0, pre_end - pre_start + 1)
    if pre_end < pre_start:
        pre_end = 0
        pre_len = 0

    # N-domain
    n_dom_start = first_start
    n_dom_end = first_end
    n_dom_len = n_dom_end - n_dom_start + 1

    # C-region
    c_start = first_end + 1
    c_end = prot_len
    c_len = max(0, c_end - c_start + 1)
    if c_start > prot_len:
        c_start = 0
        c_end = 0
        c_len = 0

    is_single = (n_dom == 1)
    src = source_lookup.get(acc, '')

    region_rows.append({
        'uniprot_accession': acc,
        'protein_length': prot_len,
        'n_domains': n_dom,
        'is_single_domain': is_single,
        'pre_domain_start': pre_start if pre_len > 0 else 0,
        'pre_domain_end': pre_end,
        'pre_domain_length': pre_len,
        'n_domain_start': n_dom_start,
        'n_domain_end': n_dom_end,
        'n_domain_length': n_dom_len,
        'c_region_start': c_start,
        'c_region_end': c_end,
        'c_region_length': c_len,
        'source_dataset': src,
    })

region_df = pd.DataFrame(region_rows)
region_df.to_csv(OUT_DIR / "region_boundaries.tsv", sep='\t', index=False)
print(f"Region boundaries saved: {len(region_df)} proteins")
print(f"  Single-domain: {region_df['is_single_domain'].sum()}")
print(f"  Multi-domain:  {(~region_df['is_single_domain']).sum()}")
print(f"  Mean pre-domain length: {region_df['pre_domain_length'].mean():.1f}")
print(f"  Mean N-domain length:   {region_df['n_domain_length'].mean():.1f}")
print(f"  Mean C-region length:   {region_df['c_region_length'].mean():.1f}")
print()

# ═══════════════════════════════════════════════════════════════════════════
# STEP F2: Sequence-derived metrics per region
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("STEP F2: Sequence-derived metrics per region")
print("=" * 72)

# Load DSSP per-residue (has residue_name and residue_number)
print("Loading DSSP per-residue data...")
dssp_df = pd.read_csv(DSSP_FILE, sep='\t')
print(f"  Total DSSP records: {len(dssp_df)}")

# Group by protein for faster access
dssp_grouped = {acc: grp for acc, grp in dssp_df.groupby('uniprot_accession')}

seq_rows = []
skipped = 0

for _, rrow in region_df.iterrows():
    acc = rrow['uniprot_accession']
    if acc not in dssp_grouped:
        skipped += 1
        continue

    pdssp = dssp_grouped[acc]
    resnum_to_name = dict(zip(pdssp['residue_number'], pdssp['residue_name']))

    # Define regions to process
    regions = []

    # Pre-domain (if exists)
    if rrow['pre_domain_length'] > 0:
        regions.append(('pre_domain', int(rrow['pre_domain_start']), int(rrow['pre_domain_end'])))

    # N-domain
    regions.append(('n_domain', int(rrow['n_domain_start']), int(rrow['n_domain_end'])))

    # C-region (if exists)
    if rrow['c_region_length'] > 0:
        regions.append(('c_region', int(rrow['c_region_start']), int(rrow['c_region_end'])))

    # For single-domain proteins, also split the domain at midpoint
    if rrow['is_single_domain']:
        mid = (rrow['n_domain_start'] + rrow['n_domain_end']) // 2
        regions.append(('n_half', int(rrow['n_domain_start']), int(mid)))
        regions.append(('c_half', int(mid + 1), int(rrow['n_domain_end'])))

    for region_name, rstart, rend in regions:
        residues_1letter = []
        for rn in range(rstart, rend + 1):
            res_name = resnum_to_name.get(rn, '')
            # DSSP file uses 1-letter codes; CIF uses 3-letter. Handle both.
            if len(res_name) == 1 and res_name in 'ACDEFGHIKLMNPQRSTVWY':
                residues_1letter.append(res_name)
            elif len(res_name) == 3:
                aa1 = AA3_TO_1.get(res_name, '')
                if aa1:
                    residues_1letter.append(aa1)

        length = len(residues_1letter)
        net_charge = compute_net_charge(residues_1letter) if length > 0 else np.nan
        fracs = compute_fractions(residues_1letter)

        seq_rows.append({
            'uniprot_accession': acc,
            'region': region_name,
            'length': length,
            'net_charge': net_charge,
            **fracs,
        })

seq_df = pd.DataFrame(seq_rows)
seq_df.to_csv(OUT_DIR / "sequence_metrics.tsv", sep='\t', index=False)
print(f"Sequence metrics saved: {len(seq_df)} region entries")
if skipped > 0:
    print(f"  Skipped {skipped} proteins (no DSSP data)")
print()

# ═══════════════════════════════════════════════════════════════════════════
# STEP F3: Structure-derived metrics per region
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("STEP F3: Structure-derived metrics per region (DSSP + pLDDT)")
print("=" * 72)

struct_rows = []
cif_missing = 0
cif_parsed = 0

for idx, rrow in region_df.iterrows():
    acc = rrow['uniprot_accession']

    # Get DSSP data
    pdssp = dssp_grouped.get(acc)
    if pdssp is None:
        continue

    resnum_to_ss = dict(zip(pdssp['residue_number'], pdssp['ss_code']))

    # Parse CIF for pLDDT (CA atoms)
    cif_path = CIF_DIR / f"AF-{acc}-F1-model_v6.cif"
    if not cif_path.exists():
        # Try v4 or other versions
        cif_candidates = list(CIF_DIR.glob(f"AF-{acc}-F1-model_v*.cif"))
        if cif_candidates:
            cif_path = cif_candidates[0]
        else:
            cif_missing += 1
            continue

    ca_atoms = parse_cif_ca_atoms(cif_path)
    cif_parsed += 1

    # Define regions
    regions = []
    if rrow['pre_domain_length'] > 0:
        regions.append(('pre_domain', int(rrow['pre_domain_start']), int(rrow['pre_domain_end'])))
    regions.append(('n_domain', int(rrow['n_domain_start']), int(rrow['n_domain_end'])))
    if rrow['c_region_length'] > 0:
        regions.append(('c_region', int(rrow['c_region_start']), int(rrow['c_region_end'])))
    if rrow['is_single_domain']:
        mid = (rrow['n_domain_start'] + rrow['n_domain_end']) // 2
        regions.append(('n_half', int(rrow['n_domain_start']), int(mid)))
        regions.append(('c_half', int(mid + 1), int(rrow['n_domain_end'])))

    for region_name, rstart, rend in regions:
        # SS fractions from DSSP
        ss_codes = [resnum_to_ss.get(rn, '-') for rn in range(rstart, rend + 1)
                     if rn in resnum_to_ss]
        ss_fracs = compute_ss_fractions(ss_codes)

        # pLDDT from CIF
        plddt_vals = [ca_atoms[rn]['bfactor'] for rn in range(rstart, rend + 1)
                      if rn in ca_atoms]
        plddt_metrics = compute_plddt_metrics(plddt_vals)

        struct_rows.append({
            'uniprot_accession': acc,
            'region': region_name,
            **ss_fracs,
            **plddt_metrics,
        })

    if (cif_parsed % 200) == 0:
        print(f"  Processed {cif_parsed} CIF files...")

struct_df = pd.DataFrame(struct_rows)
struct_df.to_csv(OUT_DIR / "structure_metrics.tsv", sep='\t', index=False)
print(f"Structure metrics saved: {len(struct_df)} region entries")
print(f"  CIF files parsed: {cif_parsed}")
if cif_missing > 0:
    print(f"  CIF files missing: {cif_missing}")
print()

# ═══════════════════════════════════════════════════════════════════════════
# STEP F4: Contact Order (Plaxco-Simons) per region
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("STEP F4: Relative contact order per region")
print("=" * 72)

co_rows = []
co_processed = 0
co_errors = 0

for idx, rrow in region_df.iterrows():
    acc = rrow['uniprot_accession']

    # Parse CIF for CA coordinates
    cif_path = CIF_DIR / f"AF-{acc}-F1-model_v6.cif"
    if not cif_path.exists():
        cif_candidates = list(CIF_DIR.glob(f"AF-{acc}-F1-model_v*.cif"))
        if cif_candidates:
            cif_path = cif_candidates[0]
        else:
            continue

    try:
        ca_atoms = parse_cif_ca_atoms(cif_path)
    except Exception as e:
        co_errors += 1
        continue

    # Define regions for contact order
    regions = []
    regions.append(('n_domain', int(rrow['n_domain_start']), int(rrow['n_domain_end'])))
    if rrow['c_region_length'] > 0:
        regions.append(('c_region', int(rrow['c_region_start']), int(rrow['c_region_end'])))
    if rrow['pre_domain_length'] > 0:
        regions.append(('pre_domain', int(rrow['pre_domain_start']), int(rrow['pre_domain_end'])))
    # Full protein
    regions.append(('full_protein', 1, int(rrow['protein_length'])))
    # For single-domain: n_half and c_half
    if rrow['is_single_domain']:
        mid = (rrow['n_domain_start'] + rrow['n_domain_end']) // 2
        regions.append(('n_half', int(rrow['n_domain_start']), int(mid)))
        regions.append(('c_half', int(mid + 1), int(rrow['n_domain_end'])))

    for region_name, rstart, rend in regions:
        # Get CA coords in this range
        resnums_in_range = sorted([rn for rn in ca_atoms if rstart <= rn <= rend])
        if len(resnums_in_range) < 2:
            co_rows.append({
                'uniprot_accession': acc,
                'region': region_name,
                'n_residues': len(resnums_in_range),
                'n_contacts': 0,
                'abs_contact_order': np.nan,
                'relative_contact_order': np.nan,
            })
            continue

        coords = np.array([[ca_atoms[rn]['x'], ca_atoms[rn]['y'], ca_atoms[rn]['z']]
                           for rn in resnums_in_range])
        resnums_arr = np.array(resnums_in_range)

        n_res, n_contacts, aco, rco = compute_contact_order(coords, resnums_arr)

        co_rows.append({
            'uniprot_accession': acc,
            'region': region_name,
            'n_residues': n_res,
            'n_contacts': n_contacts,
            'abs_contact_order': aco,
            'relative_contact_order': rco,
        })

    co_processed += 1
    if (co_processed % 100) == 0:
        print(f"  Contact order computed for {co_processed} proteins...")

co_df = pd.DataFrame(co_rows)
co_df.to_csv(OUT_DIR / "contact_order.tsv", sep='\t', index=False)
print(f"Contact order saved: {len(co_df)} region entries ({co_processed} proteins)")
if co_errors > 0:
    print(f"  Errors: {co_errors}")
print()

# ═══════════════════════════════════════════════════════════════════════════
# STEP F5: Within-protein paired comparison (N-domain vs C-region)
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("STEP F5: Within-protein paired comparison")
print("=" * 72)

# Load GroEL class data
groel_df = pd.read_csv(GROEL_FILE, sep='\t')
groel_lookup = dict(zip(groel_df['current_accession'], groel_df['groel_class']))

# Merge all metrics into a single frame, pivoted by region
# Sequence metrics
seq_pivot = seq_df.pivot(index='uniprot_accession', columns='region')
seq_pivot.columns = [f"{col[0]}_{col[1]}" for col in seq_pivot.columns]
seq_pivot = seq_pivot.reset_index()

# Structure metrics
struct_pivot = struct_df.pivot(index='uniprot_accession', columns='region')
struct_pivot.columns = [f"{col[0]}_{col[1]}" for col in struct_pivot.columns]
struct_pivot = struct_pivot.reset_index()

# Contact order metrics
co_pivot = co_df.pivot(index='uniprot_accession', columns='region')
co_pivot.columns = [f"{col[0]}_{col[1]}" for col in co_pivot.columns]
co_pivot = co_pivot.reset_index()

# Merge all
merged = region_df[['uniprot_accession', 'source_dataset', 'n_domains', 'is_single_domain',
                     'protein_length', 'pre_domain_length', 'n_domain_length', 'c_region_length']].copy()
merged = merged.merge(seq_pivot, on='uniprot_accession', how='left')
merged = merged.merge(struct_pivot, on='uniprot_accession', how='left')
merged = merged.merge(co_pivot, on='uniprot_accession', how='left')

# Add GroEL class
merged['groel_class'] = merged['uniprot_accession'].map(groel_lookup).fillna('')

# ── Multi-domain paired comparison ──
multi = merged[~merged['is_single_domain']].copy()
# Only keep those that have valid C-region (length > 0)
multi = multi[multi['c_region_length'] > 0].copy()

print(f"Multi-domain proteins with C-region: {len(multi)}")

# Define the metrics for paired comparison
paired_metrics = [
    ('length', 'length_n_domain', 'length_c_region'),
    ('net_charge', 'net_charge_n_domain', 'net_charge_c_region'),
    ('frac_hydrophobic', 'frac_hydrophobic_n_domain', 'frac_hydrophobic_c_region'),
    ('frac_helix', 'frac_helix_n_domain', 'frac_helix_c_region'),
    ('frac_strand', 'frac_strand_n_domain', 'frac_strand_c_region'),
    ('mean_plddt', 'mean_plddt_n_domain', 'mean_plddt_c_region'),
    ('frac_plddt_gt70', 'frac_plddt_gt70_n_domain', 'frac_plddt_gt70_c_region'),
    ('relative_contact_order', 'relative_contact_order_n_domain', 'relative_contact_order_c_region'),
]

# Build paired comparison table
paired_rows = []
for _, row in multi.iterrows():
    prow = {
        'uniprot_accession': row['uniprot_accession'],
        'source_dataset': row['source_dataset'],
        'n_domains': row['n_domains'],
        'groel_class': row['groel_class'],
    }
    for metric_name, n_col, c_col in paired_metrics:
        n_val = row.get(n_col, np.nan)
        c_val = row.get(c_col, np.nan)
        prow[f'{metric_name}_n_domain'] = n_val
        prow[f'{metric_name}_c_region'] = c_val
        if pd.notna(n_val) and pd.notna(c_val):
            prow[f'{metric_name}_diff'] = n_val - c_val
            prow[f'{metric_name}_ratio'] = n_val / c_val if c_val != 0 else np.nan
        else:
            prow[f'{metric_name}_diff'] = np.nan
            prow[f'{metric_name}_ratio'] = np.nan
    paired_rows.append(prow)

paired_df = pd.DataFrame(paired_rows)
paired_df.to_csv(OUT_DIR / "n_vs_c_paired.tsv", sep='\t', index=False)
print(f"Paired comparison table saved: {len(paired_df)} proteins")

# ── Summary statistics ──
print("\n" + "=" * 72)
print("SUMMARY STATISTICS")
print("=" * 72)

# Counts per dataset
print("\n--- Protein counts per dataset ---")
for ds in sorted(region_df['source_dataset'].unique()):
    sub = region_df[region_df['source_dataset'] == ds]
    n_single = sub['is_single_domain'].sum()
    n_multi = (~sub['is_single_domain']).sum()
    print(f"  {ds:30s}: {len(sub):4d} total  |  single={n_single:4d}  multi={n_multi:4d}")

# Overall
n_single_all = region_df['is_single_domain'].sum()
n_multi_all = (~region_df['is_single_domain']).sum()
print(f"  {'TOTAL':30s}: {len(region_df):4d} total  |  single={n_single_all:4d}  multi={n_multi_all:4d}")

# ── Paired comparison stats for multi-domain ──
print("\n--- Multi-domain paired comparison: N-domain vs C-region ---")
print(f"{'Metric':<25s} {'Mean(N)':>10s} {'Mean(C)':>10s} {'Mean(diff)':>11s} {'N>C':>6s} {'N<C':>6s} {'Wilcoxon p':>12s}")
print("-" * 82)

for metric_name, n_col, c_col in paired_metrics:
    diff_col = f'{metric_name}_diff'
    if diff_col not in paired_df.columns:
        continue

    diffs = paired_df[diff_col].dropna()
    if len(diffs) < 5:
        continue

    n_vals = paired_df[f'{metric_name}_n_domain'].dropna()
    c_vals = paired_df[f'{metric_name}_c_region'].dropna()

    mean_n = n_vals.mean()
    mean_c = c_vals.mean()
    mean_diff = diffs.mean()
    n_pos = (diffs > 0).sum()
    n_neg = (diffs < 0).sum()

    try:
        stat, pval = stats.wilcoxon(diffs, alternative='two-sided')
    except Exception:
        pval = np.nan

    print(f"{metric_name:<25s} {mean_n:10.4f} {mean_c:10.4f} {mean_diff:+11.4f} {n_pos:6d} {n_neg:6d} {pval:12.2e}")

# ── Breakdown by dataset ──
# Determine primary datasets
dataset_names = set()
for ds in region_df['source_dataset'].unique():
    for d in str(ds).split(','):
        d = d.strip()
        if d:
            dataset_names.add(d)

print("\n--- Breakdown by primary dataset (multi-domain only) ---")
for ds_name in sorted(dataset_names):
    sub = paired_df[paired_df['source_dataset'].str.contains(ds_name, na=False)]
    if len(sub) < 5:
        print(f"\n  [{ds_name}]: n={len(sub)}, too few for statistics")
        continue

    print(f"\n  [{ds_name}] n={len(sub)}")
    print(f"  {'Metric':<25s} {'Mean(diff)':>11s} {'N>C':>6s} {'N<C':>6s} {'Wilcoxon p':>12s}")
    print(f"  {'-'*60}")

    for metric_name, n_col, c_col in paired_metrics:
        diff_col = f'{metric_name}_diff'
        diffs = sub[diff_col].dropna()
        if len(diffs) < 5:
            continue
        n_pos = (diffs > 0).sum()
        n_neg = (diffs < 0).sum()
        try:
            stat, pval = stats.wilcoxon(diffs, alternative='two-sided')
        except Exception:
            pval = np.nan
        print(f"  {metric_name:<25s} {diffs.mean():+11.4f} {n_pos:6d} {n_neg:6d} {pval:12.2e}")

# ── Single-domain: N-half vs C-half ──
print("\n--- Single-domain proteins: N-half vs C-half of the domain ---")
single = merged[merged['is_single_domain']].copy()
print(f"Single-domain proteins: {len(single)}")

single_paired_metrics = [
    ('frac_hydrophobic', 'frac_hydrophobic_n_half', 'frac_hydrophobic_c_half'),
    ('frac_helix', 'frac_helix_n_half', 'frac_helix_c_half'),
    ('frac_strand', 'frac_strand_n_half', 'frac_strand_c_half'),
    ('mean_plddt', 'mean_plddt_n_half', 'mean_plddt_c_half'),
    ('relative_contact_order', 'relative_contact_order_n_half', 'relative_contact_order_c_half'),
]

print(f"{'Metric':<25s} {'Mean(N-half)':>12s} {'Mean(C-half)':>12s} {'Mean(diff)':>11s} {'N>C':>6s} {'N<C':>6s} {'Wilcoxon p':>12s}")
print("-" * 84)

for metric_name, n_col, c_col in single_paired_metrics:
    if n_col not in single.columns or c_col not in single.columns:
        continue
    n_vals = single[n_col]
    c_vals = single[c_col]
    valid = n_vals.notna() & c_vals.notna()
    if valid.sum() < 5:
        continue
    diffs = (n_vals[valid] - c_vals[valid])
    mean_n = n_vals[valid].mean()
    mean_c = c_vals[valid].mean()
    n_pos = (diffs > 0).sum()
    n_neg = (diffs < 0).sum()
    try:
        stat, pval = stats.wilcoxon(diffs, alternative='two-sided')
    except Exception:
        pval = np.nan
    print(f"{metric_name:<25s} {mean_n:12.4f} {mean_c:12.4f} {diffs.mean():+11.4f} {n_pos:6d} {n_neg:6d} {pval:12.2e}")

# ── Final summary ──
print("\n" + "=" * 72)
print("MODULE F COMPLETE")
print("=" * 72)
print(f"\nOutput files:")
print(f"  1. {OUT_DIR / 'region_boundaries.tsv'}")
print(f"  2. {OUT_DIR / 'sequence_metrics.tsv'}")
print(f"  3. {OUT_DIR / 'structure_metrics.tsv'}")
print(f"  4. {OUT_DIR / 'contact_order.tsv'}")
print(f"  5. {OUT_DIR / 'n_vs_c_paired.tsv'}")
print(f"\nTotal proteins analyzed: {len(region_df)}")
print(f"  Single-domain: {n_single_all}")
print(f"  Multi-domain:  {n_multi_all}")
