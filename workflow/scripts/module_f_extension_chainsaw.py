#!/usr/bin/env python3
"""
Module F Extension: Include Chainsaw-assigned proteins in N-vs-C analysis.

The original Module F analyzed 1,151 CATH-assigned proteins. This extension
processes the 236 Chainsaw-assigned proteins that previously lacked CATH domain
assignments and merges the results with the original Module F outputs.

Steps:
  1. Read Chainsaw domain predictions, define region boundaries for new proteins
  2. Compute sequence metrics (charge, hydrophobicity, aromaticity) per region
  3. Compute structure metrics (SS fractions, pLDDT) using DSSP data
  4. Compute contact order (Plaxco-Simons relative CO) from AlphaFold CIF
  5. Run paired within-protein Wilcoxon comparison for multi-domain proteins
  6. Merge with existing results and produce updated statistics
  7. Save extended output files
  8. Print comparison report: original vs extended
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
CHAINSAW_FILE  = BASE / "results/domains/chainsaw_domain_predictions.tsv"
DSSP_FILE      = BASE / "results/structures/dssp_per_residue.tsv"
STRUCT_INDEX   = BASE / "results/structures/structure_index.tsv"
GROEL_FILE     = BASE / "data/processed/groel_substrates_standardized.tsv"
CIF_DIR        = BASE / "data/raw/alphafold/pilot"
OUT_DIR        = BASE / "results/termini"

# Original Module F outputs
ORIG_REGION    = OUT_DIR / "region_boundaries.tsv"
ORIG_PAIRED    = OUT_DIR / "n_vs_c_paired.tsv"
ORIG_CO        = OUT_DIR / "contact_order.tsv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Amino acid definitions ────────────────────────────────────────────────
AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

HYDROPHOBIC = set('AVILFWM')
CHARGED     = set('KRDE')
POLAR       = set('STNQYHC')
AROMATIC    = set('FWY')
SMALL       = set('GAS')
POSITIVE    = set('KR')
NEGATIVE    = set('DE')


# ── Helper functions (same as original Module F) ──────────────────────────

def parse_cif_ca_atoms(cif_path):
    """Parse CA atom records from an AlphaFold CIF file."""
    ca_atoms = {}
    with open(cif_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
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
    pos = sum(1 for r in residues_1letter if r in POSITIVE)
    neg = sum(1 for r in residues_1letter if r in NEGATIVE)
    return pos - neg


def compute_fractions(residues_1letter):
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
    n = len(ss_codes)
    if n == 0:
        return {'frac_helix': np.nan, 'frac_strand': np.nan, 'frac_coil': np.nan}
    helix_codes = {'H', 'G', 'I'}
    strand_codes = {'E', 'B'}
    n_helix  = sum(1 for s in ss_codes if s in helix_codes)
    n_strand = sum(1 for s in ss_codes if s in strand_codes)
    n_coil   = n - n_helix - n_strand
    return {
        'frac_helix':  n_helix / n,
        'frac_strand': n_strand / n,
        'frac_coil':   n_coil / n,
    }


def compute_plddt_metrics(plddt_values):
    if len(plddt_values) == 0:
        return {'mean_plddt': np.nan, 'frac_plddt_gt70': np.nan, 'frac_plddt_gt90': np.nan}
    arr = np.array(plddt_values)
    return {
        'mean_plddt':      float(np.mean(arr)),
        'frac_plddt_gt70': float(np.mean(arr > 70)),
        'frac_plddt_gt90': float(np.mean(arr > 90)),
    }


def compute_contact_order(coords, residue_numbers, min_seq_sep=6, distance_cutoff=8.0):
    n = len(coords)
    if n < 2:
        return n, 0, np.nan, np.nan

    coords = np.array(coords, dtype=np.float64)
    resnums = np.array(residue_numbers, dtype=np.int64)

    if n > 3000:
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
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dists = np.sqrt(np.sum(diff**2, axis=2))
        seq_sep = np.abs(resnums[:, np.newaxis] - resnums[np.newaxis, :])
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
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("MODULE F EXTENSION: Adding Chainsaw-assigned proteins")
print("=" * 72)

# Load Chainsaw predictions
chainsaw_df = pd.read_csv(CHAINSAW_FILE, sep='\t')
print(f"Chainsaw domain predictions: {len(chainsaw_df)} rows")

# Load original region boundaries to determine which proteins are NEW
orig_region_df = pd.read_csv(ORIG_REGION, sep='\t')
existing_accs = set(orig_region_df['uniprot_accession'])
print(f"Original Module F proteins: {len(existing_accs)}")

# Identify Chainsaw-only proteins (not already in CATH analysis)
chainsaw_accs_all = set(chainsaw_df['uniprot_accession'])
new_accs = chainsaw_accs_all - existing_accs
print(f"Chainsaw proteins total: {len(chainsaw_accs_all)}")
print(f"New proteins (Chainsaw-only): {len(new_accs)}")

# Filter chainsaw_df to new proteins only
chainsaw_new = chainsaw_df[chainsaw_df['uniprot_accession'].isin(new_accs)].copy()

# Load structure index for protein lengths and source_dataset
struct_df = pd.read_csv(STRUCT_INDEX, sep='\t')
length_lookup = dict(zip(struct_df['uniprot_accession'], struct_df['residues_modeled']))
source_lookup = dict(zip(struct_df['uniprot_accession'], struct_df['source_dataset']))

# Load DSSP per-residue
print("Loading DSSP per-residue data...")
dssp_df = pd.read_csv(DSSP_FILE, sep='\t')
dssp_grouped = {acc: grp for acc, grp in dssp_df.groupby('uniprot_accession')}
print(f"  DSSP data for {len(dssp_grouped)} proteins")

# Load GroEL
groel_df = pd.read_csv(GROEL_FILE, sep='\t')
groel_lookup = dict(zip(groel_df['current_accession'], groel_df['groel_class']))

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: Define three regions for Chainsaw proteins
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 1: Define region boundaries for Chainsaw proteins")
print("=" * 72)

region_rows = []
skipped_no_length = 0

for acc in sorted(new_accs):
    prot_domains = chainsaw_new[chainsaw_new['uniprot_accession'] == acc].copy()
    prot_domains = prot_domains.sort_values('domain_start')

    if len(prot_domains) == 0:
        continue

    # Get protein length from structure_index
    prot_len = length_lookup.get(acc)
    if prot_len is None:
        # Try to infer from max domain_end
        prot_len = int(prot_domains['domain_end'].max())
        skipped_no_length += 1

    n_dom = int(prot_domains.iloc[0]['n_domains_total'])

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

new_region_df = pd.DataFrame(region_rows)
n_new_single = new_region_df['is_single_domain'].sum()
n_new_multi = (~new_region_df['is_single_domain']).sum()

print(f"New region boundaries: {len(new_region_df)} proteins")
print(f"  Single-domain: {n_new_single}")
print(f"  Multi-domain:  {n_new_multi}")
print(f"  Mean pre-domain length: {new_region_df['pre_domain_length'].mean():.1f}")
print(f"  Mean N-domain length:   {new_region_df['n_domain_length'].mean():.1f}")
print(f"  Mean C-region length:   {new_region_df['c_region_length'].mean():.1f}")
if skipped_no_length > 0:
    print(f"  Proteins with inferred length: {skipped_no_length}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: Sequence-derived metrics per region
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 2: Sequence-derived metrics for Chainsaw proteins")
print("=" * 72)

seq_rows = []
seq_skipped = 0

for _, rrow in new_region_df.iterrows():
    acc = rrow['uniprot_accession']
    if acc not in dssp_grouped:
        seq_skipped += 1
        continue

    pdssp = dssp_grouped[acc]
    resnum_to_name = dict(zip(pdssp['residue_number'], pdssp['residue_name']))

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
        residues_1letter = []
        for rn in range(rstart, rend + 1):
            res_name = resnum_to_name.get(rn, '')
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
print(f"Sequence metrics computed: {len(seq_df)} region entries")
if seq_skipped > 0:
    print(f"  Skipped {seq_skipped} proteins (no DSSP data)")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: Structure-derived metrics (DSSP + pLDDT)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 3: Structure-derived metrics for Chainsaw proteins")
print("=" * 72)

struct_rows = []
cif_missing = 0
cif_parsed = 0

for idx, rrow in new_region_df.iterrows():
    acc = rrow['uniprot_accession']

    pdssp = dssp_grouped.get(acc)
    if pdssp is None:
        continue

    resnum_to_ss = dict(zip(pdssp['residue_number'], pdssp['ss_code']))

    cif_path = CIF_DIR / f"AF-{acc}-F1-model_v6.cif"
    if not cif_path.exists():
        cif_candidates = list(CIF_DIR.glob(f"AF-{acc}-F1-model_v*.cif"))
        if cif_candidates:
            cif_path = cif_candidates[0]
        else:
            cif_missing += 1
            continue

    ca_atoms = parse_cif_ca_atoms(cif_path)
    cif_parsed += 1

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
        ss_codes = [resnum_to_ss.get(rn, '-') for rn in range(rstart, rend + 1)
                     if rn in resnum_to_ss]
        ss_fracs = compute_ss_fractions(ss_codes)

        plddt_vals = [ca_atoms[rn]['bfactor'] for rn in range(rstart, rend + 1)
                      if rn in ca_atoms]
        plddt_metrics = compute_plddt_metrics(plddt_vals)

        struct_rows.append({
            'uniprot_accession': acc,
            'region': region_name,
            **ss_fracs,
            **plddt_metrics,
        })

    if (cif_parsed % 50) == 0:
        print(f"  Processed {cif_parsed} CIF files...")

struct_df_new = pd.DataFrame(struct_rows)
print(f"Structure metrics computed: {len(struct_df_new)} region entries")
print(f"  CIF files parsed: {cif_parsed}")
if cif_missing > 0:
    print(f"  CIF files missing: {cif_missing}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: Contact order per region
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 4: Contact order for Chainsaw proteins")
print("=" * 72)

co_rows = []
co_processed = 0
co_errors = 0

for idx, rrow in new_region_df.iterrows():
    acc = rrow['uniprot_accession']

    cif_path = CIF_DIR / f"AF-{acc}-F1-model_v6.cif"
    if not cif_path.exists():
        cif_candidates = list(CIF_DIR.glob(f"AF-{acc}-F1-model_v*.cif"))
        if cif_candidates:
            cif_path = cif_candidates[0]
        else:
            continue

    try:
        ca_atoms = parse_cif_ca_atoms(cif_path)
    except Exception:
        co_errors += 1
        continue

    regions = []
    regions.append(('n_domain', int(rrow['n_domain_start']), int(rrow['n_domain_end'])))
    if rrow['c_region_length'] > 0:
        regions.append(('c_region', int(rrow['c_region_start']), int(rrow['c_region_end'])))
    if rrow['pre_domain_length'] > 0:
        regions.append(('pre_domain', int(rrow['pre_domain_start']), int(rrow['pre_domain_end'])))
    regions.append(('full_protein', 1, int(rrow['protein_length'])))
    if rrow['is_single_domain']:
        mid = (rrow['n_domain_start'] + rrow['n_domain_end']) // 2
        regions.append(('n_half', int(rrow['n_domain_start']), int(mid)))
        regions.append(('c_half', int(mid + 1), int(rrow['n_domain_end'])))

    for region_name, rstart, rend in regions:
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
    if (co_processed % 50) == 0:
        print(f"  Contact order computed for {co_processed} proteins...")

co_df_new = pd.DataFrame(co_rows)
print(f"Contact order computed: {len(co_df_new)} region entries ({co_processed} proteins)")
if co_errors > 0:
    print(f"  Errors: {co_errors}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Build paired comparison for new multi-domain proteins
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 5: Paired comparison for new Chainsaw multi-domain proteins")
print("=" * 72)

# Pivot new metrics
seq_pivot = seq_df.pivot(index='uniprot_accession', columns='region')
seq_pivot.columns = [f"{col[0]}_{col[1]}" for col in seq_pivot.columns]
seq_pivot = seq_pivot.reset_index()

struct_pivot = struct_df_new.pivot(index='uniprot_accession', columns='region')
struct_pivot.columns = [f"{col[0]}_{col[1]}" for col in struct_pivot.columns]
struct_pivot = struct_pivot.reset_index()

co_pivot = co_df_new.pivot(index='uniprot_accession', columns='region')
co_pivot.columns = [f"{col[0]}_{col[1]}" for col in co_pivot.columns]
co_pivot = co_pivot.reset_index()

# Merge all new metrics
merged_new = new_region_df[['uniprot_accession', 'source_dataset', 'n_domains',
                             'is_single_domain', 'protein_length',
                             'pre_domain_length', 'n_domain_length', 'c_region_length']].copy()
merged_new = merged_new.merge(seq_pivot, on='uniprot_accession', how='left')
merged_new = merged_new.merge(struct_pivot, on='uniprot_accession', how='left')
merged_new = merged_new.merge(co_pivot, on='uniprot_accession', how='left')
merged_new['groel_class'] = merged_new['uniprot_accession'].map(groel_lookup).fillna('')

# Multi-domain with valid C-region
multi_new = merged_new[~merged_new['is_single_domain']].copy()
multi_new = multi_new[multi_new['c_region_length'] > 0].copy()
print(f"New multi-domain proteins with C-region: {len(multi_new)}")

# Define paired metrics (same as original)
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

# Build paired rows for new proteins
new_paired_rows = []
for _, row in multi_new.iterrows():
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
    new_paired_rows.append(prow)

new_paired_df = pd.DataFrame(new_paired_rows)
print(f"New paired comparison rows: {len(new_paired_df)}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6: Merge with original results
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 6: Merge with original Module F results")
print("=" * 72)

# Merge region boundaries
ext_region_df = pd.concat([orig_region_df, new_region_df], ignore_index=True)
ext_region_df = ext_region_df.sort_values('uniprot_accession').reset_index(drop=True)
ext_region_df.to_csv(OUT_DIR / "region_boundaries_extended.tsv", sep='\t', index=False)
print(f"Extended region boundaries: {len(ext_region_df)} proteins (was {len(orig_region_df)})")

# Merge paired comparison
orig_paired_df = pd.read_csv(ORIG_PAIRED, sep='\t')
ext_paired_df = pd.concat([orig_paired_df, new_paired_df], ignore_index=True)
ext_paired_df = ext_paired_df.sort_values('uniprot_accession').reset_index(drop=True)
ext_paired_df.to_csv(OUT_DIR / "n_vs_c_paired_extended.tsv", sep='\t', index=False)
print(f"Extended paired comparison: {len(ext_paired_df)} proteins (was {len(orig_paired_df)})")

# Merge contact order
orig_co_df = pd.read_csv(ORIG_CO, sep='\t')
ext_co_df = pd.concat([orig_co_df, co_df_new], ignore_index=True)
ext_co_df = ext_co_df.sort_values(['uniprot_accession', 'region']).reset_index(drop=True)
ext_co_df.to_csv(OUT_DIR / "contact_order_extended.tsv", sep='\t', index=False)
print(f"Extended contact order: {len(ext_co_df)} entries (was {len(orig_co_df)})")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7: Extended statistics and Wilcoxon tests
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("STEP 7: Extended paired comparison statistics")
print("=" * 72)

# --- Stats for NEW Chainsaw-only multi-domain ---
print("\n--- NEW Chainsaw-only multi-domain proteins (N-domain vs C-region) ---")
if len(new_paired_df) >= 5:
    print(f"n = {len(new_paired_df)}")
    print(f"{'Metric':<25s} {'Mean(N)':>10s} {'Mean(C)':>10s} {'Mean(diff)':>11s} {'N>C':>6s} {'N<C':>6s} {'Wilcoxon p':>12s}")
    print("-" * 82)

    for metric_name, n_col, c_col in paired_metrics:
        diff_col = f'{metric_name}_diff'
        if diff_col not in new_paired_df.columns:
            continue
        diffs = new_paired_df[diff_col].dropna()
        if len(diffs) < 5:
            continue
        n_vals = new_paired_df[f'{metric_name}_n_domain'].dropna()
        c_vals = new_paired_df[f'{metric_name}_c_region'].dropna()
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
else:
    print(f"  Only {len(new_paired_df)} new multi-domain proteins -- too few for statistics")

# --- Stats for EXTENDED (all) multi-domain ---
print("\n--- EXTENDED dataset: all multi-domain proteins (N-domain vs C-region) ---")
print(f"n = {len(ext_paired_df)}")
print(f"{'Metric':<25s} {'Mean(N)':>10s} {'Mean(C)':>10s} {'Mean(diff)':>11s} {'N>C':>6s} {'N<C':>6s} {'Wilcoxon p':>12s}")
print("-" * 82)

for metric_name, n_col, c_col in paired_metrics:
    diff_col = f'{metric_name}_diff'
    if diff_col not in ext_paired_df.columns:
        continue
    diffs = ext_paired_df[diff_col].dropna()
    if len(diffs) < 5:
        continue
    n_vals = ext_paired_df[f'{metric_name}_n_domain'].dropna()
    c_vals = ext_paired_df[f'{metric_name}_c_region'].dropna()
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


# ═══════════════════════════════════════════════════════════════════════════
# STEP 8: Comparison report -- Original vs Extended
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("COMPARISON REPORT: Original vs Extended Module F")
print("=" * 72)

# Counts
orig_total = len(orig_region_df)
orig_single = orig_region_df['is_single_domain'].sum()
orig_multi = (~orig_region_df['is_single_domain']).sum()
orig_paired_n = len(orig_paired_df)

ext_total = len(ext_region_df)
ext_single = ext_region_df['is_single_domain'].sum()
ext_multi = (~ext_region_df['is_single_domain']).sum()
ext_paired_n = len(ext_paired_df)

print(f"\n{'':30s} {'Original':>12s} {'Extended':>12s} {'Added':>12s}")
print("-" * 70)
print(f"{'Total proteins':30s} {orig_total:12d} {ext_total:12d} {ext_total - orig_total:+12d}")
print(f"{'Single-domain':30s} {orig_single:12d} {ext_single:12d} {ext_single - orig_single:+12d}")
print(f"{'Multi-domain':30s} {orig_multi:12d} {ext_multi:12d} {ext_multi - orig_multi:+12d}")
print(f"{'Multi-domain paired (tested)':30s} {orig_paired_n:12d} {ext_paired_n:12d} {ext_paired_n - orig_paired_n:+12d}")

# Coverage
total_project_proteins = 1390
print(f"\n{'Domain coverage':30s} {orig_total}/{total_project_proteins} ({100*orig_total/total_project_proteins:.1f}%) -> {ext_total}/{total_project_proteins} ({100*ext_total/total_project_proteins:.1f}%)")

# Compare key Wilcoxon results side-by-side
print("\n--- Key Wilcoxon p-values: Original (n={}) vs Extended (n={}) ---".format(
    orig_paired_n, ext_paired_n))
print(f"{'Metric':<25s} {'p (original)':>14s} {'p (extended)':>14s} {'Direction':>10s}")
print("-" * 68)

for metric_name, n_col, c_col in paired_metrics:
    diff_col = f'{metric_name}_diff'

    # Original
    if diff_col in orig_paired_df.columns:
        diffs_orig = orig_paired_df[diff_col].dropna()
        if len(diffs_orig) >= 5:
            try:
                _, p_orig = stats.wilcoxon(diffs_orig, alternative='two-sided')
            except Exception:
                p_orig = np.nan
            mean_diff_orig = diffs_orig.mean()
        else:
            p_orig = np.nan
            mean_diff_orig = np.nan
    else:
        p_orig = np.nan
        mean_diff_orig = np.nan

    # Extended
    if diff_col in ext_paired_df.columns:
        diffs_ext = ext_paired_df[diff_col].dropna()
        if len(diffs_ext) >= 5:
            try:
                _, p_ext = stats.wilcoxon(diffs_ext, alternative='two-sided')
            except Exception:
                p_ext = np.nan
            mean_diff_ext = diffs_ext.mean()
        else:
            p_ext = np.nan
            mean_diff_ext = np.nan
    else:
        p_ext = np.nan
        mean_diff_ext = np.nan

    direction = "N>C" if mean_diff_ext > 0 else "N<C" if mean_diff_ext < 0 else "="
    p_orig_str = f"{p_orig:.2e}" if pd.notna(p_orig) else "N/A"
    p_ext_str = f"{p_ext:.2e}" if pd.notna(p_ext) else "N/A"
    print(f"{metric_name:<25s} {p_orig_str:>14s} {p_ext_str:>14s} {direction:>10s}")

# Dataset breakdown for extended
print("\n--- Extended dataset: protein counts per source ---")
for ds in sorted(ext_region_df['source_dataset'].dropna().unique()):
    sub = ext_region_df[ext_region_df['source_dataset'] == ds]
    n_s = sub['is_single_domain'].sum()
    n_m = (~sub['is_single_domain']).sum()
    print(f"  {ds:30s}: {len(sub):4d} total  |  single={n_s:4d}  multi={n_m:4d}")

total_s = ext_region_df['is_single_domain'].sum()
total_m = (~ext_region_df['is_single_domain']).sum()
print(f"  {'TOTAL':30s}: {len(ext_region_df):4d} total  |  single={total_s:4d}  multi={total_m:4d}")

# ── Final ──
print("\n" + "=" * 72)
print("MODULE F EXTENSION COMPLETE")
print("=" * 72)
print(f"\nOutput files:")
print(f"  1. {OUT_DIR / 'region_boundaries_extended.tsv'}")
print(f"  2. {OUT_DIR / 'n_vs_c_paired_extended.tsv'}")
print(f"  3. {OUT_DIR / 'contact_order_extended.tsv'}")
print(f"\nProteins analyzed:")
print(f"  Original CATH:  {orig_total}")
print(f"  New Chainsaw:   {len(new_region_df)}")
print(f"  Extended total: {ext_total}")
print(f"  Multi-domain paired tests: {orig_paired_n} -> {ext_paired_n}")
