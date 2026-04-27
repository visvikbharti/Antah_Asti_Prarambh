#!/bin/bash
# =============================================================================
# Antah Asti Prarambh — Download External Data
# =============================================================================
# Downloads all external datasets required for the analysis.
# Run from the project root directory.
#
# Usage:
#   bash scripts/download_external_data.sh
#
# Requirements: curl or wget
# =============================================================================

set -euo pipefail

# Resolve project root (parent of scripts/)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "============================================================"
echo "  Antah Asti Prarambh — External Data Download"
echo "  Project: ${PROJECT_DIR}"
echo "  Date: $(date)"
echo "============================================================"

# Create directories
mkdir -p "${PROJECT_DIR}/data/raw/uniprot"
mkdir -p "${PROJECT_DIR}/data/raw/mitocarta"

# -----------------------------------------------------------------------------
# 1. E. coli K-12 proteome (UniProt UP000000625)
# -----------------------------------------------------------------------------
ECOLI_FASTA="${PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.fasta"
ECOLI_TSV="${PROJECT_DIR}/data/raw/uniprot/ecoli_k12_proteome.tsv"

if [ -f "$ECOLI_FASTA" ]; then
    echo "[SKIP] E. coli FASTA already exists"
else
    echo "[DOWNLOAD] E. coli K-12 proteome FASTA..."
    curl -sL "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:UP000000625+AND+reviewed:true" \
        -o "$ECOLI_FASTA"
    echo "  -> $(grep -c '^>' "$ECOLI_FASTA") sequences"
fi

if [ -f "$ECOLI_TSV" ]; then
    echo "[SKIP] E. coli TSV already exists"
else
    echo "[DOWNLOAD] E. coli K-12 proteome metadata TSV..."
    curl -sL "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=proteome:UP000000625+AND+reviewed:true&fields=accession,id,protein_name,gene_names,organism_name,length,mass,go_id,go_p,go_c,go_f,ft_transit,ft_signal,cc_subcellular_location" \
        -o "$ECOLI_TSV"
    echo "  -> $(wc -l < "$ECOLI_TSV") rows"
fi

# -----------------------------------------------------------------------------
# 2. Human reference proteome (UniProt UP000005640)
# -----------------------------------------------------------------------------
HUMAN_FASTA="${PROJECT_DIR}/data/raw/uniprot/human_proteome.fasta"
HUMAN_TSV="${PROJECT_DIR}/data/raw/uniprot/human_proteome.tsv"

if [ -f "$HUMAN_FASTA" ]; then
    echo "[SKIP] Human FASTA already exists"
else
    echo "[DOWNLOAD] Human reference proteome FASTA..."
    curl -sL "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:UP000005640+AND+reviewed:true" \
        -o "$HUMAN_FASTA"
    echo "  -> $(grep -c '^>' "$HUMAN_FASTA") sequences"
fi

if [ -f "$HUMAN_TSV" ]; then
    echo "[SKIP] Human TSV already exists"
else
    echo "[DOWNLOAD] Human reference proteome metadata TSV..."
    curl -sL "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=proteome:UP000005640+AND+reviewed:true&fields=accession,id,protein_name,gene_names,organism_name,length,mass,go_id,go_p,go_c,go_f,ft_transit,ft_signal,cc_subcellular_location" \
        -o "$HUMAN_TSV"
    echo "  -> $(wc -l < "$HUMAN_TSV") rows"
fi

# -----------------------------------------------------------------------------
# 3. MitoCarta 3.0 (Broad Institute)
# -----------------------------------------------------------------------------
MITOCARTA="${PROJECT_DIR}/data/raw/mitocarta/Human.MitoCarta3.0.xls"

if [ -f "$MITOCARTA" ]; then
    echo "[SKIP] MitoCarta 3.0 already exists"
else
    echo "[DOWNLOAD] MitoCarta 3.0..."
    curl -sL "https://personal.broadinstitute.org/scalMDo/MitoCarta3.0/Human.MitoCarta3.0.xls" \
        -o "$MITOCARTA"
    echo "  -> $(du -h "$MITOCARTA" | cut -f1)"
fi

# -----------------------------------------------------------------------------
# 4. Bie et al. 2020 — HSP60/HSP10 SILAC interactome supplementary data
#    Bie AS, Cömert C, Körner R, Corydon TJ, Palmfeldt J, Hipp MS, Hartl FU, Bross P
#    "An inventory of interactors of the human HSP60/HSP10 chaperonin in the
#    mitochondrial matrix space."
#    Cell Stress and Chaperones 25(3):407-416 (2020).
#    DOI: 10.1007/s12192-020-01080-6 · PMID: 32060690 · PMC: PMC7192978
# -----------------------------------------------------------------------------
BIE_SUPPLEMENT="${PROJECT_DIR}/data/raw/custom/12192_2020_1080_MOESM4_ESM.xlsx"

if [ -f "$BIE_SUPPLEMENT" ]; then
    echo "[SKIP] Bie 2020 supplement already exists"
else
    echo "[DOWNLOAD] Bie et al. 2020 Supplementary Table S1..."
    echo "  NOTE: This file must be downloaded manually from the publisher."
    echo "  Source: Cell Stress and Chaperones 25(3):407-416 (2020) — DOI: 10.1007/s12192-020-01080-6"
    echo "  Direct link (may require institutional access):"
    echo "    https://static-content.springer.com/esm/art%3A10.1007%2Fs12192-020-01080-6/MediaObjects/12192_2020_1080_MOESM4_ESM.xlsx"
    echo "  Save to: ${BIE_SUPPLEMENT}"
fi

# -----------------------------------------------------------------------------
# 5. Kerner et al. 2005 — GroEL substrate data
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Kerner et al. 2005 data (Table S3) is included in the repository:"
echo "  -> data/raw/custom/kerner_2005_groel_interactors_table_s3.csv"
echo "  -> data/raw/custom/kerner_2005_groel_interactors_clean.csv"

# -----------------------------------------------------------------------------
# Verify checksums
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Verifying checksums"
echo "============================================================"

FAIL=0
verify_md5() {
    local file="$1"
    local expected="$2"
    local name="$3"
    if [ ! -f "$file" ]; then
        echo "  [SKIP] $name — file not found"
        return
    fi
    # macOS uses md5, Linux uses md5sum
    if command -v md5sum &>/dev/null; then
        actual=$(md5sum "$file" | awk '{print $1}')
    else
        actual=$(md5 -q "$file")
    fi
    if [ "$actual" = "$expected" ]; then
        echo "  [OK]   $name"
    else
        echo "  [WARN] $name — checksum mismatch (expected: $expected, got: $actual)"
        echo "         This may be due to a newer version from the source."
        FAIL=$((FAIL + 1))
    fi
}

# Checksums computed on 2026-03-25
verify_md5 "$ECOLI_FASTA" "d7621ce7710915756a8e230f3bc14978" "E. coli FASTA"
verify_md5 "$ECOLI_TSV"   "62362f32d93a96d9e72e6f6615cb31d7" "E. coli TSV"
verify_md5 "$HUMAN_FASTA" "f6c592be0e543cb46d2efc48e3a1e958" "Human FASTA"
verify_md5 "$HUMAN_TSV"   "26af2ecde9f5eac0a033c297a023d7e1" "Human TSV"
verify_md5 "$MITOCARTA"   "f113d5020ab277fb8869954aab773464" "MitoCarta 3.0"

echo ""
if [ $FAIL -gt 0 ]; then
    echo "WARNING: $FAIL checksum(s) did not match. Files may have been updated."
    echo "The analysis should still work, but protein counts may differ slightly."
else
    echo "All checksums verified."
fi

echo ""
echo "============================================================"
echo "  Download complete."
echo "  Next: Run Phase 1 pipeline (see Makefile or README.md)"
echo "============================================================"
