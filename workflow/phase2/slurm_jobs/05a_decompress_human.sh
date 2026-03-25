#!/bin/bash
#SBATCH --job-name=aap_gunzip
#SBATCH --output=logs/05a_decompress_%j.out
#SBATCH --error=logs/05a_decompress_%j.err
#SBATCH --partition=compute
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# =============================================================================
# Decompress Human AlphaFold CIF files and remove unneeded PDB.gz files
# =============================================================================

set -euo pipefail

HUMAN_DIR="/lustre/vishal.bharti/Antah_Asti_Prarambh_hpc/data/raw/alphafold/full/human"

echo "============================================================"
echo "  Decompress Human AlphaFold CIFs — $(date)"
echo "============================================================"

# Count files before
N_GZ=$(ls "${HUMAN_DIR}"/*.cif.gz 2>/dev/null | wc -l)
N_CIF=$(ls "${HUMAN_DIR}"/*.cif 2>/dev/null | wc -l)
N_PDB_GZ=$(ls "${HUMAN_DIR}"/*.pdb.gz 2>/dev/null | wc -l)
echo "Before: ${N_GZ} .cif.gz, ${N_CIF} .cif, ${N_PDB_GZ} .pdb.gz"

# Remove PDB.gz files (we only need CIF format)
if [ "${N_PDB_GZ}" -gt 0 ]; then
    echo "Removing ${N_PDB_GZ} .pdb.gz files (not needed)..."
    rm -f "${HUMAN_DIR}"/*.pdb.gz
    echo "Done removing .pdb.gz files."
fi

# Decompress CIF.gz files
if [ "${N_GZ}" -gt 0 ]; then
    echo "Decompressing ${N_GZ} .cif.gz files..."
    cd "${HUMAN_DIR}"
    gunzip -f *.cif.gz
    echo "Done decompressing."
fi

# Count files after
N_CIF_AFTER=$(ls "${HUMAN_DIR}"/*.cif 2>/dev/null | wc -l)
echo ""
echo "After: ${N_CIF_AFTER} .cif files"
echo "Finished: $(date)"
echo "============================================================"
