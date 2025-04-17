#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
SWISSPROT_DIR="${PS_DIR}/swissprot_database"
SWISSPROT_GZ_FILE="${SWISSPROT_DIR}/swissprot.gz"
SWISSPROT_FILE="${SWISSPROT_DIR}/swissprot"
SWISSPROT_DB_FILE="${SWISSPROT_DIR}/swissprot_db"

# Create the swissprot_database directory if needed.
if [ ! -d "${SWISSPROT_DIR}" ]; then
    mkdir "${SWISSPROT_DIR}"
fi
# Download.
wget -P "${SWISSPROT_DIR}" https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
# Decompress.
gunzip -dc "${SWISSPROT_GZ_FILE}" > "${SWISSPROT_FILE}"
# Create a DIAMOND database
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda activate ps_diamond
diamond makedb --in "${SWISSPROT_FILE}" --db "${SWISSPROT_DB_FILE}"
conda deactivate
# Remove the compressed file.
rm "${SWISSPROT_GZ_FILE}"