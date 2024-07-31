#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PFAM_DB_DIR="${PS_DIR}/pfam_database"
PFAM_DB_GZ_FILE="${PFAM_DB_DIR}/Pfam-A.hmm.gz"
PFAM_DB_FILE="${PFAM_DB_DIR}/Pfam-A.hmm"

# Create the pfam_database dir if needed.
if [ ! -d "${PFAM_DB_DIR}" ]; then
  mkdir "${PFAM_DB_DIR}"
fi
# Download.
wget -P "${PFAM_DB_DIR}" https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# Decompress.
gunzip -dc "${PFAM_DB_GZ_FILE}" > "${PFAM_DB_FILE}"
# Press the database.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda activate ps_hmmer
hmmpress "${PFAM_DB_FILE}"
conda deactivate