#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
TAXONOMY_DB_DIR="${INSTALLATION_DIR}/taxonomy_db_info"
TAXONOMY_TAR_GZ_FILE="${TAXONOMY_DB_DIR}/taxdump.tar.gz"
TAXONOMY_HOME_DB_DIR="${HOME}/.taxonkit"

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# TaxonKit
conda create -n ps_taxonkit -y
conda activate ps_taxonkit
conda install -c bioconda taxonkit=0.16.0 -y
conda install bioconda::csvtk=0.30.0 -y
conda deactivate

conda activate ps_taxonkit
if [[ $(which taxonkit) ]]; then
    echo "taxonkit was installed successfully."
else
    echo "taxonkit was not installed successfully."
    exit 1
fi
conda deactivate

# Create the taxonomy_db_info directory if needed.
if [ ! -d "${TAXONOMY_DB_DIR}" ]; then
    mkdir "${TAXONOMY_DB_DIR}"
fi
# Download.
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -P "${TAXONOMY_DB_DIR}" ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# Decompress
tar -zxvf "${TAXONOMY_TAR_GZ_FILE}" -C "${TAXONOMY_DB_DIR}"
# Create the taxonomy directory in the home directory.
mkdir -p "${TAXONOMY_HOME_DB_DIR}"
# Copy the files.
cp "${TAXONOMY_DB_DIR}/names.dmp" "${TAXONOMY_DB_DIR}/nodes.dmp" "${TAXONOMY_DB_DIR}/delnodes.dmp" "${TAXONOMY_DB_DIR}/merged.dmp" "${TAXONOMY_HOME_DB_DIR}"