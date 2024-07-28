#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
BBTOOLS_DIR="${PS_TOOLS_DIR}/bbtools"
BBTOOLS_DN_FILE="${BBTOOLS_DIR}/download"
BBTOOLS_DN_TG_FILE="${BBTOOLS_DIR}/download.tar.gz"
BBDUK_SH_FILE="${BBTOOLS_DIR}/bbmap/bbduk.sh"

# BBDuk
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create --name ps_bbtools -y
conda activate ps_bbtools
conda install bioconda::bbmap -y
# Check if bbduk.sh was installed. If not, try installation from its git repository.
if [[ $(which bbduk.sh) ]]; then
    echo "bbduk.sh was installed successfully."
fi
conda deactivate

# Trying another conda package.
conda activate ps_bbtools
if ! [[ $(which bbduk.sh) ]]; then
    echo "bbduk.sh was not installed successfully. Trying another approach."
    conda install agbiome::bbtools -y
    if [[ $(which bbduk.sh) ]]; then
        echo "bbduk.sh was installed successfully."
    fi
fi
conda deactivate

# Trying by source.
conda activate ps_bbtools
if ! [[ $(which bbduk.sh) ]]; then
    conda deactivate
    # Create the ps_tools dir if needed.
    if [ ! -d "${PS_TOOLS_DIR}" ]; then
      mkdir "${PS_TOOLS_DIR}"
    fi
    # Create the bbtools dir if needed.
    if [ ! -d "${BBTOOLS_DIR}" ]; then
      mkdir "${BBTOOLS_DIR}"
    fi
    # Download.
    wget -P "${BBTOOLS_DIR}" 'https://sourceforge.net/projects/bbmap/files/latest/download'
    # Rename.
    mv "${BBTOOLS_DN_FILE}" "${BBTOOLS_DN_TG_FILE}"
    # Decompress.
    tar -xvzf "${BBTOOLS_DN_TG_FILE}" -C "${BBTOOLS_DIR}"
    # Check for the shell script.
    if [[ -f "${BBDUK_SH_FILE}" ]]; then
      echo "bbduk.sh was installed successfully."
    else
      echo "bbduk.sh was not installed successfully."
    fi
else
    conda deactivate
fi