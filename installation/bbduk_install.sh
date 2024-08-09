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

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# BBDuk
conda create --name ps_bbtools -y
conda activate ps_bbtools
# If conda fails the following installation, this script wont exit. Hence, the subsequent installation will be attempted.
(
    conda install bioconda::bbmap=39.01 -y
)
# Check if bbduk.sh was installed. If not, try installation from its git repository.
if [[ $(which bbduk.sh) ]]; then
    echo "bbtools were installed successfully."
fi
conda deactivate

# Trying another conda package.
conda activate ps_bbtools
if ! [[ $(which bbduk.sh) ]]; then
    echo "bbtools were not installed successfully. Trying installing bbtools from the agbiome channel."
    # If conda fails the following installation, this script wont exit. Hence, the subsequent installation will be attempted.
    (
        conda install agbiome::bbtools -y
    )
    if [[ $(which bbduk.sh) ]]; then
        echo "bbtools were installed successfully."
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
    # Download the latest version from sourceforge.
    # wget -P "${BBTOOLS_DIR}" 'https://sourceforge.net/projects/bbmap/files/latest/download'
    # Download a specific version from sourceforge.
    wget -P "${BBTOOLS_DIR}" 'https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download'
    # Rename.
    mv "${BBTOOLS_DN_FILE}" "${BBTOOLS_DN_TG_FILE}"
    # Decompress.
    tar -xvzf "${BBTOOLS_DN_TG_FILE}" -C "${BBTOOLS_DIR}"
    # Check for the shell script.
    if [[ -f "${BBDUK_SH_FILE}" ]]; then
        echo "bbtools were installed successfully."
    else
        echo "bbtools were not installed successfully."
        exit 1
    fi
else
    conda deactivate
fi