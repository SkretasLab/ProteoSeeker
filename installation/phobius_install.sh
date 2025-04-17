#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
PHOBIUS_DIR="${PS_TOOLS_DIR}/phobius_files"

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Phobius environment
conda create --name ps_phobius -y
conda activate ps_phobius
conda install conda-forge::perl=5.32.1 -y
conda deactivate

# Create the ps_tools directory if needed.
if [ ! -d "${PS_TOOLS_DIR}" ]; then
    mkdir "${PS_TOOLS_DIR}"
fi
# Create the phobius_files directory if needed.
if [ ! -d "${PHOBIUS_DIR}" ]; then
    mkdir "${PHOBIUS_DIR}"
fi