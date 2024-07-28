#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
KRAKEN2_DIR="${PS_TOOLS_DIR}/kraken2"
KRAKEN2_TOOL_DIR="${KRAKEN2_DIR}/kraken_tool"
KRAKEN2_SB_DIR="${KRAKEN2_DIR}/kraken2_sb"

# Kraken2
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create -n ps_kraken -y
conda activate ps_kraken
conda install bioconda::kraken2 -y
# Check if kraken2 was installed. If not, try installation from its git repository.
if [[ $(which kraken2) ]]; then
    echo "Kraken2 was installed successfully."
fi
conda deactivate

conda activate ps_kraken
if ! [[ $(which kraken2) ]]; then
    conda deactivate
    echo "Kraken2 was not installed successfully. Trying installing Kraken2 based on its git reporisory."
    # Create the ps_tools dir if needed.
    if [ ! -d "${PS_TOOLS_DIR}" ]; then
      mkdir "${PS_TOOLS_DIR}"
    fi
    # Create the kraken2 dir if needed.
    if [ ! -d "${KRAKEN2_DIR}" ]; then
      mkdir "${KRAKEN2_DIR}"
    fi
    # Create the kraken_tool dir if needed.
    if [ ! -d "${KRAKEN2_TOOL_DIR}" ]; then
      mkdir "${KRAKEN2_TOOL_DIR}"
    fi
    # Create the kraken2_sb dir if needed.
    if [ ! -d "${KRAKEN2_SB_DIR}" ]; then
      mkdir "${KRAKEN2_SB_DIR}"
    fi
    # Download.
    git clone https://github.com/DerrickWood/kraken2.git "${KRAKEN2_TOOL_DIR}"
    # Run installer.
    conda create -n ps_kraken -y
    conda activate ps_kraken
    "${KRAKEN2_TOOL_DIR}/install_kraken2.sh" "${KRAKEN2_SB_DIR}"
    # Add the Kraken2 installation directory to the PATH and check for the kraken2 executable.
    PATH=$PATH:${KRAKEN2_SB_DIR}
    if [[ $(which kraken2) ]]; then
      echo "kraken2 was installed successfully."
    else
      echo "kraken2 was not installed successfully."
    fi
    conda deactivate
else
    conda deactivate
fi