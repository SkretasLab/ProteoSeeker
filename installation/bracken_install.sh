#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
BRACKEN_GIT_DIR="${PS_TOOLS_DIR}/Bracken"
BRACKEN_FILE="${BRACKEN_GIT_DIR}/bracken"

# Source conda.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Create the ps_tools dir if needed.
if [ ! -d "${PS_TOOLS_DIR}" ]; then
    mkdir "${PS_TOOLS_DIR}"
fi

# Create an environment.
conda create -n ps_bracken -y
conda activate ps_bracken
(
    conda install bioconda::bracken=2.9 -y
)
if [[ $(which bracken) ]]; then
    echo "Bracken was installed successfully."
else
    echo "Bracken was not installed successfully. Attempting installation based on its source code."
fi
conda deactivate

conda activate ps_bracken
if ! [[ $(which bracken) ]]; then
    conda deactivate
    # If the git directory has already been cloned delete it.
    if [ -d "${BRACKEN_GIT_DIR}" ]; then
        rm -ri "${BRACKEN_GIT_DIR}"
    fi
    # Download the head branch of the repository.
    git clone --branch v2.9 https://github.com/jenniferlu717/Bracken.git "${BRACKEN_GIT_DIR}" 
    # Moving to the git directory for bracken and installing bracken. Then, moving back
    # to the installation directory.
    cd "${BRACKEN_GIT_DIR}"
    bash install_bracken.sh
    cd "${INSTALLATION_DIR}"

    # Check for the installation.
    if [[ ! -f "${BRACKEN_FILE}" ]]; then
        echo "Bracken was not installed successfully."
        exit 1
    else
        echo "Bracken was installed successfully."
    fi
else
    conda deactivate
fi