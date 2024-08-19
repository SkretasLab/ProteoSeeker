#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
COMEBIN_GIT_DIR="${PS_TOOLS_DIR}/COMEBin"
COMEBIN_FILE="${COMEBIN_GIT_DIR}/COMEBin/run_comebin.sh"

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# COMEBin
conda create -n ps_comebin -y
conda activate ps_comebin
# If conda fails the following installation, this script wont exit. Hence, the subsequent installation will be attempted.
(
    conda install -c conda-forge -c bioconda comebin=1.0.4 -y
)
# Check if COMEBin was installed. If not, try installation from its git repository.
if [[ $(which run_comebin.sh) ]]; then
    echo "COMEBin was installed successfully."
fi
conda deactivate

# Trying by source.
conda activate ps_comebin
if ! [[ $(which run_comebin.sh) ]]; then
    conda deactivate
    echo "COMEBin was not installed successfully. Trying installation based on it git reporisory."
    # Create the ps_tools dir if needed.
    if [ ! -d "${PS_TOOLS_DIR}" ]; then
        mkdir "${PS_TOOLS_DIR}"
    fi
    # Download the head branch of the repository.
    # Download a specific branch.
    git clone --branch 1.0.4 https://github.com/ziyewang/COMEBin.git "${COMEBIN_GIT_DIR}"
    # Move to the git clone.
    cd "${COMEBIN_GIT_DIR}"
    # Install dependencies.
    conda env create -f comebin_env.yaml
    conda activate comebin_env
    if [[ -f "${COMEBIN_FILE}" ]]; then
        echo "COMEBin was installed successfully."
    else
        echo "COMEBin was not installed successfully."
        exit 1
    fi
    conda deactivate
    cd "${INSTALLATION_DIR}"
else
    conda deactivate
fi