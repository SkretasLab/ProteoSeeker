#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
METABINNER_GIT_DIR="${PS_TOOLS_DIR}/MetaBinner"
METABINNER_FILE="${METABINNER_GIT_DIR}/run_metabinner.sh"

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# MetaBinner
conda create -n ps_metabinner python=3.7.6 -y
conda activate ps_metabinner
# If conda fails the following installation, this script wont exit. Hence, the subsequent installation will be attempted.
(
    conda install -c bioconda metabinner=1.4.4 -y
)
# Check if MetaBinner was installed. If not, try installation from its git repository.
if [[ $(which run_metabinner.sh) ]]; then
    echo "MetaBinner was installed successfully."
fi
conda deactivate

conda activate ps_metabinner
if ! [[ $(which run_metabinner.sh) ]]; then
    conda deactivate
    echo "MetaBinner was not installed successfully. Trying installing MetaBinner based on its git reporisory."
    # Create the ps_tools dir if needed.
    if [ ! -d "${PS_TOOLS_DIR}" ]; then
        mkdir "${PS_TOOLS_DIR}"
    fi
    # Download the head branch of the repository.
    # git clone https://github.com/ziyewang/MetaBinner.git "${PS_TOOLS_DIR}"
    # Download a specific release.
    git clone https://github.com/ziyewang/MetaBinner.git "${PS_TOOLS_DIR}" --branch v1.4.4
    # Check for the installation.
    conda env create -f metabinner_env.yaml
    conda activate metabinner_env
    if [[ -f "${METABINNER_FILE}" ]]; then
        echo "MetaBinner was installed successfully."
    else
        echo "MetaBinner was not installed successfully."
        exit 1
    fi
    conda deactivate
else
    conda deactivate
fi