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

# COMEBin - GPU
conda create -n ps_comebin_gpu -y
conda activate ps_comebin_gpu
conda install -c conda-forge -c bioconda comebin=1.0.4 -y
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia -c conda-forge
# Check if COMEBin was installed. If not, try installation from its git repository.
if [[ $(which run_comebin.sh) ]]; then
    echo "COMEBin with GPU support was installed successfully."
else
    echo "COMEBin with GPU support was not installed successfully."
    exit 1
fi
conda deactivate