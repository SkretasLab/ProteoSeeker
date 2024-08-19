#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Diamond
conda create -n ps_diamond -y
conda activate ps_diamond
conda install bioconda::diamond=2.1.9 -y
conda deactivate

conda activate ps_diamond
if [[ $(which diamond) ]]; then
    echo "diamond was installed successfully."
else
	echo "diamond was not installed successfully."
	exit 1
fi
conda deactivate