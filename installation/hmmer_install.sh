#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# HMMER
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create -n ps_hmmer -y
conda activate ps_hmmer
conda install bioconda::hmmer=3.4 -y
conda deactivate

conda activate ps_hmmer
if [[ $(which hmmscan) ]]; then
    echo "hmmer was installed successfully."
else
	echo "hmmer was not installed successfully."
fi
conda deactivate