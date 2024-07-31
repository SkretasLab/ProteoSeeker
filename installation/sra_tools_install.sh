#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# SRA tools
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create --name ps_sra_tools -y
conda activate ps_sra_tools
conda install bioconda::sra-tools=3.1.0 -y
conda deactivate

conda activate ps_sra_tools
if [[ $(which prefetch) ]]; then
    echo "sra tools were installed successfully."
else
	echo "sra tools were not installed successfully."
fi
conda deactivate