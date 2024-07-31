#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Megahit
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create --name ps_megahit -y
conda activate ps_megahit
conda install -c bioconda megahit=1.2.9 -y
conda deactivate

conda activate ps_megahit
if [[ $(which megahit) ]]; then
    echo "megahit was installed successfully."
else
	echo "megahit was not installed successfully."
fi
conda deactivate