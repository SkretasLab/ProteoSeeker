#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Bowtie2
conda create -n ps_bowtie -y
conda activate ps_bowtie
conda install bioconda::bowtie2=2.5.3 -y
conda deactivate

conda activate ps_bowtie
if [[ $(which bowtie2) ]]; then
    echo "bowtie2 was installed successfully."
else
	echo "bowtie2 was not installed successfully."
	exit 1
fi
conda deactivate