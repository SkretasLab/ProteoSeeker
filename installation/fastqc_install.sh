#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# FastQC
conda create --name ps_fastqc -y
conda activate ps_fastqc
conda install bioconda::fastqc=0.12.1 -y
conda deactivate

conda activate ps_fastqc
if [[ $(which fastqc) ]]; then
    echo "fastqc was installed successfully."
else
	echo "fastqc was not installed successfully."
	exit 1
fi
conda deactivate