#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# FastQC
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create --name ps_fastqc -y
conda activate ps_fastqc
conda install bioconda::fastqc -y
conda deactivate
