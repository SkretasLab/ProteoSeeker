#!/bin/bash
# SRA tools
source ./find_conda.sh
source $CONDA_SH_PATH
conda create --name ps_sra_tools -y
conda activate ps_sra_tools
conda install bioconda::sra-tools -y
conda deactivate
