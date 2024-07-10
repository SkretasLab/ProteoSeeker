#!/bin/bash
# FastQC
source ./find_conda.sh
source $CONDA_SH_PATH
conda create --name ps_fastqc -y
conda activate ps_fastqc
conda install bioconda::fastqc -y
conda deactivate
