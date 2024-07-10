#!/bin/bash
# HMMER
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_hmmer -y
conda activate ps_hmmer
conda install bioconda::hmmer -y
conda deactivate
