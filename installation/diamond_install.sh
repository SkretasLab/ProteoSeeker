#!/bin/bash
# Diamond
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_diamond -y
conda activate ps_diamond
conda install bioconda::diamond -y
conda deactivate
