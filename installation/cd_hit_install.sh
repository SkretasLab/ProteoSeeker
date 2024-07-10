#!/bin/bash
# CD-HIT
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_cd_hit -y
conda activate ps_cd_hit
conda install bioconda::cd-hit -y
conda deactivate
