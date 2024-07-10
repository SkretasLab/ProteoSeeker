#!/bin/bash
# Megahit
source ./find_conda.sh
source $CONDA_SH_PATH
conda create --name ps_megahit -y
conda activate ps_megahit
conda install -c bioconda megahit -y
conda deactivate
