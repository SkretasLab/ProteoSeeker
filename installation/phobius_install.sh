#!/bin/bash
# Phobius
source ./find_conda.sh
source $CONDA_SH_PATH
conda create --name ps_phobius -y
conda activate ps_phobius
conda install conda-forge::perl -y
conda deactivate
cd ../ps_tools
if [ ! -d phobius_files ]; then
  mkdir phobius_files
fi
cd ../installation
