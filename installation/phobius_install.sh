#!/bin/bash
# TMHMM
conda create --name ps_phobius -y
source activate ps_phobius
conda install conda-forge::perl -y
conda deactivate
cd ../ps_tools
if [ ! -d phobius_files ]; then
  mkdir phobius_files
fi
cd ../installation