#!/bin/bash
# BBDuk
source ./find_conda.sh
source $CONDA_SH_PATH
conda create --name ps_bbtools -y
conda activate ps_bbtools
conda install bioconda::bbmap -y
# Check if bbduk.sh was installed. If not, try installation from its git repository.
if [[ $(which bbduk.sh) ]]; then
    echo "bbduk.sh was installed successfully."
    conda deactivate
else
    echo "bbduk.sh was not installed successfully. Trying installation based on it git reporisory."
    conda install agbiome::bbtools -y

if [[ $(which bbduk.sh) ]]; then
    echo "bbduk.sh was installed successfully."
    conda deactivate
else    
    conda deactivate
    cd ..
    if [ ! -d ps_tools ]; then
      mkdir ps_tools
    fi
    cd ps_tools
    mkdir bbtools
    cd bbtools
    # Needs testing.
    #wget 'https://sourceforge.net/projects/bbmap/files/latest/download'
    if [[ $(which bbduk.sh) ]]; then
      echo "bbduk.sh was installed successfully."
    else
      echo "bbduk.sh was not installed successfully."
    fi
    conda deactivate
    cd ..
    cd ..
    cd installation
fi