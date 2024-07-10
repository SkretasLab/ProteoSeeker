#!/bin/bash
# Kraken2
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_kraken -y
conda activate ps_kraken
conda install bioconda::kraken2 -y
# Check if kraken2 was installed. If not, try installation from its git repository.
if [[ $(which kraken2) ]]; then
    echo "kraken2 was installed successfully."
    conda deactivate
else
    echo "kraken2 was not installed successfully. Trying installation based on it git reporisory."
    conda deactivate
    cd ..
    if [ ! -d ps_tools ]; then
      mkdir ps_tools
    fi
    cd ps_tools
    mkdir kraken_tool
    cd kraken_tool
    mkdir kraken2_sb
    git clone https://github.com/DerrickWood/kraken2.git
    cd kraken2
    conda create -n ps_kraken -y
    conda activate ps_kraken
    ./install_kraken2.sh ps_tools/kraken_tool/kraken2_sb
    cd ../kraken2_sb
    cur_path=$(pwd)
    PATH=$PATH:${cur_path}
    if [[ $(which kraken2) ]]; then
      echo "kraken2 was installed successfully."
    else
      echo "kraken2 was not installed successfully."
    fi
    conda deactivate
    cd ../../installation
fi