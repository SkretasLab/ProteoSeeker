#!/bin/bash
# COMEBin
conda create -n ps_kraken -y
source activate ps_kraken
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
# In any case create a directory for the kraken database and download the minkraken database.
cd ..
if [ ! -d ps_tools ]; then
  mkdir ps_tools
fi
cd ps_tools
if [ ! -d kraken2 ]; then
  mkdir kraken2
fi
cd kraken2
if [ ! -d kraken_databases ]; then
  mkdir kraken_databases
fi
mkdir kraken_databases
cd kraken_databases
# Download the minikraken 8GB database.
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz
tar -xvzf minikraken2_v2_8GB_201904.tgz
rm minikraken2_v2_8GB_201904.tgz
cd ../../../installation
