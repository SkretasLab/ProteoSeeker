#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Installation enviroment for ProteoSeeker
conda create -n ps_install -y
conda activate ps_install
conda install anaconda::git -y
conda install anaconda::wget -y
conda install conda-forge::gzip -y
conda install conda-forge::tar -y
conda deactivate