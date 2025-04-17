#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Base environment for ProteoSeeker with: Python, Openpyxl
conda create -n ps_env python=3.9.7 -y
conda activate ps_env
python -m pip install openpyxl==3.1.2
python -c "import openpyxl"
if [ ! $? ]; then
    echo "The python package openpyxl was not installed successfully. Trying installation through anaconda."
    conda install anaconda::openpyxl=3.1.2 -y
fi
python -c "import openpyxl"
if [ $? ]; then
    echo "The python package openpyxl was installed successfully."
else
    echo "The python package openpyxl was not installed successfully."
    exit 1
fi
conda deactivate
