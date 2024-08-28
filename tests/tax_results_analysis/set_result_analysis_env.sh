#!/bin/bash
set -e

# Get the current directory.
TAX_ANALYSIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the tests directory of ProteoSeeker.
TESTS_DIR=$(dirname "${TAX_ANALYSIS_DIR}")
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${TESTS_DIR}")
# Get the path of the installation directory.
INSTALLATION_DIR="${PS_DIR}/installation"

# Source conda.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

conda create -n result_analysis python=3.9.7 -y
conda activate result_analysis
conda install anaconda::numpy=1.26.4 -y
conda install anaconda::pandas=2.2.2 -y
conda install anaconda::scikit-learn=1.5.1 -y
conda install conda-forge::matplotlib=3.9.2 -y
conda install anaconda::seaborn=0.13.2 -y
conda deactivate