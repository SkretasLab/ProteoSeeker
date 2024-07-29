#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check for the conda installation.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Create the installation enviroment.
"${INSTALLATION_DIR}/set_install_env.sh" &> "${INSTALLATION_DIR}/output_errors/set_install_stdoe.txt"

# Run the installers.
conda activate ps_install
"${INSTALLATION_DIR}/install_envs.sh"
"${INSTALLATION_DIR}/install_dbs.sh"
"${INSTALLATION_DIR}/parameter_files.sh"
conda deactivate