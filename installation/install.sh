#!/bin/bash
set -e

# Get the current directory.
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INSTALL_OUTPUT_DIR="${INSTALLATION_DIR}/output_errors"
PS_DIR=$(dirname "${INSTALLATION_DIR}")
INPUT_DIR="${PS_DIR}/input"
RESULTS_DIR="${PS_DIR}/results"
FPD_PD_DIR="${PS_DIR}/profile_protein_dbs"
FPD_DIR="${FPD_PD_DIR}/filtered_protein_dbs"
PD_DIR="${FPD_PD_DIR}/profile_dbs"
SRA_DIR="${PS_DIR}/sra_files"

# Create directories, if needed.
if [ ! -d "${INSTALL_OUTPUT_DIR}" ]; then
    mkdir "${INSTALL_OUTPUT_DIR}"
fi
if [ ! -d "${INPUT_DIR}" ]; then
    mkdir "${INPUT_DIR}"
fi
if [ ! -d "${RESULTS_DIR}" ]; then
    mkdir "${RESULTS_DIR}"
fi
if [ ! -d "${SRA_DIR}" ]; then
    mkdir "${SRA_DIR}"
fi
if [ ! -d "${FPD_DIR}" ]; then
    mkdir "${FPD_DIR}"
fi
if [ ! -d "${PD_DIR}" ]; then
    mkdir "${PD_DIR}"
fi

# Check for the conda installation.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Create the installation environment.
echo -e "\nSetting up the installation environment for ProteoSeeker..."
"${INSTALLATION_DIR}/set_install_env.sh" &> "${INSTALLATION_DIR}/output_errors/set_install_stdoe.txt"

# Run the installers.
conda activate ps_install
"${INSTALLATION_DIR}/install_envs.sh"
"${INSTALLATION_DIR}/install_dbs.sh"
conda deactivate