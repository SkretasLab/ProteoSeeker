#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
KRAKENTOOLS_GIT_DIR="${PS_TOOLS_DIR}/KrakenTools"

# Source conda.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Create the ps_tools dir if needed.
if [ ! -d "${PS_TOOLS_DIR}" ]; then
    mkdir "${PS_TOOLS_DIR}"
fi
# If the git directory has already been cloned delete it.
if [ -d "${KRAKENTOOLS_GIT_DIR}" ]; then
    rm -ri "${KRAKENTOOLS_GIT_DIR}"
fi
# Download the branch for a specific tag of the repository.
git clone --branch v1.2 https://github.com/jenniferlu717/KrakenTools.git "${KRAKENTOOLS_GIT_DIR}"