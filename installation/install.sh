#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "${INSTALLATION_DIR}/find_conda.sh"
# Run the installers.
"${INSTALLATION_DIR}/install_envs.sh"
"${INSTALLATION_DIR}/install_dbs.sh"
"${INSTALLATION_DIR}/parameter_files.sh"