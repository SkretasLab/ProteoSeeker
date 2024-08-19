#!/bin/bash
set -e

# Check if conda is activated in the shell.
if ! command -v conda &> /dev/null
then
    echo "Conda was not found. Exiting."
    exit 1
fi
# Find the conda.sh full path.
CONDA_SCRIPTS=$(which conda)
if [ -z "${CONDA_SCRIPTS}" ]; then
    echo "Conda is not activated in the shell. Exiting."
    exit 1
fi
# Find the path to the installation directory of conda.
CONDA_INST_DIR=$(dirname $(dirname "${CONDA_SCRIPTS}"))
# Determine the conda.sh path.
CONDA_SH_PATH="${CONDA_INST_DIR}/etc/profile.d/conda.sh"
if [ ! -f "${CONDA_SH_PATH}" ]; then
    echo "The path to conda.sh could not be determined."
    exit 1
fi
