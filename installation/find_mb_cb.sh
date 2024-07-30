#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
COMEBIN_GIT_DIR="${PS_TOOLS_DIR}/COMEBin"
COMEBIN_FILE="${COMEBIN_GIT_DIR}/COMEBin/run_comebin.sh"
METABINNER_GIT_DIR="${PS_TOOLS_DIR}/MetaBinner"
METABINNER_FILE="${METABINNER_GIT_DIR}/run_metabinner.sh"

source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# If when the enviroment is activated the shell script is found, then the enviroment was installed successfully. If not,
# then the installation based on the git repository is checked.
if [[ $(which run_comebin.sh) ]]; then
    CTBPATH=$(which run_comebin.sh)
	CBINPATH="$(dirname "${CTBPATH}")" ; FILE="$(basename "${CTBPATH}")"
	conda deactivate
elif [[ -f "${COMEBIN_FILE}" ]]; then
    CTBPATH="${COMEBIN_FILE}"
	CBINPATH="${COMEBIN_GIT_DIR}"
fi

# If when the enviroment is activated the shell script is found, then the enviroment was installed successfully. If not,
# then the installation based on the git repository is checked.
conda activate ps_metabinner
if [[ $(which run_metabinner.sh) ]]; then
    MTBPATH=$(which run_metabinner.sh)
	MBINPATH="$(dirname "${MTBPATH}")" ; FILE="$(basename "${MTBPATH}")"
	conda deactivate
elif [[ -f "${METABINNER_FILE}" ]]; then
    MTBPATH="${METABINNER_FILE}"
	MBINPATH="${METABINNER_GIT_DIR}"
fi