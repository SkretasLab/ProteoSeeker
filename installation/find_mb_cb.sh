#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
COMEBIN_GIT_DIR="${PS_TOOLS_DIR}/COMEBin"
COMEBIN_FILE="${COMEBIN_GIT_DIR}/COMEBin/run_comebin.sh"
METABINNER_GIT_DIR="${PS_TOOLS_DIR}/MetaBinner"
METABINNER_FILE="${METABINNER_GIT_DIR}/run_metabinner.sh"
CTBPATH=""
CBINPATH=""
MTBENV=""
MTBPATH=""
MBINPATH=""

# Source conda
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# MetaBinner enviroment
SCMET_ENV_CHECK="${CONDA_INST_DIR}/envs/metabinner_env"
PSMET_ENV_CHECK="${CONDA_INST_DIR}/envs/ps_metabinner"

# If when the enviroment is activated the shell script is found, then the enviroment was installed successfully. If not,
# then the installation based on the git repository is checked.
conda activate ps_comebin
if [[ $(which run_comebin.sh) ]]; then
    CTBPATH=$(which run_comebin.sh)
	CBINPATH="$(dirname "${CTBPATH}")" ; FILE="$(basename "${CTBPATH}")"
	CBINPATH="${CBINPATH}/COMEBin"
fi
conda deactivate
conda activate ps_comebin
if ! [[ $(which run_comebin.sh) ]]; then
	conda deactivate
	if [[ -f "${COMEBIN_FILE}" ]]; then
	    CTBPATH="${COMEBIN_FILE}"
		CBINPATH="${COMEBIN_GIT_DIR}/COMEBin"
	fi
else
	conda deactivate
fi

# If the metabinner_env enviroment exists.
if [ -d "${SCMET_ENV_CHECK}" ]; then
	# If the metabinner Bash script exists in the git directory.
	if [[ -f "${METABINNER_FILE}" ]]; then
		MTBENV="metabinner_env"
		MTBPATH="${METABINNER_FILE}"
		MBINPATH="${METABINNER_GIT_DIR}"
	else
		# If the metabinner Bash script does not exist in the git directory.
		# Check if the ps_metabinner enviroment exists.
		if [ -d "${PSMET_ENV_CHECK}" ]; then
			# If the metabinner Bash script exists in the ps_metabinner enviroment.
			if [[ $(which run_metabinner.sh) ]]; then
			    MTBPATH=$(which run_metabinner.sh)
				MBINPATH="$(dirname "${MTBPATH}")" ; FILE="$(basename "${MTBPATH}")"
			fi
		fi
	fi
else
	# If the metabinner_env enviroment does not exist
	# Check if the ps_metabinner enviroment exists.
	if [ -d "${PSMET_ENV_CHECK}" ]; then
		# If the metabinner Bash script exists in the ps_metabinner enviroment.
		if [[ $(which run_metabinner.sh) ]]; then
		  	MTBPATH=$(which run_metabinner.sh)
			MBINPATH="$(dirname "${MTBPATH}")" ; FILE="$(basename "${MTBPATH}")"
		fi
	fi
fi