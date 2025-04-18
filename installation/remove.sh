#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
TAXONOMY_DB_DIR="${PS_DIR}/taxonomy_db_info"
TAXONOMY_HOME_DB_DIR="${HOME}/.taxonkit"
PFAM_DB_DIR="${PS_DIR}/pfam_database"
SWISSPROT_DIR="${PS_DIR}/swissprot_database"
UNIREF_DBS_DIR="${PS_DIR}/uniref_db"

# An array of the installation environments.
install_envs=("ps_bbtools" "ps_bowtie" "ps_cd_hit" "ps_comebin" "ps_diamond" "ps_env" "ps_fastqc" "ps_fraggenescan" "ps_hmmer" "ps_install" "ps_kraken" "ps_bracken" "ps_megahit" "metabinner_env" "ps_metabinner" "ps_phobius" "ps_samtools" "ps_sra_tools" "ps_taxonkit", "ps_result_analysis")

# An array of the installation directories.
install_dirs=("${PS_TOOLS_DIR}" "${TAXONOMY_DB_DIR}" "${TAXONOMY_HOME_DB_DIR}" "${PFAM_DB_DIR}" "${SWISSPROT_DIR}" "${UNIREF_DBS_DIR}")

# Check for the conda installation.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# User input.
read -p "Delete all installation environments for ProteoSeeker? Type y/yes to delete or n/no otherwise: " ENV_DEL_STATUS
if [ ${ENV_DEL_STATUS} = "y" ] || [ ${ENV_DEL_STATUS} = "yes" ]; then
	echo "Deleting ProteoSeeker environments..."
	# Check if each installation environment exists. If yes, delete it.
	for ITEM_ENV in "${install_envs[@]}"; do
		if conda info --envs | grep -q "${ITEM_ENV}"; then
			echo "Deleting ProteoSeeker environment ${ITEM_ENV}..."
			conda remove --name "${ITEM_ENV}" --all -y
		fi
	done
fi

# User input.
read -p "Delete all installation directories for ProteoSeeker? Type y/yes to delete or n/no otherwise: " DIR_DEL_STATUS
# Check if each installation directory exists. If yes delete it.
if [ ${DIR_DEL_STATUS} = "y" ] || [ ${DIR_DEL_STATUS} = "yes" ]; then
	echo "Deleting ProteoSeeker directories..."
	for ITEM_DIR in "${install_dirs[@]}"; do
		if [ -d "${ITEM_DIR}" ]; then
			echo "Deleting ProteoSeeker directory ${ITEM_DIR}..."
			rm -r "${ITEM_DIR}"
		fi
	done
fi