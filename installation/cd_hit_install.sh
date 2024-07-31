#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# CD-HIT
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH
conda create -n ps_cd_hit -y
conda activate ps_cd_hit
conda install bioconda::cd-hit=4.8.1 -y
conda deactivate

conda activate ps_cd_hit
if [[ $(which cd-hit) ]]; then
    echo "cd-hit was installed successfully."
else
	echo "cd-hit was not installed successfully."
fi
conda deactivate