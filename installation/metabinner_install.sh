#!/bin/bash
# MetaBinner
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_metabinner python=3.7.6 -y
conda activate ps_metabinner
conda install -c bioconda metabinner -y
# Check if MetaBinner was installed. If not, try installation from its git repository.
if [[ $(which run_metabinner.sh) ]]; then
    echo "MetaBinner was installed successfully."
    conda deactivate
else
    echo "MetaBinner was not installed successfully. Trying installation based on it git reporisory."
    conda deactivate
    cd ..
    if [ ! -d ps_tools ]; then
      mkdir ps_tools
    fi
    cd ps_tools
    mkdir metabinner
    cd metabinner
    git clone https://github.com/ziyewang/MetaBinner.git
    cd MetaBinner
    conda env create -f metabinner_env.yaml
    conda activate metabinner_env
    # Need to check for the specific path of ps_tools/metabinner/MetaBinner/run_metabinner.sh if the file exists or not.
    # echo "MetaBinner was installed successfully."
    # echo "MetaBinner was not installed successfully."
    conda deactivate
    cd ..
    cd ..
    cd installation
fi
