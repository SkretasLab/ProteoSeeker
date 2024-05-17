#!/bin/bash
# COMEBin
conda create -n ps_comebin -y
source activate ps_comebin
conda install -c conda-forge -c bioconda comebin -y
# Check if COMEBin was installed. If not, try installation from its git repository.
if [[ $(which run_comebin.sh) ]]; then
    echo "COMEBin was installed successfully."
    conda deactivate
else
    echo "COMEBin was not installed successfully. Trying installation based on it git reporisory."
    conda deactivate
    cd ..
    if [ ! -d ps_tools ]; then
      mkdir ps_tools
    fi
    cd ps_tools
    mkdir comebin
    cd comebin
    git clone https://github.com/ziyewang/COMEBin.git
    cd COMEBin
    conda env create -f comebin_env.yaml
    conda activate comebin_env
    if [[ $(which run_comebin.sh) ]]; then
      echo "COMEBin was installed successfully."
    else
      echo "COMEBin was not installed successfully."
    fi
    conda deactivate
    cd ..
    cd ..
    cd installation
fi

