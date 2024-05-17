#!/bin/bash
# Diamond
conda create -n ps_diamond -y
source activate ps_diamond
conda install bioconda::diamond -y
conda deactivate
