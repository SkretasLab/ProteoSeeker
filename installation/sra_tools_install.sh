#!/bin/bash
# SRA tools
conda create --name ps_sra_tools -y
source activate ps_sra_tools
conda install bioconda::sra-tools -y
conda deactivate
