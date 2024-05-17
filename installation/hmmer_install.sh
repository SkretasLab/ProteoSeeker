#!/bin/bash
# HMMER
conda create -n ps_hmmer -y
source activate ps_hmmer
conda install bioconda::hmmer -y
conda deactivate
