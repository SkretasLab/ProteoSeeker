#!/bin/bash
# Megahit
conda create --name ps_megahit -y
source activate ps_megahit
conda install -c bioconda megahit -y
conda deactivate
