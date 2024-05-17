#!/bin/bash
# Bowtie2
conda create -n ps_bowtie -y
source activate ps_bowtie
conda install bioconda::bowtie2 -y
conda deactivate
