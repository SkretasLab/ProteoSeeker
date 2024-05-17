#!/bin/bash
# FastQC
conda create --name ps_fastqc -y
source activate ps_fastqc
conda install bioconda::fastqc -y
conda deactivate
