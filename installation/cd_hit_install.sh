#!/bin/bash
# CD-HIT
conda create -n ps_cd_hit -y
source activate ps_cd_hit
conda install bioconda::cd-hit -y
conda deactivate
