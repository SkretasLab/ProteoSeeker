#!/bin/bash
# Moving to the main folder.
cd ..
# Swiss-Prot/UniprotKB database
if [ ! -d swissprot_database ]; then
  mkdir swissprot_database
fi
cd swissprot_database
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
gunzip -dk swissprot.gz
# Create a DIAMOND database
source ./find_conda.sh
source $CONDA_SH_PATH
conda activate ps_diamond
diamond makedb --in swissprot --db swissprot_db
conda deactivate
cd ..