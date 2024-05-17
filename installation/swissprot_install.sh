#!/bin/bash
# Moving to the main folder.
cd ..
# Swiss-Prot/UniprotKB database
if [ ! -d ps_tools ]; then
  mkdir swissprot_database
fi
cd swissprot_database
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
# DIAMOND database
if test -d swissprot_database; then
  rm -r swissprot_dtb
fi
mkdir swissprot_dtb
cd swissprot_dtb
gunzip -dkc ../swissprot.gz > swissprot
source activate ps_diamond
diamond makedb --in swissprot --db swissprot_db
conda deactivate
cd ..
cd ..
