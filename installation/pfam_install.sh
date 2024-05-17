#!/bin/bash
# Moving to the main folder.
cd ..
# Pfam database
if test -d pfam_database; then
  rm -r pfam_database
fi
mkdir pfam_database
cd pfam_database
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip -dk Pfam-A.hmm.gz
source activate ps_hmmer
hmmpress Pfam-A.hmm
conda deactivate
cd ..
