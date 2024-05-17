#!/bin/bash
source /home/gfilis/anaconda3_2024_02_1/bin/activate ps_diamond
cat profile_nr_families/fpr_databases/fpr_cas_als/*.fasta > profile_nr_families/fpr_databases/fpr_cas_als/fpr_cas_als.fasta
diamond makedb --threads 6 --in profile_nr_families/fpr_databases/fpr_cas_als/fpr_cas_als.fasta --db profile_nr_families/fpr_databases/fpr_cas_als/fpr_cas_als_database &> profile_nr_families/fpr_databases/fpr_cas_als/diamond_makedb_stdoe.txt
echo "No version available for diamond makedb." > profile_nr_families/fpr_databases/fpr_cas_als/diamond_makedb_version.txt
conda deactivate
