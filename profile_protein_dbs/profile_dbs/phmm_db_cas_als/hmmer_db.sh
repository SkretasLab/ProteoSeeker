#!/bin/bash
source /home/gfilis/anaconda3_2024_02_1/bin/activate ps_hmmer
hmmfetch -o profile_nr_families/phmms_databases/phmm_db_cas_als/profiles -f pfam_database/Pfam-A.hmm profile_nr_families/phmms_databases/phmm_db_cas_als/profile_names.txt &> profile_nr_families/phmms_databases/phmm_db_cas_als/hmmfetch_stdoe.txt
hmmpress profile_nr_families/phmms_databases/phmm_db_cas_als/profiles &> profile_nr_families/phmms_databases/phmm_db_cas_als/hmmpress_stdoe.txt
hmmfetch -h > profile_nr_families/phmms_databases/phmm_db_cas_als/hmmfetch_version.txt
hmmpress -h > profile_nr_families/phmms_databases/phmm_db_cas_als/hmmpress_version.txt
conda deactivate
