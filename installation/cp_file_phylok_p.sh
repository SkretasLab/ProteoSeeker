#!/bin/bash
# Phylo - Kraken - Paired

# Sourcing
source ./find_conda.sh
source ./find_mb_cb.sh
source ./find_pd.sh

# Moving to main directory.
cd ..

# Creating the file for the input paths and option and writing in it.
echo "#---------Input and output options---------" > parameters_phylok_p.txt
echo "input_folder=\"input\"" >> parameters_phylok_p.txt
echo "sra_code=\"SRR12829170\"" >> parameters_phylok_p.txt
echo "contigs=\"\"" >> parameters_phylok_p.txt
echo "protein_input=\"\"" >> parameters_phylok_p.txt
echo "adapters_path=\"adapters.fa\"" >> parameters_phylok_p.txt
echo "protein_db_path=\"${PD_PATH}\"" >> parameters_phylok_p.txt
echo "kraken_db_path=\"ps_tools/kraken2/kraken_databases/minikraken2_v2_8GB_201904_UPDATE\"" >> parameters_phylok_p.txt
echo "profiles_path=\"\"" >> parameters_phylok_p.txt
echo "profiles_phylo_path=\"\"" >> parameters_phylok_p.txt
echo "profiles_broad_path=\"pfam_database/Pfam-A.hmm\"" >> parameters_phylok_p.txt
echo "swissprot_path=\"swissprot_dtb/swissprot\"" >> parameters_phylok_p.txt
echo "motifs_path=\"motifs.txt\"" >> parameters_phylok_p.txt
echo "output_path=\"/home/ps_data/output/SRR12829170_results\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------Protein family options---------" >> parameters_phylok_p.txt
echo "family_code=\"\"" >> parameters_phylok_p.txt
echo "family_code_phylo=\"\"" >> parameters_phylok_p.txt
echo "db_name=\"\"" >> parameters_phylok_p.txt
echo "db_name_phylo=\"\"" >> parameters_phylok_p.txt
echo "input_protein_names_status=\"\"" >> parameters_phylok_p.txt
echo "input_protein_names=\"\"" >> parameters_phylok_p.txt
echo "input_protein_names_phylo_status=\"\"" >> parameters_phylok_p.txt
echo "input_protein_names_phylo=\"\"" >> parameters_phylok_p.txt
echo "name_threshold=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------General options---------" >> parameters_phylok_p.txt
echo "#---Pipeline---" >> parameters_phylok_p.txt
echo "seek_route=\"\"" >> parameters_phylok_p.txt
echo "paired_end=\"\"" >> parameters_phylok_p.txt
echo "compressed=\"\"" >> parameters_phylok_p.txt
echo "create_nr_db_status=\"\"" >> parameters_phylok_p.txt
echo "prefetch_size=\"\"" >> parameters_phylok_p.txt
echo "adapters_status=\"\"" >> parameters_phylok_p.txt
echo "add_seek_info=\"\"" >> parameters_phylok_p.txt
echo "add_taxonomy_info=\"\"" >> parameters_phylok_p.txt
echo "#---FastQC---" >> parameters_phylok_p.txt
echo "skip_fastqc=\"\"" >> parameters_phylok_p.txt
echo "#---BBDuk---" >> parameters_phylok_p.txt
echo "bbduk_max_ram=\"\"" >> parameters_phylok_p.txt
echo "clear_space=\"\"" >> parameters_phylok_p.txt
echo "#---Megahit---" >> parameters_phylok_p.txt
echo "k_list=\"\"" >> parameters_phylok_p.txt
echo "#---Kraken---" >> parameters_phylok_p.txt
echo "kraken_mode=\"\"" >> parameters_phylok_p.txt
echo "kraken_threshold=\"\"" >> parameters_phylok_p.txt
echo "kraken_memory_mapping=\"\"" >> parameters_phylok_p.txt
echo "#---Binning---" >> parameters_phylok_p.txt
echo "binning_tool=\"\"" >> parameters_phylok_p.txt
echo "binning_max_ram=\"\"" >> parameters_phylok_p.txt
echo "bin_contig_len=\"200\"" >> parameters_phylok_p.txt
echo "bin_kmer=\"\"" >> parameters_phylok_p.txt
echo "comebin_batch_size=\"\"" >> parameters_phylok_p.txt
echo "#---CD-HIT---" >> parameters_phylok_p.txt
echo "cd_hit_t=\"\"" >> parameters_phylok_p.txt
echo "cd_hit_max_mem=\"\"" >> parameters_phylok_p.txt
echo "#---Gene prediction---" >> parameters_phylok_p.txt
echo "gene_encoding=\"\"" >> parameters_phylok_p.txt
echo "genetic_code=\"\"" >> parameters_phylok_p.txt
echo "#---HMMER---" >> parameters_phylok_p.txt
echo "val_type=\"\"" >> parameters_phylok_p.txt
echo "second_dom_search=\"\"" >> parameters_phylok_p.txt
echo "e_value_nodom_thr=\"\"" >> parameters_phylok_p.txt
echo "#---Annotation---" >> parameters_phylok_p.txt
echo "add_type=\"organism bias  number of organisms biosample name\"" >> parameters_phylok_p.txt
echo "add_info=\"synthetic metagenome no  10  SAMN16441671\"" >> parameters_phylok_p.txt
echo "#---Threads---" >> parameters_phylok_p.txt
echo "thread_num=\"4\"" >> parameters_phylok_p.txt
echo "filtering_threads=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------Processes performed after---------" >> parameters_phylok_p.txt
echo "after_preprocessing=\"\"" >> parameters_phylok_p.txt
echo "after_assembly=\"\"" >> parameters_phylok_p.txt
echo "after_gene_pred=\"\"" >> parameters_phylok_p.txt
echo "after_binning=\"\"" >> parameters_phylok_p.txt
echo "after_db=\"\"" >> parameters_phylok_p.txt
echo "after_tm=\"\"" >> parameters_phylok_p.txt
echo "after_ap=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------Processes performed up to---------" >> parameters_phylok_p.txt
echo "up_to_sra=\"\"" >> parameters_phylok_p.txt
echo "up_to_databases=\"\"" >> parameters_phylok_p.txt
echo "up_to_preprocessing_com=\"\"" >> parameters_phylok_p.txt
echo "up_to_preprocessing_uncom=\"\"" >> parameters_phylok_p.txt
echo "up_to_alignment=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------Tool enviroments---------" >> parameters_phylok_p.txt
echo "sra_env=\"\"" >> parameters_phylok_p.txt
echo "fastqc_env=\"\"" >> parameters_phylok_p.txt
echo "bbduk_env=\"\"" >> parameters_phylok_p.txt
echo "megahit_env=\"\"" >> parameters_phylok_p.txt
echo "kraken_env=\"\"" >> parameters_phylok_p.txt
echo "metabinner_env=\"\"" >> parameters_phylok_p.txt
echo "comebin_env=\"\"" >> parameters_phylok_p.txt
echo "cdhit_env=\"\"" >> parameters_phylok_p.txt
echo "genepred_env=\"\"" >> parameters_phylok_p.txt
echo "hmmer_env=\"\"" >> parameters_phylok_p.txt
echo "diamond_env=\"\"" >> parameters_phylok_p.txt
echo "taxonkit_env=\"\"" >> parameters_phylok_p.txt
echo "phobius_env=\"\"" >> parameters_phylok_p.txt
echo "bowtie_env=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt

echo "#---------Tool paths---------" >> parameters_phylok_p.txt
echo "conda_bin=\"\"" >> parameters_phylok_p.txt
echo "conda_sh=\"${CONDA_SH_PATH}\"" >> parameters_phylok_p.txt
echo "prefetch_path=\"\"" >> parameters_phylok_p.txt
echo "vdb_validate_path=\"\"" >> parameters_phylok_p.txt
echo "fastq_dump_path=\"\"" >> parameters_phylok_p.txt
echo "fastqc_path=\"\"" >> parameters_phylok_p.txt
echo "gzip_path=\"\"" >> parameters_phylok_p.txt
echo "cat_path=\"\"" >> parameters_phylok_p.txt
echo "bbduk_path=\"\"" >> parameters_phylok_p.txt
echo "megahit_path=\"\"" >> parameters_phylok_p.txt
echo "kraken_path=\"\"" >> parameters_phylok_p.txt
echo "metabinner_bin_path=\"${MBINPATH}\"" >> parameters_phylok_p.txt
echo "comebin_bin_path=\"${CBINPATH}\"" >> parameters_phylok_p.txt
echo "cd_hit_path=\"\"" >> parameters_phylok_p.txt
echo "fraggenescanrs_path=\"ps_tools/fgsrs/FragGeneScanRs\"" >> parameters_phylok_p.txt
echo "hmmscan_path=\""\" >> parameters_phylok_p.txt
echo "hmmpress_path=\"\"" >> parameters_phylok_p.txt
echo "hmmfetch_path=\"\"" >> parameters_phylok_p.txt
echo "diamond_path=\"\"" >> parameters_phylok_p.txt
echo "taxonkit_path=\"\"" >> parameters_phylok_p.txt
echo "phobius_path=\"\"" >> parameters_phylok_p.txt
echo "bowtie_build_path=\"\"" >> parameters_phylok_p.txt
echo "bowtie_path=\"\"" >> parameters_phylok_p.txt
echo "" >> parameters_phylok_p.txt