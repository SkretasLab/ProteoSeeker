#!/bin/bash
# Seek,Phylo - Metabinner/COMEBin - Contigs

# Sourcing
source ./find_conda.sh
source ./find_mb_cb.sh
source ./find_pd.sh

# Moving to main directory.
cd ..

# Creating the file for the input paths and option and writing in it.
echo "#---------Input and output options---------" > parameters_seek_phylomc_c.txt
echo "input_folder=\"input\"" >> parameters_seek_phylomc_c.txt
echo "sra_code=\"SRR12829170\"" >> parameters_seek_phylomc_c.txt
echo "contigs=\"True\"" >> parameters_seek_phylomc_c.txt
echo "protein_input=\"\"" >> parameters_seek_phylomc_c.txt
echo "adapters_path=\"adapters.fa\"" >> parameters_seek_phylomc_c.txt
echo "protein_db_path=\"${PD_PATH}\"" >> parameters_seek_phylomc_c.txt
echo "kraken_db_path=\"ps_tools/kraken2/kraken_databases/minikraken2_v2_8GB_201904_UPDATE\"" >> parameters_seek_phylomc_c.txt
echo "profiles_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "profiles_phylo_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "profiles_broad_path=\"pfam_database/Pfam-A.hmm\"" >> parameters_seek_phylomc_c.txt
echo "swissprot_path=\"swissprot_dtb/swissprot\"" >> parameters_seek_phylomc_c.txt
echo "motifs_path=\"motifs.txt\"" >> parameters_seek_phylomc_c.txt
echo "output_path=\"/home/ps_data/output/SRR12829170_results\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------Protein family options---------" >> parameters_seek_phylomc_c.txt
echo "family_code=\"1191,3957\"" >> parameters_seek_phylomc_c.txt
echo "family_code_phylo=\"9050,9051,9052,9053,9054,9055,9056,9057,9058,9496\"" >> parameters_seek_phylomc_c.txt
echo "db_name=\"cas_als\"" >> parameters_seek_phylomc_c.txt
echo "db_name_phylo=\"phylo_1\"" >> parameters_seek_phylomc_c.txt
echo "input_protein_names_status=\"\"" >> parameters_seek_phylomc_c.txt
echo "input_protein_names=\"\"" >> parameters_seek_phylomc_c.txt
echo "input_protein_names_phylo_status=\"\"" >> parameters_seek_phylomc_c.txt
echo "input_protein_names_phylo=\"\"" >> parameters_seek_phylomc_c.txt
echo "name_threshold=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------General options---------" >> parameters_seek_phylomc_c.txt
echo "#---Pipeline---" >> parameters_seek_phylomc_c.txt
echo "seek_route=\"\"" >> parameters_seek_phylomc_c.txt
echo "paired_end=\"False\"" >> parameters_seek_phylomc_c.txt
echo "compressed=\"False\"" >> parameters_seek_phylomc_c.txt
echo "create_nr_db_status=\"True\"" >> parameters_seek_phylomc_c.txt
echo "prefetch_size=\"\"" >> parameters_seek_phylomc_c.txt
echo "adapters_status=\"\"" >> parameters_seek_phylomc_c.txt
echo "add_seek_info=\"\"" >> parameters_seek_phylomc_c.txt
echo "add_taxonomy_info=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---FastQC---" >> parameters_seek_phylomc_c.txt
echo "skip_fastqc=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---BBDuk---" >> parameters_seek_phylomc_c.txt
echo "bbduk_max_ram=\"\"" >> parameters_seek_phylomc_c.txt
echo "clear_space=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---Megahit---" >> parameters_seek_phylomc_c.txt
echo "k_list=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---Kraken---" >> parameters_seek_phylomc_c.txt
echo "kraken_mode=\"False\"" >> parameters_seek_phylomc_c.txt
echo "kraken_threshold=\"\"" >> parameters_seek_phylomc_c.txt
echo "kraken_memory_mapping=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---Binning---" >> parameters_seek_phylomc_c.txt
echo "binning_tool=\"\"" >> parameters_seek_phylomc_c.txt
echo "binning_max_ram=\"\"" >> parameters_seek_phylomc_c.txt
echo "bin_contig_len=\"200\"" >> parameters_seek_phylomc_c.txt
echo "bin_kmer=\"\"" >> parameters_seek_phylomc_c.txt
echo "comebin_batch_size=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---CD-HIT---" >> parameters_seek_phylomc_c.txt
echo "cd_hit_t=\"\"" >> parameters_seek_phylomc_c.txt
echo "cd_hit_max_mem=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---Gene prediction---" >> parameters_seek_phylomc_c.txt
echo "gene_encoding=\"\"" >> parameters_seek_phylomc_c.txt
echo "genetic_code=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---HMMER---" >> parameters_seek_phylomc_c.txt
echo "val_type=\"\"" >> parameters_seek_phylomc_c.txt
echo "second_dom_search=\"\"" >> parameters_seek_phylomc_c.txt
echo "e_value_nodom_thr=\"\"" >> parameters_seek_phylomc_c.txt
echo "#---Annotation---" >> parameters_seek_phylomc_c.txt
echo "add_type=\"organism bias  number of organisms biosample name\"" >> parameters_seek_phylomc_c.txt
echo "add_info=\"synthetic metagenome no  10  SAMN16441671\"" >> parameters_seek_phylomc_c.txt
echo "#---Threads---" >> parameters_seek_phylomc_c.txt
echo "thread_num=\"4\"" >> parameters_seek_phylomc_c.txt
echo "filtering_threads=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------Processes performed after---------" >> parameters_seek_phylomc_c.txt
echo "after_preprocessing=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_assembly=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_gene_pred=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_binning=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_db=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_tm=\"\"" >> parameters_seek_phylomc_c.txt
echo "after_ap=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------Processes performed up to---------" >> parameters_seek_phylomc_c.txt
echo "up_to_sra=\"\"" >> parameters_seek_phylomc_c.txt
echo "up_to_databases=\"\"" >> parameters_seek_phylomc_c.txt
echo "up_to_preprocessing_com=\"\"" >> parameters_seek_phylomc_c.txt
echo "up_to_preprocessing_uncom=\"\"" >> parameters_seek_phylomc_c.txt
echo "up_to_alignment=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------Tool enviroments---------" >> parameters_seek_phylomc_c.txt
echo "sra_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "fastqc_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "bbduk_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "megahit_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "kraken_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "metabinner_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "comebin_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "cdhit_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "genepred_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "hmmer_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "diamond_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "taxonkit_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "phobius_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "bowtie_env=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt

echo "#---------Tool paths---------" >> parameters_seek_phylomc_c.txt
echo "conda_bin=\"${ABINPATH}\"" >> parameters_seek_phylomc_c.txt
echo "prefetch_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "vdb_validate_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "fastq_dump_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "fastqc_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "gzip_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "cat_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "bbduk_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "megahit_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "kraken_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "metabinner_bin_path=\"${MBINPATH}\"" >> parameters_seek_phylomc_c.txt
echo "comebin_bin_path=\"${CBINPATH}\"" >> parameters_seek_phylomc_c.txt
echo "cd_hit_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "fraggenescanrs_path=\"ps_tools/fgsrs/FragGeneScanRs\"" >> parameters_seek_phylomc_c.txt
echo "hmmscan_path=\""\" >> parameters_seek_phylomc_c.txt
echo "hmmpress_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "hmmfetch_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "diamond_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "taxonkit_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "phobius_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "bowtie_build_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "bowtie_path=\"\"" >> parameters_seek_phylomc_c.txt
echo "" >> parameters_seek_phylomc_c.txt