import os
import sys
import time
import shutil
import string
import random
# ProteoSeeker modules
sys.path.append('modules')
import sra_process
import path_process
import gene_process
import bbduk_process
import hmmer_process
import motif_process
import family_process
import fastqc_process
import contig_process
import cd_hit_process
import bowtie_process
import adapter_process
import megahit_process
import diamond_process
import phobius_process
import help_message_show
import bin_general_process
import information_collect
import txt_results_process
import tsv_results_process
import supportive_functions
import bin_taxonomy_process
import fastqc_results_process
import kraken_bracken_process
import kraken_binning_process
import domain_coverage_compute
import profile_database_process
import protein_database_process
import predicted_family_process
import comebin_metabinner_process


def proteoseek(input_folder=None, sra_code=False, contigs=False, protein_input=False, adapters_path="adapters.fa", protein_db_path="", kraken_db_path="", profiles_path="", profiles_taxonomy_path="", profiles_broad_path="", swissprot_path="", motifs_path="motifs.txt", options_file_path="", output_path="", seek_mode=True, taxonomy_mode=False, seek_route=1, taxonomy_route=1, seek_family_code=None, taxonomy_family_code=None, seek_db_name="", taxonomy_db_name="", input_seek_protein_names_status=False, input_seek_protein_names=None, input_taxonomy_protein_names_status=False, input_taxonomy_protein_names=None, name_thr=0.5, paired_end=True, compressed=True, prefetch_size=20, adapters_status="pre", add_seek_info=True, add_taxonomy_info=True, skip_fastqc=False, bbduk_max_ram=4, clear_space=False, k_list=None, kraken_threshold="", kraken_memory_mapping=True, bracken_length=150, bracken_level="S", bracken_threshold=10, binning_tool=1, bin_ram_ammount=4, bin_num_contig_len=500, bin_num_kmer=4, comebin_batch_size=256, cd_hit_t=0.99, cd_hit_mem=4000, prs_source=1, genetic_code=11, val_type="--cut_ga ", second_dom_search=True, e_value_nodom_thr=1e-70, add_type="", add_info="", thread_num=4, pdf_threads=None, after_trimming=False, after_assembly=False, after_gene_pred=False, after_binning=False, after_mapping=False, after_db=False, after_tm=False, after_ap=False, up_to_sra=False, up_to_databases=False, up_to_trimming_com=False, up_to_trimming_uncom=False, up_to_alignment=False, sra_env="ps_sra_tools", fastqc_env="ps_fastqc", bbduk_env="ps_bbtools", megahit_env="ps_megahit", kraken_env="ps_kraken", bracken_env="ps_bracken", metabinner_env="ps_metabinner", comebin_env="ps_comebin", cdhit_env="ps_cd_hit", genepred_env="", hmmer_env="ps_hmmer", diamond_env="ps_diamond", taxonkit_env="ps_taxonkit", phobius_env="ps_phobius", bowtie_env="ps_bowtie", conda_bin="", conda_sh="", prefetch_path="", vdb_validate_path="", fastq_dump_path="", fastqc_path="", gzip_path="", cat_path="", bbduk_path="", megahit_path="", kraken_path="", bracken_path="", alpha_diversity_path="", metabinner_bin_path="", comebin_bin_path="", fraggenescanrs_path="", hmmscan_path="", hmmpress_path="", hmmfetch_path="", diamond_path="", cd_hit_path="", taxonkit_path="", phobius_path="", bowtie_build_path="", bowtie_path="", input_command=None):
    # Start time
    start_time = time.time()
    
    # Print a new line in the output.
    print("\n----------------")
    print("Pipeline Start")
    print("----------------")
    print()

    if options_file_path:
        print("\nOptions file path: {}".format(options_file_path))
    else:
        print("\nOptions file path: None")
    
    # Read the values for the options from a file, if its was specified.
    if options_file_path:
        if not os.path.exists(options_file_path):
            print("\nThe path to the options file is wrong. Exiting.")
            exit()
        print("Option values set from the input file:")
        option_lines = supportive_functions.read_file(options_file_path)
        for line in option_lines:
            if line:
                if line[0] == "#":
                    continue
                if "=" in line:
                    line_splited = line.split("=")
                    option_type = line_splited[0]
                    option_value = line_splited[1]
                    if option_value and option_value != "\"\"":
                        print("{}: {}".format(option_type, option_value))
                        # Procesing the option value.
                        if option_value[0] == "\"":
                            option_value = option_value[1:]
                        if option_value[-1] == "\"":
                            option_value = option_value[:-1]
                        # Passing the options value to the corresponding option types.
                        # Input and output options
                        if option_type == "input_folder":
                            input_folder = option_value
                        if option_type == "sra_code":
                            sra_code = option_value
                        if option_type == "contigs":
                            contigs = option_value
                            arg_str = "contigs"
                            contigs = supportive_functions.check_i_value(contigs, arg_str)
                        if option_type == "protein_input":
                            protein_input = option_value
                            arg_str = "protein_input"
                            protein_input = supportive_functions.check_i_value(protein_input, arg_str)
                        if option_type == "adapters_path":
                            adapters_path = option_value
                        if option_type == "protein_db_path":
                            protein_db_path = option_value
                        if option_type == "kraken_db_path":
                            kraken_db_path = option_value
                        if option_type == "profiles_path":
                            profiles_path = option_value
                        if option_type == "profiles_taxonomy_path":
                            profiles_taxonomy_path = option_value
                        if option_type == "profiles_broad_path":
                            profiles_broad_path = option_value
                        if option_type == "swissprot_path":
                            swissprot_path = option_value
                        if option_type == "motifs_path":
                            motifs_path = option_value
                        if option_type == "output_path":
                            output_path = option_value
                            if output_path[-1] == "/":
                                output_path = output_path[:-1]
                        # Modes and Routes
                        if option_type == "seek_mode":
                            seek_mode = option_value
                            arg_str = "seek_mode"
                            seek_mode = supportive_functions.check_i_value(seek_mode, arg_str)
                        if option_type == "taxonomy_mode":
                            taxonomy_mode = option_value
                            arg_str = "taxonomy_mode"
                            taxonomy_mode = supportive_functions.check_i_value(taxonomy_mode, arg_str)
                        if option_type == "seek_route":
                            seek_route = int(option_value)
                        if option_type == "taxonomy_route":
                            taxonomy_route = int(option_value)
                        # Protein family options
                        if option_type == "seek_family_code":
                            seek_family_code = option_value
                        if option_type == "taxonomy_family_code":
                            taxonomy_family_code = option_value
                        if option_type == "seek_db_name":
                            seek_db_name = option_value
                        if option_type == "taxonomy_db_name":
                            taxonomy_db_name = option_value
                        if option_type == "input_seek_protein_names_status":
                            input_seek_protein_names_status = option_value
                            arg_str = "input_seek_protein_names_status"
                            input_seek_protein_names_status = supportive_functions.check_i_value(input_seek_protein_names_status, arg_str)
                        if option_type == "input_seek_protein_names":
                            input_seek_protein_names = option_value
                            input_seek_protein_names = input_seek_protein_names.replace("_", " ")
                            input_seek_protein_names = input_seek_protein_names.split(",")
                        if option_type == "input_taxonomy_protein_names_status":
                            input_taxonomy_protein_names_status = option_value
                            arg_str = "input_taxonomy_protein_names_status"
                            input_taxonomy_protein_names_status = supportive_functions.check_i_value(input_taxonomy_protein_names_status, arg_str)
                        if option_type == "input_taxonomy_protein_names":
                            input_taxonomy_protein_names = option_value
                            input_taxonomy_protein_names = input_taxonomy_protein_names.replace("_", " ")
                            input_taxonomy_protein_names = input_taxonomy_protein_names.split(",")
                        if option_type == "name_threshold":
                            name_thr = option_value
                            if "." in name_thr:
                                name_thr = float(name_thr)
                            else:
                                name_thr = int(name_thr)
                        # General options
                        # Pipeline
                        if option_type == "paired_end":
                            paired_end = option_value
                            arg_str = "paired_end"
                            paired_end = supportive_functions.check_i_value(paired_end, arg_str)
                        if option_type == "compressed":
                            compressed = option_value
                            arg_str = "compressed"
                            compressed = supportive_functions.check_i_value(compressed, arg_str)
                        if option_type == "prefetch_size":
                            prefetch_size = int(option_value)
                        if option_type == "adapters_status":
                            adapters_status = option_value
                        if option_type == "add_seek_info":
                            add_seek_info = option_value
                            arg_str = "-add_seek_info"
                            add_seek_info = supportive_functions.check_i_value(add_seek_info, arg_str)
                        if option_type == "add_taxonomy_info":
                            add_taxonomy_info = option_value
                            arg_str = "-add_taxonomy_info"
                            add_taxonomy_info = supportive_functions.check_i_value(add_taxonomy_info, arg_str)
                        # FastQC
                        if option_type == "skip_fastqc":
                            skip_fastqc = option_value
                            arg_str = "skip_fastqc"
                            skip_fastqc = supportive_functions.check_i_value(skip_fastqc, arg_str)
                        # BBDuk
                        if option_type == "bbduk_max_ram":
                            bbduk_max_ram = int(option_value)
                        if option_type == "clear_space":
                            clear_space = option_value
                            arg_str = "clear_space"
                            clear_space = supportive_functions.check_i_value(clear_space, arg_str)
                        # Megahit
                        if option_type == "k_list":
                            k_list = option_value
                        # Kraken
                        if option_type == "kraken_threshold":
                            kraken_threshold = option_value
                        if option_type == "kraken_memory_mapping":
                            kraken_memory_mapping = option_value
                            arg_str = "kraken_memory_mapping"
                            kraken_memory_mapping = supportive_functions.check_i_value(kraken_memory_mapping, arg_str)
                        if option_type == "bracken_length":
                            bracken_length = int(option_value)
                        if option_type == "bracken_level":
                            bracken_level = option_value
                        if option_type == "bracken_threshold":
                            bracken_threshold = int(option_value)
                        # Binning
                        if option_type == "binning_tool":
                            binning_tool = int(option_value)                        
                        if option_type == "binning_max_ram":
                            bin_ram_ammount = int(option_value)
                        if option_type == "bin_contig_len":
                            bin_num_contig_len = int(option_value)
                        if option_type == "bin_kmer":
                            bin_num_kmer = int(option_value)
                        if option_type == "comebin_batch_size":
                            comebin_batch_size = int(option_value)
                        # CD-HIT
                        if option_type == "cd_hit_max_mem":
                            cd_hit_mem = int(option_value)
                        if option_type == "cd_hit_t":
                            cd_hit_t = float(option_value)
                        # Gene prediction
                        if option_type == "gene_encoding":
                            prs_source = int(option_value)
                        if option_type == "genetic_code":
                            genetic_code = int(option_value)
                        # HMMER
                        if option_type == "val_type":
                            val_type = option_value
                            if val_type == "cut_ga":
                                val_type = "--{} ".format(val_type)
                            elif val_type == "default":
                                val_type = ""
                        if option_type == "second_dom_search":
                            arg_second_dom_search = option_value
                            arg_str = "-second_dom_search"
                            arg_second_dom_search = supportive_functions.check_i_value(arg_second_dom_search, arg_str)
                        if option_type == "e_value_nodom_thr":
                            e_value_nodom_thr = int(option_value)
                            e_value_nodom_thr = float(1*(10**(-e_value_nodom_thr)))                            
                        # Annotation
                        if option_type == "add_type":
                            add_type = option_value
                        if option_type == "add_info":
                            add_info = option_value
                        # Threads
                        if option_type == "thread_num":
                            thread_num = int(option_value)
                        if option_type == "filtering_threads":
                            pdf_threads = int(option_value)
                        # Processes performed after
                        if option_type == "after_preprocessing":
                            after_trimming = option_value
                            arg_str = "-after_preprocessing"
                            after_trimming = supportive_functions.check_i_value(after_trimming, arg_str)
                        if option_type == "after_assembly":
                            after_assembly = option_value
                            arg_str = "after_assembly"
                            after_assembly = supportive_functions.check_i_value(after_assembly, arg_str)
                        if option_type == "after_gene_pred":
                            after_gene_pred = option_value
                            arg_str = "after_gene_pred"
                            after_gene_pred = supportive_functions.check_i_value(after_gene_pred, arg_str)
                        if option_type == "after_binning":
                            after_binning = option_value
                            arg_str = "after_binning"
                            after_binning = supportive_functions.check_i_value(after_binning, arg_str)
                        if option_type == "after_mapping":
                            after_mapping = option_value
                            arg_str = "after_mapping"
                            after_mapping = supportive_functions.check_i_value(after_mapping, arg_str)
                        if option_type == "after_db":
                            after_db = option_value
                            arg_str = "after_db"
                            after_db = supportive_functions.check_i_value(after_db, arg_str)
                        if option_type == "after_tm":
                            after_tm = option_value
                            arg_str = "after_tm"
                            after_tm = supportive_functions.check_i_value(after_tm, arg_str)
                        if option_type == "after_ap":
                            after_ap = option_value
                            arg_str = "after_ap"
                            after_ap = supportive_functions.check_i_value(after_ap, arg_str)
                        # Processes performed up to
                        if option_type == "up_to_sra":
                            up_to_sra = option_value
                            arg_str = "up_to_sra"
                            up_to_sra = supportive_functions.check_i_value(up_to_sra, arg_str)
                        if option_type == "up_to_databases":
                            up_to_databases = option_value
                            arg_str = "up_to_databases"
                            up_to_databases = supportive_functions.check_i_value(up_to_databases, arg_str)
                        if option_type == "up_to_preprocessing_com":
                            up_to_preprocessing_com = option_value
                            arg_str = "up_to_preprocessing_com"
                            up_to_preprocessing_com = supportive_functions.check_i_value(up_to_preprocessing_com, arg_str)
                        if option_type == "up_to_preprocessing_uncom":
                            up_to_preprocessing_uncom = option_value
                            arg_str = "up_to_preprocessing_uncom"
                            up_to_preprocessing_uncom = supportive_functions.check_i_value(up_to_preprocessing_uncom, arg_str)
                        if option_type == "up_to_assembly":
                            up_to_alignment = option_value
                            arg_str = "up_to_assembly"
                            up_to_alignment = supportive_functions.check_i_value(up_to_alignment, arg_str)
                        # Tool environments
                        if option_type == "sra_env":
                            sra_env = option_value
                            if sra_env in ["None", "none"]:
                                sra_env = ""
                        if option_type == "fastqc_env":
                            fastqc_env = option_value
                            if fastqc_env in ["None", "none"]:
                                fastqc_env = ""
                        if option_type == "bbduk_env":
                            bbduk_env = option_value
                            if bbduk_env in ["None", "none"]:
                                bbduk_env = ""
                        if option_type == "megahit_env":
                            megahit_env = option_value
                            if megahit_env in ["None", "none"]:
                                megahit_env = ""
                        if option_type == "kraken_env":
                            kraken_env = option_value
                            if kraken_env in ["None", "none"]:
                                kraken_env = ""
                        if option_type == "bracken_env":
                            bracken_env = option_value
                            if bracken_env in ["None", "none"]:
                                bracken_env = ""
                        if option_type == "metabinner_env":
                            metabinner_env = option_value
                            if metabinner_env in ["None", "none"]:
                                metabinner_env = ""
                        if option_type == "comebin_env":
                            comebin_env = option_value
                            if comebin_env in ["None", "none"]:
                                comebin_env = ""
                        if option_type == "cdhit_env":
                            cdhit_env = option_value
                            if cdhit_env in ["None", "none"]:
                                cdhit_env = ""
                        if option_type == "genepred_env":
                            genepred_env = option_value
                            if genepred_env in ["None", "none"]:
                                genepred_env = ""
                        if option_type == "hmmer_env":
                            hmmer_env = option_value
                            if hmmer_env in ["None", "none"]:
                                hmmer_env = ""
                        if option_type == "diamond_env":
                            diamond_env = option_value
                            if diamond_env in ["None", "none"]:
                                diamond_env = ""
                        if option_type == "taxonkit_env":
                            taxonkit_env = option_value
                            if taxonkit_env in ["None", "none"]:
                                taxonkit_env = ""
                        if option_type == "phobius_env":
                            phobius_env = option_value
                            if phobius_env in ["None", "none"]:
                                phobius_env = ""
                        if option_type == "bowtie_env":
                            bowtie_env = option_value
                            if bowtie_env in ["None", "none"]:
                                bowtie_env = ""
                        # Tool paths
                        if option_type == "conda_bin":
                            conda_bin = option_value
                        if option_type == "conda_sh":
                            conda_sh = option_value
                        if option_type == "prefetch_path":
                            prefetch_path = option_value
                        if option_type == "vdb_validate_path":
                            vdb_validate_path = option_value
                        if option_type == "fastq_dump_path":
                            fastq_dump_path = option_value
                        if option_type == "fastqc_path":
                            fastqc_path = option_value
                        if option_type == "gzip_path":
                            gzip_path = option_value
                        if option_type == "cat_path":
                            cat_path = option_value
                        if option_type == "bbduk_path":
                            bbduk_path = option_value
                        if option_type == "megahit_path":
                            megahit_path = option_value
                        if option_type == "kraken_path":
                            kraken_path = option_value
                        if option_type == "bracken_path":
                            bracken_path = option_value
                        if option_type == "alpha_diversity_path":
                            alpha_diversity_path = option_value
                        if option_type == "metabinner_bin_path":
                            metabinner_bin_path = option_value
                        if option_type == "comebin_bin_path":
                            comebin_bin_path = option_value
                        if option_type == "cd_hit_path":
                            cd_hit_path = option_value
                        if option_type == "fraggenescanrs_path":
                            fraggenescanrs_path = option_value
                        if option_type == "hmmscan_path":
                            hmmscan_path = option_value
                        if option_type == "hmmpress_path":
                            hmmpress_path = option_value
                        if option_type == "hmmfetch_path":
                            hmmfetch_path = option_value
                        if option_type == "diamond_path":
                            diamond_path = option_value
                        if option_type == "taxonkit_path":
                            taxonkit_path = option_value
                        if option_type == "phobius_path":
                            phobius_path = option_value
                        if option_type == "bowtie_build_path":
                            bowtie_build_path = option_value
                        if option_type == "bowtie_path":
                            bowtie_path = option_value

    # Check if at least one mode is enabled. If one mode is disabled the other one must be enabled.
    if (not seek_mode) and (not taxonomy_mode):
        print("Both the seek mode and the taxonomy mode are disabled. Exiting.")
        exit()

    # Determine the filtering threads.
    if pdf_threads is None:
        pdf_threads = thread_num

    # The path to the anaconda bin folder is mandatory.
    if ((not conda_bin) or (conda_bin == "") or (not os.path.exists(conda_bin))) and ((not conda_sh) or (conda_sh == "") or (not os.path.exists(conda_sh))):
        print("\nThe path to the anaconda installation directory and the path to the conda.sh file were not found. At least one must be set. Exiting.")
        exit()
    if conda_sh:
        conda_sh_path = conda_sh
    else:
        conda_sh_path = "{}/etc/profile.d/conda.sh".format(conda_bin)

    # If an output path has not been provided then a random string of 10 characters in length is generated again and again, until
    # the output path does not correspond to an existing folder.
    output_path_suffix = None
    if (not output_path) or output_path == "":
        random_str_length = 10
        output_path_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=random_str_length))
        output_path = "results_{}".format(output_path_suffix)
        while(os.path.exists(output_path)):
            output_path_suffix = ''.join(random.choices(string.ascii_uppercase + string.digits, k=random_str_length))
            output_path = "results_{}".format(output_path_suffix)

    print("\nOutput path: {}".format(output_path))
    if (not os.path.exists(output_path)) and (after_ap or after_tm or after_trimming or after_assembly or after_gene_pred or after_binning or  after_mapping or after_db):
        print("\nThe output folder does not exist and the processes have been set to be continued based on an already existing analysis.")
        exit()

    # Protein input
    if protein_input:
        paired_end = False
        compressed = False

    # Binning tool name
    if binning_tool == 1:
        binning_tool_name = "metabinner"
    else:
        binning_tool_name = "comebin"

    # Input paths and output paths.
    output_path_summaries_errors = "{}/summaries_errors".format(output_path)
    output_path_fastqc = "{}/fastqc_results".format(output_path)
    output_path_trimmed = "{}/trimmed_results".format(output_path)
    output_path_megahit = "{}/megahit_results".format(output_path)
    output_path_conitgs = "{}/contigs".format(output_path)
    if binning_tool == 1:
        output_path_bin = "{}/metabinner_results".format(output_path)
    elif binning_tool == 2:
        output_path_bin = "{}/comebin_results".format(output_path)
    else:
        print("\nWrong selection for the binning tool. Changing to the default selection (1) for the binning tool.")
        binning_tool == 1
        output_path_bin = "{}/metabinner_results".format(output_path)
    output_path_genepred = "{}/gene_results".format(output_path)
    output_path_cdhit = "{}/cd_hit".format(output_path)
    output_path_blastp = "{}/blastp_results".format(output_path)
    output_path_hmmer = "{}/hmmer_results".format(output_path)
    output_path_topology = "{}/topology_results".format(output_path)
    output_path_motifs = "{}/motifs_results".format(output_path)
    output_path_bowtie = "{}/bowtie_results".format(output_path)
    output_path_kraken = "{}/kraken_results".format(output_path)
    output_path_annotation = "{}/annotation_results".format(output_path)
    # SRA file paths
    sra_folder = "sra_files"
    sra_run_folder = "{}/{}_files".format(sra_folder, sra_code)
    sra_file = "{}/{}/{}.sra".format(sra_run_folder, sra_code, sra_code)
    fastq_folder = "{}/{}_fastq".format(sra_run_folder, sra_code)
    fastq_single_end = "{}/{}_fastq_single_end".format(sra_run_folder, sra_code)
    # Enzyme-Profile/Domain file paths
    pr_domains_folder = "profile_protein_dbs"
    fpd_gen_folder = "{}/filtered_protein_dbs".format(pr_domains_folder)
    fam_nums_file_name = "{}/prfamilies_numbered.tsv".format(pr_domains_folder)
    fam_pfam_file_name = "{}/prfamilies_pfamdomains.tsv".format(pr_domains_folder)
    pfam_domains_names = "{}/pfam_accs_names.tsv".format(pr_domains_folder)
    hmmer_profile_lengths_file = "{}/profiles_lengths.tsv".format(pr_domains_folder)
    # Adapters path
    adapters_ap_path = "{}/identified_adapters.fa".format(output_path)
    file_read_info_path = "{}/read_information.txt".format(output_path)
    # Contigs paths
    output_megahit_contigs = "{}/megahit_contigs/final.contigs.fa".format(output_path_megahit)
    output_final_contigs = "{}/contigs.fa".format(output_path_conitgs)
    output_final_contigs_formated = "{}/contigs_formated.fna".format(output_path_conitgs)
    output_final_contigs_formated_name = "contigs_formated.fna"
    bowtie_contigs_basename = "{}/contigs_formated_bwtb".format(output_path_conitgs)
    # Binning paths.
    # Get the path of the current working directory.
    enz_dir = os.getcwd()
    # Binning.
    relpath_contigs = "{}/{}".format(output_path_conitgs, output_final_contigs_formated_name)
    fullpath_contigs_formated = os.path.abspath(relpath_contigs)
    output_fullpath_trimmed = os.path.abspath(output_path_trimmed)
    fullpath_tr_fastq_files = "\"{}/\"*fastq".format(output_fullpath_trimmed)
    full_binning_folder = os.path.abspath(output_path_bin)
    relpath_bin_results_folder = "{}/bins".format(output_path_bin)
    full_bin_results_folder = os.path.abspath(relpath_bin_results_folder)
    full_coverage_folder = "{}/coverages".format(full_binning_folder)
    bin_bam_folder = "{}/bam_results".format(full_binning_folder)
    if binning_tool == 1:
        binning_results_path = "{}/metabinner_res/metabinner_result.tsv".format(full_bin_results_folder)
    else:
        binning_results_path = "{}/comebin_res/comebin_res.tsv".format(full_bin_results_folder)
    full_coverage_profile_path = "{}/coverage_profile_f{}.tsv".format(full_coverage_folder, bin_num_contig_len)
    if taxonomy_route == 1:
        binphylo_path_prsids = "{}/binphylo_prsids.tsv".format(output_path_kraken)
    else:
        binphylo_path_prsids = "{}/binphylo_prsids.tsv".format(output_path_bin)
    binphylo_path_prstax = "{}/binphylo_prstax.tsv".format(output_path_bin)
    binphylo_path_binstax = "{}/binphylo_binstax.tsv".format(output_path_bin)
    binphylo_path_binstax_max = "{}/binphylo_binstax_max.tsv".format(output_path_bin)
    binphylo_path_binstax_names = "{}/binphylo_binstax_names.txt".format(output_path_bin)
    binphylo_path_binstax_max_names = "{}/binphylo_binstax_max_names.txt".format(output_path_bin)
    binphylo_freq_taxids_path = "{}/binphylo_freq_taxids.tsv".format(output_path_bin)
    binphylo_maxfreq_taxids_path = "{}/binphylo_maxfreq_taxids.tsv".format(output_path_bin)
    freq_taxids_path = "{}/freq_taxids.txt".format(output_path_bin)
    maxfreq_taxids_path = "{}/maxfreq_taxids.txt".format(output_path_bin)
    freq_lineage_path = "{}/freq_lineages.tsv".format(output_path_bin)
    maxfreq_lineage_path = "{}/maxfreq_lineages.tsv".format(output_path_bin)
    freq_lineage_form_path = "{}/freq_lineages_form.txt".format(output_path_bin)
    maxfreq_lineage_form_path = "{}/maxfreq_lineages_form.txt".format(output_path_bin)
    family_profile_path = "{}/family_profile_frequencies.tsv".format(output_path_bin)
    bin_summary_info_path = "{}/b_summary_info_{}.tsv".format(output_path_bin, binning_tool_name)
    contig_read_summary_path = "{}/cr_summary_{}.tsv".format(output_path_bin, binning_tool_name)
    binner_bin_info_path = "{}/{}_bin_summary.tsv".format(output_path_bin, binning_tool_name)
    # FragGeneScan paths
    output_path_genepred_base = "{}/GP".format(output_path_genepred)
    output_path_genepred_faa = "{}.faa".format(output_path_genepred_base)
    output_path_genepred_ffn = "{}.ffn".format(output_path_genepred_base)
    output_fgs_protein_formatted_1 = "{}/gp_proteins_formatted.fasta".format(output_path_genepred)
    output_fgs_protein_formatted_2 = "{}/gc_proteins_formatted.fasta".format(output_path_genepred)
    gene_contig_dist_path = "{}/gene_contig_distance.txt".format(output_path_genepred)
    # CD_HIT paths
    cd_hit_results_path = "{}/cd_hit_results".format(output_path_cdhit)
    cd_hit_results_fasta_path = "{}/cd_hit_results_formatted.fasta".format(output_path_cdhit)
    # Hmmer paths
    hmmer_dmbl_results = "{}/hmmer_domtblout.txt".format(output_path_hmmer)
    hmmer_dmbl_results_phylo = "{}/hmmer_domtblout_phylo.txt".format(output_path_hmmer)
    hmmer_simple_results = "{}/hmmer_simple.txt".format(output_path_hmmer)
    hmmer_simple_results_phylo = "{}/hmmer_simple_phylo.txt".format(output_path_hmmer)
    hmmer_enz_domains_all_proteins = "{}/hmmer_enz_domains_all_proteins_info.tsv".format(output_path_hmmer)
    hmmer_enz_domains_all_proteins_phylo = "{}/hmmer_enz_domains_all_proteins_info_phylo.tsv".format(output_path_hmmer)
    file_prs_seq_enz_domains_name = "{}/protein_enz_domains.fasta".format(output_path_hmmer)
    file_prs_seq_enz_domains_phylo_name = "{}/protein_enz_domains_phylo.fasta".format(output_path_hmmer)
    file_prs_seq_no_enzs_name = "{}/proteins_no_enz_domains.fasta".format(output_path_hmmer)
    comb_dmbl_results = "{}/comb_domtblout.txt".format(output_path_hmmer)
    comb_simple_results = "{}/comb_simple.txt".format(output_path_hmmer)
    # File with the information for all the domains found in the combined proteins.
    comb_all_domains_proteins = "{}/comb_all_domains_info.tsv".format(output_path_hmmer)
    # Hmmer Profiles paths
    hmmer_common_lengths_file = "{}/profile_similarity_perc.tsv".format(output_path_hmmer)
    # Motifs paths
    input_motifs_results = "{}/input_motifs_output.tsv".format(output_path_motifs)
    # Blastp paths
    # Blastp - SwissProt
    blastp_results_swissprot_file = "{}/blastp_results_swissprot.tsv".format(output_path_blastp)
    blastp_results_swissprot_phylo_file = "{}/blastp_results_phylo_swissprot.tsv".format(output_path_blastp)
    blastp_info_swissprot_file = "{}/blastp_swissprot_info.tsv".format(output_path_blastp)
    blastp_info_swissprot_phylo_file = "{}/blastp_phylo_swissprot_info.tsv".format(output_path_blastp)
    # Blastp - Domains - Nr
    blastp_doms_nr_file = "{}/blastp_doms_nr.tsv".format(output_path_blastp)
    blastp_doms_nr_phylo_file = "{}/blastp_doms_nr_phylo.tsv".format(output_path_blastp)
    blastp_doms_info_nr_file = "{}/blastp_doms_nr_info.tsv".format(output_path_blastp)
    blastp_doms_info_nr_phylo_file = "{}/blastp_doms_nr_phylo_info.tsv".format(output_path_blastp)
    # Blastp - No domains - Nr
    blastp_results_no_doms_nr_file = "{}/blastp_results_no_doms_nr.tsv".format(output_path_blastp)
    blastp_info_no_doms_nr_file = "{}/blastp_no_doms_nr_info.tsv".format(output_path_blastp)
    blastp_no_doms_below_threshold = "{}/blastp_no_doms_nr_below_threshold.fasta".format(output_path_blastp)
    blastp_info_no_doms_below_threshold = "{}/blastp_no_doms_nr_below_threshold_info.tsv".format(output_path_blastp)
    # Blastp - Combined Info - Nr
    blastp_info_comb_nr_file = "{}/blastp_comb_nr_info.tsv".format(output_path_blastp)
    # Proteins found from the two methods combined in a file. This file is located in the HMMER folder, because the
    # next analysis is based on findind all domains of those proteins.
    proteins_combined_file = "{}/proteins_combined.fasta".format(output_path_hmmer)
    # Gene paths
    gene_info_file = "{}/gene_info.tsv".format(output_path_genepred)
    # The path for the file below must change to "PrFamilies_Length.txt". The information about the family of the Swiss-Prot protein which is the first hit of a putative protein
    # against the Swiss-Prot database with the mean and median length of that family should be drawn for this new file instead of the one below. Then a comparison should be made between
    # this fmaily and the one selected at the beginning. A score should be set for the match or mismatch of these familes.
    pr_fams_path = "{}/prfamilies_length.tsv".format(pr_domains_folder)
    family_info_path = "{}/family_pred_info.tsv".format(output_path_blastp)
    # Bowtie2 analysis
    bowtie_single_unaligned_path = "{}/unaligned_singe_end.fastq".format(output_path_bowtie)
    bowtie_paired_unaligned_con_path = "{}/unaligned_con_paired_end.fastq".format(output_path_bowtie)
    bowtie_stats_path = "{}/bowtie_stats.txt".format(output_path_bowtie)
    mapped_reads_path = "{}/mapped_reads.sam".format(output_path_bowtie)
    # kraken2
    kraken_results_path = "{}/taxon_results.kraken2".format(output_path_kraken)
    kraken_report_path = "{}/report_info.k2report".format(output_path_kraken)
    kraken_species_path = "{}/kraken_species.txt".format(output_path_kraken)
    kraken_species_thr_path = "{}/kraken_species_thr".format(output_path_kraken)
    kraken_reads_path = "{}/kraken_reads.tsv".format(output_path_kraken)
    kraken_taxname_path = "{}/kraken_taxa_names.tsv".format(output_path_kraken)
    binned_ctb_path = "{}/contigs_to_bins.tsv".format(output_path_kraken)
    binned_btc_path = "{}/bins_to_contigs.tsv".format(output_path_kraken)
    binned_taxa_path = "{}/binned_taxa.txt".format(output_path_kraken)
    kraken_bin_info_path = "{}/kraken_bin_summary.txt".format(output_path_kraken)
    bracken_filters_path = "{}/bracken_filters.txt".format(output_path_kraken)
    bracken_output_path = "{}/bracken_output.bracken".format(output_path_kraken)
    bracken_report_path = "{}/bracken_report.breport".format(output_path_kraken)
    # Phobius paths
    phobius_input_file_name = "{}/phobius_input.txt".format(output_path_topology)
    phobius_output_file_name = "{}/phobius_results.txt".format(output_path_topology)
    # Phobius information file path
    topology_info_path = "{}/topology_info.tsv".format(output_path_topology)
    # Annotation paths
    if output_path_suffix is not None:
        annotation_file_txt_name = ("{}/annotation_info_{}.txt".format(output_path_annotation, output_path_suffix))
        annotation_file_tsv_name = ("{}/annotation_info_{}.tsv".format(output_path_annotation, output_path_suffix))
    else:
        annotation_file_txt_name = ("{}/annotation_info.txt".format(output_path_annotation))
        annotation_file_tsv_name = ("{}/annotation_info.tsv".format(output_path_annotation))
    # Bash scripts
    diamond_db_bash_name = "diamond_db.sh"
    sra_bash_script = "{}/sra_run.sh".format(sra_run_folder)
    fastqc_1_bash_script = "{}/fastqc_initial.sh".format(output_path_fastqc)
    fastqc_2_bash_script = "{}/fastqc_trimmed.sh".format(output_path_fastqc)
    bbduk_bash_script = "{}/bbduk.sh".format(output_path_trimmed)
    megahit_bash_script = "{}/megahit.sh".format(output_path_megahit)
    metabinner_bash_script = "{}/metabinner.sh".format(output_path_bin)
    comebin_bash_script = "{}/comebin.sh".format(full_binning_folder, output_path_bin)
    gene_bash_script = "{}/genepred.sh".format(output_path_genepred)
    cdhit_bash_script = "{}/cdhit.sh".format(output_path_cdhit)
    hmmer_spec_bash_script = "{}/hmmer_spec.sh".format(output_path_hmmer)
    hmmer_spec_phylo_bash_script = "{}/hmmer_spec_phylo.sh".format(output_path_hmmer)
    blastp_nodoms_script = "{}/blast_nr_nodoms.sh".format(output_path_blastp)
    hmmer_broad_bash_script = "{}/hmmer_broad.sh".format(output_path_hmmer)
    blast_swissprot_bash_script = "{}/dimaond_blast_swissprot.sh".format(output_path_blastp)
    blast_swissprot_phylo_bash_script = "{}/diamond_blast_swissprot_phylo.sh".format(output_path_blastp)
    blast_nr_bash_script = "{}/diamond_blast_nr_doms.sh".format(output_path_blastp)
    blast_nr_phylo_bash_script = "{}/diamond_blast_nr_phylo.sh".format(output_path_blastp)
    tax_freq_bash_script = "{}/taxonkit_freq_names.sh".format(output_path_bin)
    tax_maxfreq_bash_script = "{}/taxonkit_maxfreq_names.sh".format(output_path_bin)
    taxonkit_freq_line_bash_script = "{}/taxonkit_freq_lineages.sh".format(output_path_bin)
    taxonkit_maxfreq_line_bash_script = "{}/taxonkit_maxfreq_lineages.sh".format(output_path_bin)
    kraken_bash_script = "{}/kraken.sh".format(output_path_kraken)
    bracken_bash_script = "{}/bracken.sh".format(output_path_kraken)
    krakentools_bash_script = "{}/krakentools.sh".format(output_path_kraken)
    phobius_bash_script = "{}/phobius.sh".format(output_path_topology)
    bowtie_bash_script = "{}/bowtie.sh".format(output_path_bowtie)
    # TXT stdout and stderr
    fastqc_stdoe_bt_path = "{}/fastqc_bt_stdoe.txt".format(output_path_fastqc)
    fastqc_stdoe_at_path = "{}/fastqc_at_stdoe.txt".format(output_path_fastqc)
    bbduk_stdoe_path = "{}/bbduk_stdoe.txt".format(output_path_trimmed)
    megahit_stdoe_path = "{}/megahit_stdoe.txt".format(output_path_megahit)
    fraggenescanrs_stdoe_path = "{}/fraggenescanrs_stdoe.txt".format(output_path_genepred)
    cd_hit_stdoe_path = "{}/cd_hit_stdoe.txt".format(output_path_cdhit)
    hmmscan_spec_stdoe_path = "{}/hmmscan_spec_stdoe.txt".format(output_path_hmmer)
    hmmscan_spec_phylo_stdoe_path = "{}/hmmscan_spec_phylo_stdoe.txt".format(output_path_hmmer)
    diamond_blastp_nr1_stdoe_path = "{}/blastp_nr1_stdoe.txt".format(output_path_blastp)
    hmmscan_broad_stdoe_path = "{}/hmmscan_broad_stdoe.txt".format(output_path_hmmer)
    diamond_blastp_swissprot_stdoe_path = "{}/blastp_swissprot_stdoe.txt".format(output_path_blastp)
    diamond_blastp_nr2_stdoe_path = "{}/blastp_nr2_stdoe.txt".format(output_path_blastp)
    diamond_blastp_swissprot_phylo_stdoe_path = "{}/blastp_swissprot_phylo_stdoe.txt".format(output_path_blastp)
    diamond_blastp_nr2_phylo_stdoe_path = "{}/blastp_nr2_phylo_stdoe.txt".format(output_path_blastp)
    bin_gen_coverage_stdoe_path = "{}/gen_coverage_stdoe.txt".format(full_binning_folder)
    bin_gen_kmer_stdoe_path = "{}/gen_kmer_stdoe.txt".format(full_binning_folder)
    bin_filter_tooshort_stdoe_path = "{}/filter_tooshort_stdoe.txt".format(full_binning_folder)
    binning_stdoe_path = "{}/binning_stdoe.txt".format(full_binning_folder)
    taxoknit_freq_stdoe_path = "{}/taxonkit_freq_taxids_stdoe.txt".format(output_path_bin)
    taxoknit_maxfreq_stdoe_path = "{}/taxonkit_maxfreq_taxids_stdoe.txt".format(output_path_bin)
    taxonkit_freq_line_stdoe_path = "{}/taxonkit_freq_lineage_stdoe.txt".format(output_path_bin)
    taxonkit_maxfreq_line_stdoe_path = "{}/taxonkit_maxfreq_lineage_stdoe.txt".format(output_path_bin)
    csvtk_freq_stdoe_path = "{}/csvtk_freq_stdoe.txt".format(output_path_bin)
    csvtk_maxfreq_stdoe_path = "{}/csvtk_maxfreq_stdoe.txt".format(output_path_bin)
    phobius_stde_path = "{}/phobius_stdoe.txt".format(output_path_topology)
    bowtie_build_stdoe_path = "{}/bowtie_build_stdoe.txt".format(output_path_bowtie)
    bowtie_stdoe_path = "{}/bowtie_stdoe.txt".format(output_path_bowtie)
    kraken_stde_path = "{}/kraken_stde.txt".format(output_path_kraken)
    bracken_stde_path = "{}/bracken_stdoe.txt".format(output_path_kraken)
    alpha_diversity_stdoe_path = "{}/alpha_diversity_stdoe.txt".format(output_path_kraken)
    # TXT files for versions
    fastqc_version_bt_path = "{}/fastqc_bt_version.txt".format(output_path_fastqc)
    fastqc_version_at_path = "{}/fastqc_at_version.txt".format(output_path_fastqc)
    bbduk_version_path = "{}/bbduk_version.txt".format(output_path_trimmed)
    megahit_version_path = "{}/megahit_version.txt".format(output_path_megahit)
    fraggenescanrs_version_path = "{}/fraggenescanrs_version.txt".format(output_path_genepred)
    cd_hit_version_path = "{}/cd_hit_version.txt".format(output_path_cdhit)
    hmmscan_spec_version_path = "{}/hmmscan_spec_version.txt".format(output_path_hmmer)
    hmmscan_spec_phylo_version_path = "{}/hmmscan_spec_phylo_version.txt".format(output_path_hmmer)
    hmmscan_broad_version_path = "{}/hmmscan_broad_version.txt".format(output_path_hmmer)
    taxoknit_freq_version_path = "{}/taxonkit_freq_taxids_version.txt".format(output_path_bin)
    taxoknit_maxfreq_version_path = "{}/taxonkit_maxfreq_taxids_version.txt".format(output_path_bin)
    taxonkit_freq_line_version_path = "{}/taxonkit_freq_lineage_version.txt".format(output_path_bin)
    taxonkit_maxfreq_line_version_path = "{}/taxonkit_maxfreq_lineage_version.txt".format(output_path_bin)
    csvtk_freq_version_path = "{}/csvtk_freq_version.txt".format(output_path_bin)
    csvtk_maxfreq_version_path = "{}/csvtk_maxfreq_version.txt".format(output_path_bin)
    phobius_version_path = "{}/phobius_version.txt".format(output_path_topology)
    bowtie_build_version_path = "{}/bowtie_build_version.txt".format(output_path_bowtie)
    bowtie_version_path = "{}/bowtie_version.txt".format(output_path_bowtie)
    kraken_version_path = "{}/kraken_version.txt".format(output_path_kraken)
    bracken_version_path = "{}/bracken_version.txt".format(output_path_kraken)
    alpha_diversity_version_path = "{}/alpha_diversity_version.txt".format(output_path_kraken)
    # Time
    time_analyis_path = "{}/time_analysis.tsv".format(output_path)

    # Initialization of variables.
    pr_names_dict = {}
    dict_contigs_bins = {}
    dict_prs_bins = {}
    prs_with_enz_domains = []
    prs_with_enz_domains_phylo = []
    prs_blast_thr = []
    protein_ids_below_thr = []
    analysis_fam_names = []
    dict_hm = {}
    dict_sp = {}
    dict_nr = {}
    dict_top = {}
    bins_prs_tax_dict = {}
    bin_tax_dict = {}
    bin_tax_max_dict = {}
    bin_group_tax_dict = {}
    pr_tax_info_dcit = {}
    dict_bins_prs_sorted = {}
    dict_input_motifs = {}
    dict_genes = {}
    contig_gene_dist_dict = {}
    swiss_fams_len_comp_dict = {}
    reads_to_contigs_dict = {}
    contigs_to_reads_dict = {}
    kraken_species_dict = {}
    bracken_species_thr_dict = {}
    read_to_species_dict = {}
    taxid_to_species_dict = {}
    tr_ex_file_paths = []
    tr_ex_file_paths_p = {}
    dict_seqs = {}
    time_dict = {
        "tool_time": "None",
        "sra_time": "None",
        "dbs_time": "None",
        "fastqc_initial_time": "None",
        "preprocessing_time": "None",
        "fastqc_final_time": "None",
        "assembly_time": "None",
        "gene_prediction_time": "None",
        "gene_annotation_time": "None",
        "cd_hit_time": "None",
        "kraken_time": "None",
        "kraken_specific_time": "None",
        "binning_time": "None",
        "bowtie_time": "None",
        "hmmer_spec_time": "None",
        "hmmer_spec_taxonomy_time": "None",
        "blastp_fpd_no_doms_time": "None",
        "blastp_fpd_swiss_doms_time": "None",
        "blastp_fpd_swiss_taxonomy_time": "None",
        "bin_analysis_cm_time": "None",
        "bin_taxonomy_cm_time": "None",
        "kraken_binning_time": "None",
        "hmmer_broad": "None",
        "topology_time": "None",
        "motifs_time": "None",
        "family_prediction_time": "None",
        "info_collection_time": "None",
        "results_time": "None"
    }
    if input_seek_protein_names is None:
        input_seek_protein_names = []
    if input_taxonomy_protein_names is None:
        input_taxonomy_protein_names = []

    # Initialize dictionaries for the different kinds of information.
    # Determine whether the output folder will be overwritten.
    if (not after_trimming) and (not after_assembly) and (not after_gene_pred) and (not protein_input) and (not after_binning) and (not after_mapping) and (not after_db) and (not after_tm) and (not after_ap):
        # The output folder is created. If it already exists then it is deleted and recreated.
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        # The output folder of the summary and error files.
        if os.path.exists(output_path_summaries_errors):
            shutil.rmtree(output_path_summaries_errors)
        os.mkdir(output_path_summaries_errors)
    elif protein_input:
        # Create the folder for the results and a folder for the CD-HIT results
        # Copy the input file in the CD-HIT results as the formatted file that would be generated by CD-HIT
        # The output folder is created. If it already exists then it is deleted and recreated.
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        # The output_path_summaries_errors is only needed for trimming, thus not for steps after CD-HIT.
        # Create the CD-HIT folder
        if os.path.exists(output_path_cdhit):
            shutil.rmtree(output_path_cdhit)
        os.mkdir(output_path_cdhit)
    
    # Create the input and output log files in the output folder.
    input_log_file = open("{}/log_input.txt".format(output_path), "a+")
    output_log_file = open("{}/log_output.txt".format(output_path), "a+")
    input_log_file.write("Input command:\n")
    input_log_file.write("{}\n\n".format(input_command))

    # Convert the input kraken thresholds to a list of kraken thresholds.
    if taxonomy_route == 1:
        if "," in kraken_threshold:
            kraken_threshold = kraken_threshold.split(",")
        elif (not kraken_threshold) or (kraken_threshold == ""):
            kraken_threshold = ["0"]
        else:
            kraken_threshold = [kraken_threshold]

    # If no seek family code and no database name have been provided then only taxonomic analysis will be performed.
    # If no taxonomy family code and no database phylo have been provided then only search for a specific protein famlily is performed.
    if (seek_family_code is None) and (not seek_db_name):
        if seek_mode:
            print("Error. The seek mode is enabled and no seek family codes nor a seek database name have been provided. Exiting.")
            exit()
        family_group_name = None
    elif (seek_family_code is not None) and seek_db_name:
        if not seek_mode:
            print("Warning. The seek mode is disabled and seek family codes and a seek database name have been provided.")
        if "," in seek_family_code:
            seek_family_code = seek_family_code.split(",")
        else:
            seek_family_code = [seek_family_code]
        family_group_name = seek_db_name
    elif (seek_family_code is not None):
        if not seek_mode:
            print("Warning. The seek mode is disabled and seek family codes have been provided.")
        if "," in seek_family_code:
            seek_family_code = seek_family_code.split(",")
        else:
            seek_family_code = [seek_family_code]
        family_group_name = "_".join(seek_family_code)
    elif seek_db_name:
        if not seek_mode:
            print("Warning. The seek mode is disabled and a seek family datbase name has been provided.")
        family_group_name = seek_db_name

    # Databases for taxonomic analysis.
    if (taxonomy_family_code is None) and (not taxonomy_db_name):
        if taxonomy_mode and seek_route == 2:
            print("Error. The taxonomy mode is enabled and no taxonomy family codes nor a taxonomy database name have been provided. Exiting.")
            exit()
        family_group_name_phylo = None
    elif (taxonomy_family_code is not None) and taxonomy_db_name:
        if not taxonomy_mode:
            print("Warning. The taxonomy mode is disabled and taxonomy family codes and a taxonomy database name have been provided.")
        else:
            if taxonomy_route == 1:
                print("Warning. The Kraken2 taxonomy route has been selected and taxonomy family codes and a taxonomy database name have been provided.")
        if "," in taxonomy_family_code:
            taxonomy_family_code = taxonomy_family_code.split(",")
        else:
            taxonomy_family_code = [taxonomy_family_code]
        family_group_name_phylo = taxonomy_db_name
    elif (taxonomy_family_code is not None):
        if not taxonomy_mode:
            print("Warning. The taxonomy mode is disabled and taxonomy family codes have been provided.")
        else:
            if taxonomy_route == 1:
                print("Warning. The Kraken2 taxonomy route has been selected and taxonomy family codes have been provided.")
        if "," in taxonomy_family_code:
            taxonomy_family_code = taxonomy_family_code.split(",")
        else:
            taxonomy_family_code = [taxonomy_family_code]
        family_group_name_phylo = "_".join(taxonomy_family_code)
    elif taxonomy_db_name:
        if not taxonomy_mode:
            print("Warning. The taxonomy mode is disabled and a taxonomy family datbase name has been provided.")
        else:
            if taxonomy_route == 1:
                print("Warning. The Kraken2 taxonomy route has been selected and a taxonomy family datbase name has been provided.")
        family_group_name_phylo = taxonomy_db_name

    # SRA file collection, if an SRA code was selected.
    if sra_code:
        start_time_sra = time.time()
        paired_end, input_folder = sra_process.collect_sra(sra_env, sra_folder, prefetch_size, sra_run_folder, sra_file, fastq_folder, fastq_single_end, sra_code, input_log_file, output_log_file, prefetch_path, vdb_validate_path, fastq_dump_path, sra_bash_script, conda_sh_path)
        label = "Process of SRA processing:"
        elpased_time = supportive_functions.end_time_analysis(label, start_time_sra, output_log_file)
        time_dict["sra_time"] = elpased_time
        compressed = True
        contigs = False
    
    # If chosen, the analysis stops here.
    if up_to_sra:
        label = "Process of SRA processing:"
        elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
        time_dict["tool_time"] = elpased_time
        supportive_functions.write_time_dict(time_dict, time_analyis_path)
        exit()

    # Type and info provided to add in the annotation results.
    if add_type:
        if "\t" in add_type:
            add_type = add_type.split("\t")
        else:
            add_type = [add_type]
    if add_info:
        if "\t" in add_info:
            add_info = add_info.split("\t")
        else:
            add_info = [add_info]
    if add_type and add_info:
        if len(add_type) != len(add_info):
            print("\nThe types and information provided to add in the annotation results differ. Exiting.")
            exit()
    elif add_type and (not add_info):
        print("\nType(s) was(were) provided to be added in the annotation results but no information was provided. Exiting.")
        exit()
    elif add_info and (not add_type):
        print("\nInformation was provided to be added in the annotation results but no type(s) was(were) provided. Exiting.")
        exit()

    # The pHMM database.
    profiles_path, analysis_fam_names, family_to_profile_seek_dict = profile_database_process.craete_phmm_db(seek_family_code, family_group_name, pr_domains_folder, fam_nums_file_name, fam_pfam_file_name, pfam_domains_names, hmmer_env, hmmfetch_path, hmmpress_path, profiles_broad_path, input_log_file, output_log_file, conda_sh_path)

    # Determine the seek protein names for the seek mode and seek route 2 or 3.
    if seek_mode and (seek_route == 2 or seek_route == 3):
        fpd_fasta, fpd_name, fpd_folder, pr_names_dict = family_process.find_fam_names(seek_family_code, family_group_name, protein_db_path, fpd_gen_folder, fam_nums_file_name, name_thr, input_seek_protein_names_status, input_seek_protein_names, pr_names_dict, conda_sh_path)
    
    # The taxonomy pHMM database.
    profiles_taxonomy_path, analysis_fam_names_phylo, family_to_profile_phylo_dict = profile_database_process.craete_phmm_db(taxonomy_family_code, family_group_name_phylo, pr_domains_folder, fam_nums_file_name, fam_pfam_file_name, pfam_domains_names, hmmer_env, hmmfetch_path, hmmpress_path, profiles_broad_path, input_log_file, output_log_file, conda_sh_path)

    # Determine the taxonomy protein names for the taxonomy mode and taxonomy route 2.
    if taxonomy_mode and taxonomy_route == 2:
        fpd_fasta_phylo, fpd_name_phylo, pfpd_gen_folder_phylo_name, pr_names_dict = family_process.find_fam_names(taxonomy_family_code, family_group_name_phylo, protein_db_path, fpd_gen_folder, fam_nums_file_name, name_thr, input_taxonomy_protein_names_status, input_taxonomy_protein_names, pr_names_dict, conda_sh_path)
    
    # Create the filtered protein database for the seek mode and seek route 2 or 3 or the taxonomy mode and taxonomy route 2.
    if (seek_mode and (seek_route == 2 or seek_route == 3)) or (taxonomy_mode and taxonomy_route == 2):
        start_time_dbs = time.time()
        protein_database_process.create_pr_db(pr_names_dict, protein_db_path, diamond_db_bash_name, diamond_env, diamond_path, thread_num, pdf_threads, input_log_file, output_log_file, conda_sh_path)    
        # Time for databases
        label = "Process of creating the pHMM and fnr databases (including phylo):"
        elpased_time = supportive_functions.end_time_analysis(label, start_time_dbs, output_log_file)
        time_dict["dbs_time"] = elpased_time

    # If chosen, the analysis stops here.
    if up_to_databases:
        label = "Databases creation:"
        elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
        time_dict["tool_time"] = elpased_time
        supportive_functions.write_time_dict(time_dict, time_analyis_path)
        exit()

    # Determine the input file paths based on whether the input files are paired-end or not, only if no variable has
    # been set to indicate omitting early steps of the pipeline.
    if (not after_gene_pred) and (not after_binning) and (not after_mapping) and (not after_db) and (not after_tm) and (not after_ap):
        file_paths, file_paths_p = path_process.file_paths_creation(input_folder, paired_end)
        # If file paths should be determined but no input files were found, then stop the pipeline.
        if (not file_paths) and (not file_paths_p):
            print("No input files were found.")
            exit()
        # If proteins have been given as input, then copy the input file as the file for CD-HIT.
        if protein_input:
            input_protein_file_path = file_paths[0]
            shutil.copyfile(input_protein_file_path, cd_hit_results_fasta_path)
    else:
        print("\nInput file check skipped.")

    # Folder for the contigs.
    if (not after_assembly) and (not after_gene_pred) and (not after_binning) and (not after_mapping) and (not after_db) and (not after_tm) and (not after_ap):
        # The folder that contains the contigs right untill the point before the process of binning, in the case where binning is implemented. Otherwise, this folder contais the final contigs.
        if os.path.exists(output_path_conitgs):
            shutil.rmtree(output_path_conitgs)
        os.mkdir(output_path_conitgs)
    
    # Quality control
    if (not contigs) and (not after_assembly) and (not after_gene_pred) and (not after_binning) and (not after_mapping) and (not after_db) and (not protein_input) and (not after_tm) and (not after_ap):
        if not after_trimming:
            # FastQC has no specific option for analyzing paired-end reads together. It analyzes them seperatly as it does with each single-end read.
            start_time_fastqc_initial = time.time()
            clear = True
            fastqc_process.fastqc(fastqc_env, output_path_fastqc, file_paths, fastqc_path, clear, thread_num, input_log_file, output_log_file, fastqc_1_bash_script, fastqc_version_bt_path, fastqc_stdoe_bt_path, conda_sh_path)
            label = "Process of initial fastqc:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_fastqc_initial, output_log_file)
            time_dict["fastqc_initial_time"] = elpased_time

            # Process the FastQC results.
            fastqc_results_process.file_reads(output_path_fastqc, paired_end, file_paths, file_paths_p, file_read_info_path)

            # Fix the path to the proper adapters file based on the choice of the user.
            adapter_process.adapter_file(output_path_fastqc, adapters_status, adapters_path, adapters_ap_path)
            
            # BBDuk.
            start_time_trimming = time.time()
            bbduk_process.bbduk(bbduk_env, output_path_trimmed, file_paths, adapters_ap_path, paired_end, file_paths_p, bbduk_path, thread_num, bbduk_max_ram, output_path_summaries_errors, input_log_file, output_log_file, bbduk_bash_script, bbduk_version_path, bbduk_stdoe_path, conda_sh_path)
            label = "Process of trimming:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_trimming, output_log_file)
            time_dict["preprocessing_time"] = elpased_time

            # Delete input files.
            supportive_functions.reduce_volume(clear_space, input_folder, output_path_trimmed, input_log_file, [1])
            # If chosen, the analysis stops here.
            if up_to_trimming_com:
                label = "Up to trimming (1):"
                elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
                time_dict["tool_time"] = elpased_time
                supportive_functions.write_time_dict(time_dict, time_analyis_path)
                exit()
            
            # FastQC has no specific option for analyzing paired-end reads together. It analyzes them seperatly as it does with each single-end read.
            ca_file_paths = []
            pre_enz_file_paths = os.listdir(output_path_trimmed)
            for i in pre_enz_file_paths:
                if "fastq" in i:
                    local_file_path = "{}/{}".format(output_path_trimmed, i)
                    ca_file_paths.append(local_file_path)
            if not skip_fastqc:
                start_time_fastqc_final = time.time()
                clear = False
                fastqc_process.fastqc(fastqc_env, output_path_fastqc, ca_file_paths, fastqc_path, clear, thread_num, input_log_file, output_log_file, fastqc_2_bash_script, fastqc_version_at_path, fastqc_stdoe_at_path, conda_sh_path)
                label = "Process of final fastqc:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_fastqc_final, output_log_file)
                time_dict["fastqc_final_time"] = elpased_time

            # If the original files were compressed then the files which are generated as results by the trimming tool are also compressed.
            # Megahit accepts only non-compressed files. Therefore, if the files were originally compressed they must be uncompressed in order
            # to be used by Megahit.
            if compressed:
                supportive_functions.unzip_files(ca_file_paths, None, False, input_log_file, output_log_file, gzip_path)
        
        # Non-compressed trimmed files are located.
        tr_ex_file_paths_p, tr_ex_file_paths = path_process.locate_nc_files(paired_end, output_path_trimmed)

        fastqc_results_process.file_reads(output_path_fastqc, paired_end, tr_ex_file_paths, tr_ex_file_paths_p, file_read_info_path)

        # Delete compressed trimmed files (if any) before the alignment.
        supportive_functions.reduce_volume(clear_space, input_folder, output_path_trimmed, input_log_file, [2])
        
        # If chosen, the analysis stops here.
        if up_to_trimming_uncom:
            label = "Up to trimming (2):"
            elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
            time_dict["tool_time"] = elpased_time
            supportive_functions.write_time_dict(time_dict, time_analyis_path)
            exit()

        start_time_assembly = time.time()
        megahit_process.megahit(megahit_env, output_path_megahit, tr_ex_file_paths, paired_end, tr_ex_file_paths_p, megahit_path, k_list, thread_num, input_log_file, output_log_file, megahit_bash_script, megahit_version_path, megahit_stdoe_path, conda_sh_path)
        label = "Process of read assembly:"
        elpased_time = supportive_functions.end_time_analysis(label, start_time_assembly, output_log_file)
        time_dict["assembly_time"] = elpased_time

        # The function below leads to the final contigs being saved in the proper folder (regardless of whether contigs were the initial input files or theyr were formed by reads).
        contig_process.contig_formation(contigs, file_paths, output_final_contigs, output_megahit_contigs, output_final_contigs_formated, input_log_file, output_log_file, cat_path)
        
        # If chosen, the analysis stops here.
        if up_to_alignment:
            label = "Up to alignment:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
            time_dict["tool_time"] = elpased_time
            supportive_functions.write_time_dict(time_dict, time_analyis_path)
            exit()

    if (not after_db) and (not protein_input) and (not after_binning) and (not after_mapping) and (not after_tm) and (not after_ap):
        if not after_gene_pred:
            # Gene scanning.
            start_time_genepred = time.time()
            gene_process.gene_prediction(genepred_env, output_path_genepred, output_path_genepred_base, output_path_genepred_faa, output_fgs_protein_formatted_1, fullpath_contigs_formated, fraggenescanrs_path, thread_num, input_log_file, output_log_file, gene_bash_script, fraggenescanrs_version_path, fraggenescanrs_stdoe_path, conda_sh_path)
            label = "Process of gene prediction:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_genepred, output_log_file)
            time_dict["gene_prediction_time"] = elpased_time
            
        # Characterize the gene sequences of the putative proteins.
        if  seek_mode:
            start_time_geneann = time.time()
            gene_sequences_dict, dict_genes, contig_gene_dist_dict = gene_process.gene_annotation(output_path_genepred, output_path_genepred_ffn, output_fgs_protein_formatted_2, output_final_contigs_formated, gene_contig_dist_path, gene_info_file, genetic_code)
            label = "Process of gene annotation:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_geneann, output_log_file)
            time_dict["gene_annotation_time"] = elpased_time

        # Protein redundancy.
        start_time_cd_hit = time.time()
        cd_hit_process.cd_hit(cdhit_env, output_path_cdhit, prs_source, output_fgs_protein_formatted_1, output_fgs_protein_formatted_2, cd_hit_results_path, cd_hit_results_fasta_path, cd_hit_t, cd_hit_path, cd_hit_mem, thread_num, input_log_file, output_log_file, cdhit_bash_script, cd_hit_version_path, cd_hit_stdoe_path, conda_sh_path)
        label = "Process of cd-hit:"
        elpased_time = supportive_functions.end_time_analysis(label, start_time_cd_hit, output_log_file)
        time_dict["cd_hit_time"] = elpased_time

    if (not after_mapping) and (not after_db) and (not protein_input) and (not after_tm) and (not after_ap):
        if not after_binning:
            # Non-compressed trimmed files are located.
            tr_ex_file_paths_p, tr_ex_file_paths = path_process.locate_nc_files(paired_end, output_path_trimmed)

            # Taxonmic analysis from kraken2.
            if taxonomy_mode and (taxonomy_route == 1):
                start_time_kraken = time.time()
                kraken_status = True
                bracken_status = True
                kraken_species_dict, bracken_species_thr_dict, read_to_species_dict, taxid_to_species_dict, time_dict = kraken_bracken_process.kraken(paired_end, tr_ex_file_paths_p, tr_ex_file_paths, output_path_kraken, kraken_db_path, kraken_results_path, kraken_report_path, kraken_threshold, kraken_species_path, kraken_species_thr_path, kraken_reads_path, kraken_taxname_path, conda_sh_path, kraken_env, kraken_path, kraken_bash_script, kraken_stde_path, kraken_version_path, kraken_status, kraken_memory_mapping, bracken_status, bracken_bash_script, bracken_path, bracken_env, bracken_output_path, bracken_report_path, bracken_length, bracken_level, bracken_threshold, bracken_stde_path, bracken_version_path, krakentools_bash_script, alpha_diversity_path, alpha_diversity_version_path, alpha_diversity_stdoe_path, thread_num, time_dict, bracken_filters_path, input_log_file, output_log_file)
                label = "Process of kraken2:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_kraken, output_log_file)
                time_dict["kraken_time"] = elpased_time

            # If selected and if FastQC were given as input then run MetaBinner2
            if taxonomy_mode and (taxonomy_route == 2):
                start_time_binning = time.time()
                comebin_metabinner_process.binning(binning_tool, metabinner_env, comebin_env, contigs, protein_input, output_path_bin, metabinner_bin_path, comebin_bin_path, fullpath_contigs_formated, bin_bam_folder, enz_dir, output_path_conitgs, full_coverage_folder, fullpath_tr_fastq_files, full_bin_results_folder, full_coverage_profile_path, thread_num, input_log_file, output_log_file, metabinner_bash_script, comebin_bash_script, bin_ram_ammount, bin_num_contig_len, bin_num_kmer, bin_gen_coverage_stdoe_path, bin_gen_kmer_stdoe_path, bin_filter_tooshort_stdoe_path, binning_stdoe_path, comebin_batch_size, conda_sh_path)
                label = "Process of binning:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_binning, output_log_file)
                time_dict["binning_time"] = elpased_time
        
        # Bowtie
        if taxonomy_mode:
            # Non-compressed trimmed files are located.
            tr_ex_file_paths_p, tr_ex_file_paths = path_process.locate_nc_files(paired_end, output_path_trimmed)
                    
            # Mapping reads to contigs.
            start_time_bowtie = time.time()
            bowtie_status = True
            reads_to_contigs_dict, contigs_to_reads_dict = bowtie_process.bowtie(output_path_bowtie, paired_end, tr_ex_file_paths, tr_ex_file_paths_p, conda_sh_path, bowtie_build_path, bowtie_build_version_path, bowtie_build_stdoe_path, bowtie_env, bowtie_bash_script, bowtie_path, bowtie_version_path, bowtie_stdoe_path, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, mapped_reads_path, thread_num, output_final_contigs_formated, bowtie_contigs_basename, bowtie_status, input_log_file, output_log_file)
            label = "Process of read mapping:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_bowtie, output_log_file)
            time_dict["bowtie_time"] = elpased_time

    # Analysis modes.
    # 1: Proteins identified by HMMER with any of the specified domains.
    # 2: Proteins with e-score value below the set (or default) threshold from the results of Blastp against the
    # filtered nr adtabase.
    # 3: Proteins combined from both modes of analyis.
    if (not after_tm) and (not after_ap):
        if not after_db:
            if (seek_route == 1 or seek_route == 3):
                if seek_mode:
                    # HMMER Specified for selected protein family.
                    start_time_hmmer_spec = time.time()
                    phylo_analysis = False
                    prs_with_enz_domains = hmmer_process.hmmer_spec(hmmer_env, output_path_hmmer, hmmer_dmbl_results, hmmer_simple_results, hmmer_enz_domains_all_proteins, cd_hit_results_path, hmmscan_path, profiles_path, val_type, file_prs_seq_enz_domains_name, thread_num, input_log_file, output_log_file, hmmer_spec_bash_script, phylo_analysis, hmmscan_spec_version_path, hmmscan_spec_stdoe_path, conda_sh_path)
                    label = "Process of specified HMMER:"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time_hmmer_spec, output_log_file)
                    time_dict["hmmer_spec_time"] = elpased_time

                # HMMER specified for phylogenetic analysis.
                if profiles_taxonomy_path:
                    start_time_hmmer_spec_phylo = time.time()
                    phylo_analysis = True
                    prs_with_enz_domains_phylo = hmmer_process.hmmer_spec(hmmer_env, output_path_hmmer, hmmer_dmbl_results_phylo, hmmer_simple_results_phylo, hmmer_enz_domains_all_proteins_phylo, cd_hit_results_path, hmmscan_path, profiles_taxonomy_path, val_type, file_prs_seq_enz_domains_phylo_name, thread_num, input_log_file, output_log_file, hmmer_spec_phylo_bash_script, phylo_analysis, hmmscan_spec_phylo_version_path, hmmscan_spec_phylo_stdoe_path, conda_sh_path)
                    label = "Process of specified taxonomy HMMER:"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time_hmmer_spec_phylo, output_log_file)
                    time_dict["hmmer_spec_taxonomy_time"] = elpased_time

            if (seek_route == 2 or seek_route == 3) and seek_mode:
                # Blast proteins with no domains from the selected protein family.
                start_time_nr_no_doms = time.time()
                prs_blast_thr, protein_ids_below_thr = diamond_process.diamond_first_round_fpd(diamond_env, output_path_hmmer, output_path_blastp, prs_with_enz_domains, blastp_results_no_doms_nr_file, blastp_info_no_doms_nr_file, blastp_no_doms_below_threshold, blastp_info_no_doms_below_threshold, cd_hit_results_path, file_prs_seq_no_enzs_name, fpd_fasta, fpd_name, e_value_nodom_thr, thread_num, input_log_file, output_log_file, diamond_path, blastp_nodoms_script, diamond_blastp_nr1_stdoe_path, conda_sh_path)
                label = "Process of BLASTP against fnr for proteins without domains:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_nr_no_doms, output_log_file)
                time_dict["blastp_fpd_no_doms_time"] = elpased_time

            if (not taxonomy_mode) or seek_mode:
                # Combine both the proteins from both lists in a file and process to their annotation.
                supportive_functions.combine_predictions(file_prs_seq_enz_domains_name, blastp_no_doms_below_threshold, proteins_combined_file)

                # Counters for the number of proteins with hits against the seek profiles and for the proteins with not hits but with low e-value scores against the fpd.
                pr_doms_num = 0
                prs_nr_num = 0
                prs_nr_thr_num = 0
                if prs_with_enz_domains:
                    pr_doms_num = len(prs_with_enz_domains)
                if prs_blast_thr:
                    prs_nr_num = len(prs_blast_thr)
                if protein_ids_below_thr:
                    prs_nr_thr_num = len(protein_ids_below_thr)
                all_prs_num = pr_doms_num + prs_nr_thr_num
                print("\nNumber of proteins with specified domains: {}.".format(pr_doms_num))
                print("Number of proteins without specified domains: {}.".format(prs_nr_num))
                print("Number of proteins without specified domains and below threshold: {}.".format(prs_nr_thr_num))
                print("Number of all proteins for further procesing: {}.".format(all_prs_num))
                if all_prs_num == 0:
                    print("\nNo proteins were found for annotation.")
                    label = "Up to combining predictions (after the analysis of proteins without domains of interest):"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
                    time_dict["tool_time"] = elpased_time
                    supportive_functions.write_time_dict(time_dict, time_analyis_path)
                    exit()

                # Domain coverages percentages.
                domain_coverage_compute.dom_percs(comb_all_domains_proteins, hmmer_profile_lengths_file, hmmer_common_lengths_file)
                
                # Blastp. The proteins coming from nr-threshold detection should not be blasted again against the nr database.
                # Combined proteins: SwissProt
                start_time_nr_swiss_doms = time.time()
                phylo_analysis = False
                dict_sp = diamond_process.diamond_second_round_swissprot_fpd(diamond_env, output_path_blastp, blastp_results_swissprot_file, blastp_doms_nr_file, blastp_info_swissprot_file, blastp_doms_info_nr_file, file_prs_seq_enz_domains_name, proteins_combined_file, swissprot_path, fpd_fasta, fpd_name, thread_num, input_log_file, output_log_file, diamond_path, blast_swissprot_bash_script, blast_nr_bash_script, phylo_analysis, conda_sh_path, diamond_blastp_swissprot_stdoe_path, diamond_blastp_nr2_stdoe_path)
                label = "Process of BLASTP against fnr and Swiss-Prot for proteins with domains:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_nr_swiss_doms, output_log_file)
                time_dict["blastp_fpd_swiss_doms_time"] = elpased_time

            # Blastp. The proteins of specifies HMMER for the phylogenetic analysis are blasted against the pnr protein database for the phylogenetic analysis.
            if profiles_taxonomy_path:
                start_time_nr_doms_phyo = time.time()
                phylo_analysis = True
                dict_sp_phylo = diamond_process.diamond_second_round_swissprot_fpd(diamond_env, output_path_blastp, blastp_results_swissprot_phylo_file, blastp_doms_nr_phylo_file, blastp_info_swissprot_phylo_file, blastp_doms_info_nr_phylo_file, file_prs_seq_enz_domains_phylo_name, proteins_combined_file, swissprot_path, fpd_fasta_phylo, fpd_name_phylo, thread_num, input_log_file, output_log_file, diamond_path, blast_swissprot_phylo_bash_script, blast_nr_phylo_bash_script, phylo_analysis, conda_sh_path, diamond_blastp_swissprot_phylo_stdoe_path, diamond_blastp_nr2_phylo_stdoe_path)
                label = "Process of BLASTP against fnr and Swiss-Prot for proteins with taxonomy domains:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_nr_doms_phyo, output_log_file)
                time_dict["blastp_fpd_swiss_taxonomy_time"] = elpased_time

            if seek_mode:
                # Combine the results of blastp against the nr database from analysis modes 1 or/and 2.
                dict_nr = diamond_process.combine_fpd_results(blastp_info_no_doms_below_threshold, blastp_doms_info_nr_file, blastp_info_comb_nr_file)

    if not after_ap:
        if not after_tm:
            if taxonomy_mode:
                if taxonomy_route == 2:
                    # Taxonomy of bins
                    start_time_binpals = time.time()
                    dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict = bin_taxonomy_process.bin_tax_analysis(fpd_fasta_phylo, blastp_doms_info_nr_phylo_file, dict_contigs_bins, binning_results_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prstax, binphylo_path_binstax, binphylo_path_binstax_max, binphylo_path_prsids, binphylo_path_binstax_names, binphylo_path_binstax_max_names, binphylo_freq_taxids_path, binphylo_maxfreq_taxids_path, cd_hit_results_fasta_path, taxonkit_env, taxonkit_path, tax_freq_bash_script, tax_maxfreq_bash_script, taxoknit_freq_version_path, taxoknit_freq_stdoe_path, taxoknit_maxfreq_version_path, taxoknit_maxfreq_stdoe_path, thread_num, conda_sh_path, freq_taxids_path, maxfreq_taxids_path, taxonkit_freq_line_bash_script, taxonkit_maxfreq_line_bash_script, taxonkit_freq_line_version_path, taxonkit_freq_line_stdoe_path, taxonkit_maxfreq_line_version_path, taxonkit_maxfreq_line_stdoe_path, freq_lineage_path, maxfreq_lineage_path, freq_lineage_form_path, maxfreq_lineage_form_path, csvtk_freq_version_path, csvtk_freq_stdoe_path, csvtk_maxfreq_version_path, csvtk_maxfreq_stdoe_path, family_to_profile_phylo_dict, hmmer_enz_domains_all_proteins_phylo, family_profile_path, input_log_file, output_log_file)
                    label = "Process of bin analysis:"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time_binpals, output_log_file)
                    time_dict["bin_analysis_cm_time"] = elpased_time

                    # Analysis of the bins
                    start_time_metacome = time.time()
                    bin_general_process.binning_analysis(mapped_reads_path, binning_results_path, file_read_info_path, reads_to_contigs_dict, contigs_to_reads_dict, output_final_contigs_formated, dict_contigs_bins, contig_read_summary_path, binphylo_path_binstax_max, binphylo_maxfreq_taxids_path, maxfreq_lineage_path, bin_summary_info_path, bin_group_tax_dict, binner_bin_info_path)
                    label = "Process of bin taxonomy:"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time_metacome, output_log_file)
                    time_dict["bin_taxonomy_cm_time"] = elpased_time
                    
                if taxonomy_route == 1:
                    start_time_kraken_binning = time.time()
                    bin_group_tax_dict = kraken_binning_process.kraken_binning(taxonomy_route, contigs_to_reads_dict, read_to_species_dict, taxid_to_species_dict, kraken_taxname_path, kraken_reads_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prsids, binned_ctb_path, binned_btc_path, binned_taxa_path, kraken_bin_info_path)
                    label = "Process of binning based on kraken:"
                    elpased_time = supportive_functions.end_time_analysis(label, start_time_kraken_binning, output_log_file)
                    time_dict["kraken_binning_time"] = elpased_time

    if seek_mode:
        if not after_ap:
            if not after_tm:
                # HMMER Broad.
                start_time_hmmer_broad = time.time()
                prs_with_enz_domains, dict_hm = hmmer_process.hmmer_broad(hmmer_env, prs_with_enz_domains, hmmer_dmbl_results, hmmer_simple_results, hmmer_enz_domains_all_proteins, proteins_combined_file, comb_dmbl_results, comb_simple_results, comb_all_domains_proteins, profiles_broad_path, hmmscan_path, val_type, second_dom_search, thread_num, input_log_file, output_log_file, hmmer_broad_bash_script, hmmscan_broad_version_path, hmmscan_broad_stdoe_path, conda_sh_path)
                label = "Process of broad HMMER:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_hmmer_broad, output_log_file)
                time_dict["hmmer_broad"] = elpased_time
                
                # Topology predictions
                start_time_tm = time.time()
                dict_top = phobius_process.phobius(phobius_env, proteins_combined_file, output_path_topology, phobius_input_file_name, phobius_output_file_name, topology_info_path, phobius_path, input_log_file, output_log_file, phobius_bash_script, phobius_version_path, phobius_stde_path, conda_sh_path)
                label = "Process of topology prediction:"
                elpased_time = supportive_functions.end_time_analysis(label, start_time_tm, output_log_file)
                time_dict["topology_time"] = elpased_time

            # Input motifs search.
            start_time_motifs = time.time()
            dict_input_motifs = motif_process.input_motif_search(motifs_path, output_path_motifs, input_motifs_results, proteins_combined_file, input_log_file)
            label = "Process of motif search:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_motifs, output_log_file)
            time_dict["motifs_time"] = elpased_time
            
            # Protein family type prediction and mean length comparison.
            start_time_swiss_fam_pred = time.time()
            swiss_fams_len_comp_dict, dict_seqs = predicted_family_process.pr_fam_pred(analysis_fam_names, blastp_info_swissprot_file, proteins_combined_file, pr_fams_path, family_info_path)
            label = "Process of Swiss-Prot family prediction:"
            elpased_time = supportive_functions.end_time_analysis(label, start_time_swiss_fam_pred, output_log_file)
            time_dict["family_prediction_time"] = elpased_time

    # Collect information. Information is collected and not used directly from the functions above because if processes have been ommited, then it would be possible to gather
    # the information provided by already processes already run in previous analysis from the files generated and present in the folder which will hold the results of an analysis.
    start_time_info_collection = time.time()
    dict_seqs, dict_hm, dict_top, dict_sp, dict_nr, dict_genes, contig_gene_dist_dict, dict_input_motifs, swiss_fams_len_comp_dict, dict_contigs_bins, dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict = information_collect.info_collection(comb_all_domains_proteins, proteins_combined_file, blastp_info_swissprot_file, blastp_info_comb_nr_file, topology_info_path, gene_info_file, gene_contig_dist_path, input_motifs_results, family_info_path, binning_results_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prstax, binphylo_path_binstax, binphylo_path_binstax_max, binphylo_path_prsids, dict_hm, dict_sp, dict_nr, dict_top, swiss_fams_len_comp_dict, dict_input_motifs, dict_genes, contig_gene_dist_dict, dict_contigs_bins, dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict, taxonomy_route, binned_taxa_path, dict_seqs, seek_mode, taxonomy_mode, cd_hit_results_fasta_path, add_seek_info, add_taxonomy_info)
    label = "Process of information collection:"
    elpased_time = supportive_functions.end_time_analysis(label, start_time_info_collection, output_log_file)
    time_dict["info_collection_time"] = elpased_time

    # Convert the taxonomy information to be protein based as the key of a dictionary.
    pr_tax_info_dcit = supportive_functions.tax_info_convert(bin_group_tax_dict)

    # Writing results
    start_time_results = time.time()
    # TXT results.
    txt_results_process.txt_results(annotation_file_txt_name, output_path_annotation, dict_hm, dict_seqs, dict_top, dict_sp, dict_nr, dict_input_motifs, dict_genes, swiss_fams_len_comp_dict, contig_gene_dist_dict, pr_tax_info_dcit, add_type, add_info)
    # EXCEL results.
    tsv_results_process.tsv_results(annotation_file_tsv_name, dict_hm, dict_seqs, dict_top, dict_sp, dict_nr, dict_input_motifs, dict_genes, swiss_fams_len_comp_dict, contig_gene_dist_dict, pr_tax_info_dcit, add_type, add_info)
    label = "Process of writing results:"
    elpased_time = supportive_functions.end_time_analysis(label, start_time_results, output_log_file)
    time_dict["results_time"] = elpased_time

    # End time
    label = "End of analysis:"
    elpased_time = supportive_functions.end_time_analysis(label, start_time, output_log_file)
    time_dict["tool_time"] = elpased_time
    supportive_functions.write_time_dict(time_dict, time_analyis_path)

    # Close log files.
    input_log_file.close()
    output_log_file.close()
    print()
    print("----------------")
    print("Pipeline End")
    print("----------------")
    print()


if __name__ == "__main__":
    # Input and output options
    arg_input_folder = None
    arg_sra_code = False
    arg_prs_source = 1
    arg_adapters_path = "adapters.fa"
    arg_protein_db_path = ""
    arg_kraken_db_path = ""
    arg_profiles_path = ""
    arg_profiles_taxonomy_path = ""
    arg_profiles_broad_path = ""
    arg_swissprot_path = ""
    arg_motifs_path = ""
    arg_output_path = ""
    # Modes and Routes
    arg_seek_mode = True
    arg_taxonomy_mode = False
    arg_seek_route = 1
    arg_taxonomy_route = 1
    # Protein family options
    arg_seek_family_code = None
    arg_taxonomy_family_code = None
    arg_seek_db_name = ""
    arg_taxonomy_db_name = ""
    arg_input_seek_protein_names_status = False
    arg_input_seek_protein_names = None
    arg_input_taxonomy_protein_names_status = False
    arg_input_taxonomy_protein_names = None
    arg_name_thr = 0.5
    # General options
    # Pipeline
    arg_paired_end = True
    arg_compressed = True
    arg_contigs = False
    arg_protein_input = False
    arg_prefetch_size = 20
    arg_adapters_status = "pre"
    arg_add_seek_info = True
    arg_add_taxonomy_info = True
    # FastQC
    arg_skip_fastqc = False
    # BBDuk
    arg_bbduk_max_ram = 4
    arg_clear_space = False
    # Megahit
    arg_k_list = None
    # Kraken
    arg_kraken_threshold = ""
    arg_kraken_memory_mapping = True
    arg_bracken_length = 100
    arg_bracken_level = "S"
    arg_bracken_threshold = 10
    # Binning
    arg_binning_tool = 1
    arg_bin_ram_ammount = 4
    arg_bin_num_contig_len = 500
    arg_bin_num_kmer = 4
    arg_comebin_batch_size = 256
    # CD-HIT
    arg_cd_hit_t = 0.99
    arg_cd_hit_mem = 4000
    # Gene prediction
    arg_genetic_code = 11
    # HMMER
    arg_val_type = "--cut_ga "
    arg_second_dom_search = True
    arg_e_value_nodom_thr = 1e-70
    # Annotation
    arg_add_type = ""
    arg_add_info = ""
    # Threads
    arg_thread_num = 4
    arg_pdf_threads = None
    # Processes performed after
    arg_after_trimming = False
    arg_after_assembly = False
    arg_after_binning = False
    arg_after_mapping = False
    arg_after_gene_pred = False
    arg_after_db = False
    arg_after_tm = False
    arg_after_ap = False
    # Processes performed up to
    arg_up_to_sra = False
    arg_up_to_databases = False
    arg_up_to_trimming_com = False
    arg_up_to_trimming_uncom = False
    arg_up_to_alignment = False
    # Tool environments
    arg_sra_env = "ps_sra_tools"
    arg_fastqc_env = "ps_fastqc"
    arg_bbduk_env = "ps_bbtools"
    arg_megahit_env = "ps_megahit"
    arg_kraken_env = "ps_kraken"
    arg_bracken_env = "ps_bracken"
    arg_metabinner_env = "ps_metabinner"
    arg_comebin_env = "ps_comebin"
    arg_cdhit_env = "ps_cd_hit"
    arg_genepred_env = ""
    arg_hmmer_env = "ps_hmmer"
    arg_diamond_env = "ps_diamond"
    arg_taxonkit_env = "ps_taxonkit"
    arg_phobius_env = "ps_phobius"
    arg_bowtie_env = "ps_bowtie"
    # Tool paths
    arg_conda_bin = ""
    arg_conda_sh = ""
    arg_prefetch_path = ""
    arg_vdb_validate_path = ""
    arg_fastq_dump_path = ""
    arg_fastqc_path = ""
    arg_gzip_path = ""
    arg_cat_path = ""
    arg_bbduk_path = ""
    arg_megahit_path = ""
    arg_kraken_path = ""
    arg_bracken_path = ""
    arg_alpha_diversity_path = ""
    arg_metabinner_bin_path = ""
    arg_comebin_bin_path = ""
    arg_cd_hit_path = ""
    arg_fraggenescanrs_path = ""
    arg_hmmscan_path = ""
    arg_hmmpress_path = ""
    arg_hmmfetch_path = ""
    arg_diamond_path = ""
    arg_taxonkit_path = ""
    arg_phobius_path = ""
    arg_bowtie_build_path = ""
    arg_bowtie_path = ""
    # Options file path
    arg_options_file_path = ""
    # Variable to store the input command
    arg_input_command = ""
    # The input file is the last argument, therefore all following items in the list are parts of the filename.
    # If more than one item is in the list then join them with spaces.
    if len(sys.argv) > 1:
        arg_input_command = "{}".format(sys.argv[0])
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input":
                arg_input_folder = sys.argv[i+1]
                if len(sys.argv[i:]) > 1:
                    arg_input_folder = " ".join(sys.argv[i+1:])
                elif len(sys.argv[i:]) == 0:
                    print("Error. No input file was found. The input file name should be the last argument in the command. Exiting.")
                    exit()
                else:
                    arg_input_folder = sys.argv[i + 1]
            elif sys.argv[i] == "-sc" or sys.argv[i] == "--sra-code":
                arg_sra_code = sys.argv[i+1]
            elif sys.argv[i] == "-c" or sys.argv[i] == "--contigs":
                arg_contigs = sys.argv[i+1]
                arg_str = "-c"
                arg_contigs = supportive_functions.check_i_value(arg_contigs, arg_str)
            elif sys.argv[i] == "-pi" or sys.argv[i] == "--protein-input":
                arg_protein_input = sys.argv[i+1]
                arg_str = "-pi"
                arg_protein_input = supportive_functions.check_i_value(arg_protein_input, arg_str)
            elif sys.argv[i] == "-a" or sys.argv[i] == "--adapters":
                arg_adapters_path = sys.argv[i+1]
            elif sys.argv[i] == "-pdp" or sys.argv[i] == "--protein-database-path":
                arg_protein_db_path = sys.argv[i+1]
            elif sys.argv[i] == "-kdp" or sys.argv[i] == "--kraken-database-path":
                arg_kraken_db_path = sys.argv[i+1]
            elif sys.argv[i] == "-psp" or sys.argv[i] == "--profiles-seek-path":
                arg_profiles_path = sys.argv[i+1]
            elif sys.argv[i] == "-ptp" or sys.argv[i] == "--profiles-taxonomy-path":
                arg_profiles_taxonomy_path = sys.argv[i+1]
            elif sys.argv[i] == "-pbp" or sys.argv[i] == "--profiles-broad-path":
                arg_profiles_broad_path = sys.argv[i+1]
            elif sys.argv[i] == "-sp" or sys.argv[i] == "--swissprot-path":
                arg_swissprot_path = sys.argv[i+1]
            elif sys.argv[i] == "-mop" or sys.argv[i] == "--motifs-path":
                arg_motifs_path = sys.argv[i+1]
            elif sys.argv[i] == "-pfp" or sys.argv[i] == "--parameters-file-path":
                arg_options_file_path = sys.argv[i+1]
            elif sys.argv[i] == "-o" or sys.argv[i] == "--output":
                arg_output_path = sys.argv[i+1]
                if arg_output_path[-1] == "/":
                    arg_output_path = arg_output_path[:-1]
            elif sys.argv[i] == "-sm" or sys.argv[i] == "--seek-mode":
                arg_seek_mode = sys.argv[i+1]
                arg_str = "-sm"
                arg_seek_mode = supportive_functions.check_i_value(arg_seek_mode, arg_str)
            elif sys.argv[i] == "-tm" or sys.argv[i] == "--taxonomy-mode":
                arg_taxonomy_mode = sys.argv[i+1]
                arg_str = "-tm"
                arg_taxonomy_mode = supportive_functions.check_i_value(arg_taxonomy_mode, arg_str)
            elif sys.argv[i] == "-sr" or sys.argv[i] == "--seek-route":
                arg_seek_route = int(sys.argv[i+1])
            elif sys.argv[i] == "-tr" or sys.argv[i] == "--taxonomy-route":
                arg_taxonomy_route = int(sys.argv[i+1])
            elif sys.argv[i] == "-sfc" or sys.argv[i] == "--family-code":
                arg_seek_family_code = sys.argv[i+1]
            elif sys.argv[i] == "-tfc" or sys.argv[i] == "--family-code-taxonomy":
                arg_taxonomy_family_code = sys.argv[i+1]
            elif sys.argv[i] == "-sdn" or sys.argv[i] == "--seek-database-name":
                arg_seek_db_name = sys.argv[i+1]
            elif sys.argv[i] == "-tdn" or sys.argv[i] == "--taxonomy-database-name":
                arg_taxonomy_db_name = sys.argv[i+1]
            elif sys.argv[i] == "-sns" or sys.argv[i] == "--seek-names-status":
                arg_input_seek_protein_names_status = sys.argv[i+1]
                if arg_input_seek_protein_names_status == "False":
                    arg_input_seek_protein_names_status = False
            elif sys.argv[i] == "-spn" or sys.argv[i] == "--seek-protein-names":
                arg_input_seek_protein_names = sys.argv[i+1]
                arg_input_seek_protein_names = arg_input_seek_protein_names.replace("_", " ")
                arg_input_seek_protein_names = arg_input_seek_protein_names.split(",")
            elif sys.argv[i] == "-tns" or sys.argv[i] == "--taxonomy-names-status":
                arg_input_taxonomy_protein_names_status = sys.argv[i+1]
                if arg_input_taxonomy_protein_names_status == "False":
                    arg_input_taxonomy_protein_names_status = False
            elif sys.argv[i] == "-tpn" or sys.argv[i] == "--taxonomy-protein-names":
                arg_input_taxonomy_protein_names = sys.argv[i+1]
                arg_input_taxonomy_protein_names = arg_input_taxonomy_protein_names.replace("_", " ")
                arg_input_taxonomy_protein_names = arg_input_taxonomy_protein_names.split(",")
            elif sys.argv[i] == "-nt" or sys.argv[i] == "--name-threshold":
                arg_name_thr = sys.argv[i+1]
                if "." in arg_name_thr:
                    arg_name_thr = float(arg_name_thr)
                else:
                    arg_name_thr = int(arg_name_thr)
            elif sys.argv[i] == "-p" or sys.argv[i] == "--paired-end":
                arg_paired_end = sys.argv[i+1]
                arg_str = "-p"
                arg_paired_end = supportive_functions.check_i_value(arg_paired_end, arg_str)
            elif sys.argv[i] == "-k" or sys.argv[i] == "--compressed":
                arg_compressed = sys.argv[i+1]
                arg_str = "-k"
                arg_compressed = supportive_functions.check_i_value(arg_compressed, arg_str)
            elif sys.argv[i] == "-ps" or sys.argv[i] == "--prefetch-size":
                arg_prefetch_size = int(sys.argv[i+1])
            elif sys.argv[i] == "-as" or sys.argv[i] == "--adapters-status":
                arg_adapters_status = sys.argv[i+1]                
            elif sys.argv[i] == "-asi" or sys.argv[i] == "--add-seek-info":
                arg_add_seek_info = sys.argv[i+1]
                arg_str = "-osi"
                arg_add_seek_info = supportive_functions.check_i_value(arg_add_seek_info, arg_str)
            elif sys.argv[i] == "-ati" or sys.argv[i] == "--add-taxonomy-info":
                arg_add_taxonomy_info = sys.argv[i+1]
                arg_str = "-ati"
                arg_add_taxonomy_info = supportive_functions.check_i_value(arg_add_taxonomy_info, arg_str)
            elif sys.argv[i] == "-sf" or sys.argv[i] == "--skip-fastqc":
                arg_skip_fastqc = sys.argv[i+1]
                arg_str = "-sf"
                arg_skip_fastqc = supportive_functions.check_i_value(arg_skip_fastqc, arg_str)
            elif sys.argv[i] == "-umr" or sys.argv[i] == "--bbduk-max-ram":
                arg_bbduk_max_ram = int(sys.argv[i+1])
            elif sys.argv[i] == "-cs" or sys.argv[i] == "--clear-space":
                arg_clear_space = sys.argv[i+1]
                arg_str = "-cs"
                arg_clear_space = supportive_functions.check_i_value(arg_clear_space, arg_str)
            elif sys.argv[i] == "-kl" or sys.argv[i] == "--k-list":
                arg_k_list = sys.argv[i+1]
            elif sys.argv[i] == "-kt" or sys.argv[i] == "--kraken-threshold":
                arg_kraken_threshold = sys.argv[i+1]
            elif sys.argv[i] == "-kmm" or sys.argv[i] == "--kraken-memory-mapping":
                arg_kraken_memory_mapping = sys.argv[i+1]
                arg_str = "-kmm"
                arg_kraken_memory_mapping = supportive_functions.check_i_value(arg_kraken_memory_mapping, arg_str)
            elif sys.argv[i] == "-bl" or sys.argv[i] == "--bracken-length":
                arg_bracken_length = int(sys.argv[i+1])
            elif sys.argv[i] == "-bv" or sys.argv[i] == "--bracken-level":
                arg_bracken_level = sys.argv[i+1]
            elif sys.argv[i] == "-bh" or sys.argv[i] == "--bracken-threshold":
                arg_bracken_threshold = int(sys.argv[i+1])
            elif sys.argv[i] == "-bt" or sys.argv[i] == "--binning-tool":
                arg_binning_tool = int(sys.argv[i+1])
            elif sys.argv[i] == "-bmr" or sys.argv[i] == "--binning-max-ram":
                arg_bin_ram_ammount = int(sys.argv[i+1])
            elif sys.argv[i] == "-bc" or sys.argv[i] == "--bin-contig-len":
                arg_bin_num_contig_len = int(sys.argv[i+1])
            elif sys.argv[i] == "-bk" or sys.argv[i] == "--bin-kmer":
                arg_bin_num_kmer = int(sys.argv[i+1])
            elif sys.argv[i] == "-cbs" or sys.argv[i] == "--comebin-batch-size":
                arg_comebin_batch_size = int(sys.argv[i+1])
            elif sys.argv[i] == "-ct" or sys.argv[i] == "--cdhit-threshold":
                arg_cd_hit_t = float(sys.argv[i+1])
            elif sys.argv[i] == "-cmr" or sys.argv[i] == "--cd-hit-max-ram":
                arg_cd_hit_mem = int(sys.argv[i+1])
            elif sys.argv[i] == "-ge" or sys.argv[i] == "--gene-encoding":
                arg_prs_source = int(sys.argv[i+1])
            elif sys.argv[i] == "-gc" or sys.argv[i] == "--genetic-code":
                arg_genetic_code = int(sys.argv[i+1])
            elif sys.argv[i] == "-st" or sys.argv[i] == "--score-type":
                arg_val_type = sys.argv[i+1]
                if arg_val_type == "cut_ga":
                    arg_val_type = "--{} ".format(arg_val_type)
                elif arg_val_type == "default":
                    arg_val_type = ""
            elif sys.argv[i] == "-sds" or sys.argv[i] == "--second-domain-search":
                arg_second_dom_search = sys.argv[i+1]
                arg_str = "-sds"
                arg_second_dom_search = supportive_functions.check_i_value(arg_second_dom_search, arg_str)
            elif sys.argv[i] == "-ndt" or sys.argv[i] == "--no-domains-thr":
                arg_e_value_nodom_thr = int(sys.argv[i+1])
                arg_e_value_nodom_thr = float(1*(10**(-arg_e_value_nodom_thr)))
            elif sys.argv[i] == "-at" or sys.argv[i] == "--add-type":
                arg_add_type = sys.argv[i+1]
            elif sys.argv[i] == "-ai" or sys.argv[i] == "--add-info":
                arg_add_info = sys.argv[i+1]
            elif sys.argv[i] == "-t" or sys.argv[i] == "--threads":
                arg_thread_num = int(sys.argv[i+1])
            elif sys.argv[i] == "-ft" or sys.argv[i] == "--filtering-threads":
                arg_pdf_threads = int(sys.argv[i+1])
            elif sys.argv[i] == "-afp" or sys.argv[i] == "--after-preprocessing":
                arg_after_trimming = sys.argv[i+1]
                arg_str = "-afp"
                arg_after_trimming = supportive_functions.check_i_value(arg_after_trimming, arg_str)
            elif sys.argv[i] == "-afa" or sys.argv[i] == "--after-assembly":
                arg_after_assembly = sys.argv[i + 1]
                arg_str = "-afa"
                arg_after_assembly = supportive_functions.check_i_value(arg_after_assembly, arg_str)
            elif sys.argv[i] == "-afg" or sys.argv[i] == "--after-gene-pred":
                arg_after_gene_pred = sys.argv[i+1]
                arg_str = "-afg"
                arg_after_gene_pred = supportive_functions.check_i_value(arg_after_gene_pred, arg_str)
            elif sys.argv[i] == "-afb" or sys.argv[i] == "--after-binning":
                arg_after_binning = sys.argv[i+1]
                arg_str = "-afb"
                arg_after_binning = supportive_functions.check_i_value(arg_after_binning, arg_str)
            elif sys.argv[i] == "-afm" or sys.argv[i] == "--after-mapping":
                arg_after_mapping = sys.argv[i+1]
                arg_str = "-afb"
                arg_after_mapping = supportive_functions.check_i_value(arg_after_mapping, arg_str)
            elif sys.argv[i] == "-adb" or sys.argv[i] == "--after-db":
                arg_after_db = sys.argv[i+1]
                arg_str = "-adb"
                arg_after_db = supportive_functions.check_i_value(arg_after_db, arg_str)
            elif sys.argv[i] == "-atp" or sys.argv[i] == "--after-topology-prediction":
                arg_after_tm = sys.argv[i+1]
                arg_str = "-atp"
                arg_after_tm = supportive_functions.check_i_value(arg_after_tm, arg_str)
            elif sys.argv[i] == "-afr" or sys.argv[i] == "--after-analysis-processes":
                arg_after_ap = sys.argv[i+1]
                arg_str = "-afr"
                arg_after_ap = supportive_functions.check_i_value(arg_after_ap, arg_str)
            elif sys.argv[i] == "-uts" or sys.argv[i] == "--up-to-sra":
                arg_up_to_sra = sys.argv[i+1]
                arg_str = "-uts"
                arg_up_to_sra = supportive_functions.check_i_value(arg_up_to_sra, arg_str)
            elif sys.argv[i] == "-utd" or sys.argv[i] == "--up-to-databases":
                arg_up_to_databases = sys.argv[i+1]
                arg_str = "-utd"
                arg_up_to_databases = supportive_functions.check_i_value(arg_up_to_databases, arg_str)
            elif sys.argv[i] == "-utpc" or sys.argv[i] == "--up-to-preprocessing-com":
                arg_up_to_trimming_com = sys.argv[i+1]
                arg_str = "-uttc"
                arg_up_to_trimming_com = supportive_functions.check_i_value(arg_up_to_trimming_com, arg_str)
            elif sys.argv[i] == "-utpu" or sys.argv[i] == "--up-to-preprocessing-uncom":
                arg_up_to_trimming_uncom = sys.argv[i+1]
                arg_str = "-uttu"
                arg_up_to_trimming_uncom = supportive_functions.check_i_value(arg_up_to_trimming_uncom, arg_str)
            elif sys.argv[i] == "-uta" or sys.argv[i] == "--up-to-assembly":
                arg_up_to_alignment = sys.argv[i+1]
                arg_str = "-uta"
                arg_up_to_alignment = supportive_functions.check_i_value(arg_up_to_alignment, arg_str)
            elif sys.argv[i] == "-sen" or sys.argv[i] == "--sra-env":
                arg_sra_env = sys.argv[i+1]
                if arg_sra_env in ["None", "none"]:
                    arg_sra_env = ""
            elif sys.argv[i] == "-fen" or sys.argv[i] == "--fastqc-env":
                arg_fastqc_env = sys.argv[i+1]
                if arg_fastqc_env in ["None", "none"]:
                    arg_fastqc_env = ""
            elif sys.argv[i] == "-uen" or sys.argv[i] == "--bbtools-env":
                arg_bbduk_env = sys.argv[i+1]
                if arg_bbduk_env in ["None", "none"]:
                    arg_bbduk_env = ""
            elif sys.argv[i] == "-men" or sys.argv[i] == "--megahit-env":
                arg_megahit_env = sys.argv[i+1]
                if arg_megahit_env in ["None", "none"]:
                    arg_megahit_env = ""
            elif sys.argv[i] == "-ken" or sys.argv[i] == "--kraken-env":
                arg_kraken_env = sys.argv[i+1]
                if arg_kraken_env in ["None", "none"]:
                    arg_kraken_env = ""
            elif sys.argv[i] == "-ren" or sys.argv[i] == "--bracken-env":
                arg_bracken_env = sys.argv[i+1]
                if arg_bracken_env in ["None", "none"]:
                    arg_bracken_env = ""
            elif sys.argv[i] == "-nen" or sys.argv[i] == "--metabinner-env":
                arg_metabinner_env = sys.argv[i+1]
                if arg_metabinner_env in ["None", "none"]:
                    arg_metabinner_env = ""
            elif sys.argv[i] == "-cen" or sys.argv[i] == "--comebin-env":
                arg_comebin_env = sys.argv[i+1]
                if arg_comebin_env in ["None", "none"]:
                    arg_comebin_env = ""
            elif sys.argv[i] == "-ien" or sys.argv[i] == "--cdhit-env":
                arg_cdhit_env = sys.argv[i+1]
                if arg_cdhit_env in ["None", "none"]:
                    arg_cdhit_env = ""
            elif sys.argv[i] == "-gen" or sys.argv[i] == "--genepred-env":
                arg_genepred_env = sys.argv[i+1]
                if arg_genepred_env in ["None", "none"]:
                    arg_genepred_env = ""
            elif sys.argv[i] == "-hen" or sys.argv[i] == "--hmmer-env":
                arg_hmmer_env = sys.argv[i+1]
                if arg_hmmer_env in ["None", "none"]:
                    arg_hmmer_env = ""
            elif sys.argv[i] == "-den" or sys.argv[i] == "--dimaond-env":
                arg_diamond_env = sys.argv[i+1]
                if arg_diamond_env in ["None", "none"]:
                    arg_diamond_env = ""
            elif sys.argv[i] == "-ten" or sys.argv[i] == "--taxonkit-env":
                arg_taxonkit_env = sys.argv[i+1]
                if arg_taxonkit_env in ["None", "none"]:
                    arg_taxonkit_env = ""
            elif sys.argv[i] == "-pen" or sys.argv[i] == "--phobius-env":
                arg_phobius_env = sys.argv[i+1]
                if arg_phobius_env in ["None", "none"]:
                    arg_phobius_env = ""
            elif sys.argv[i] == "-ben" or sys.argv[i] == "--bowtie-env":
                arg_bowtie_env = sys.argv[i+1]
                if arg_bowtie_env in ["None", "none"]:
                    arg_bowtie_env = ""
            elif sys.argv[i] == "-adp" or sys.argv[i] == "--anaconda-dir-path":
                arg_conda_bin = sys.argv[i+1]
            elif sys.argv[i] == "-asp" or sys.argv[i] == "--anaconda-sh-path":
                arg_conda_sh = sys.argv[i+1]
            elif sys.argv[i] == "-rfp" or sys.argv[i] == "--prefetch-path":
                arg_prefetch_path = sys.argv[i+1]
            elif sys.argv[i] == "-vvp" or sys.argv[i] == "--vdb_validate-path":
                arg_vdb_validate_path = sys.argv[i+1]
            elif sys.argv[i] == "-fdp" or sys.argv[i] == "--fastq-dump-path":
                arg_fastq_dump_path = sys.argv[i+1]
            elif sys.argv[i] == "-fp" or sys.argv[i] == "--fastqc-path":
                arg_fastqc_path = sys.argv[i+1]
            elif sys.argv[i] == "-gzp" or sys.argv[i] == "--gzip-path":
                arg_gzip_path = sys.argv[i+1]
            elif sys.argv[i] == "-ctp" or sys.argv[i] == "--cat-path":
                arg_cat_path = sys.argv[i+1]
            elif sys.argv[i] == "-bdp" or sys.argv[i] == "--bbduk-path":
                arg_bbduk_path = sys.argv[i+1]
            elif sys.argv[i] == "-mp" or sys.argv[i] == "--megahit-path":
                arg_megahit_path = sys.argv[i+1]
            elif sys.argv[i] == "-kp" or sys.argv[i] == "--kraken-path":
                arg_kraken_path = sys.argv[i+1]
            elif sys.argv[i] == "-bp" or sys.argv[i] == "--bracken-path":
                arg_bracken_path = sys.argv[i+1]
            elif sys.argv[i] == "-ap" or sys.argv[i] == "--alpha-diversity-path":
                arg_alpha_diversity_path = sys.argv[i+1]
            elif sys.argv[i] == "-bfp" or sys.argv[i] == "--binner-folder-path":
                arg_metabinner_bin_path = sys.argv[i+1]
            elif sys.argv[i] == "-cfp" or sys.argv[i] == "--comebin-folder-path":
                arg_comebin_bin_path = sys.argv[i+1]
            elif sys.argv[i] == "-chp" or sys.argv[i] == "--cd-hit-path":
                arg_cd_hit_path = sys.argv[i+1]
            elif sys.argv[i] == "-fgp" or sys.argv[i] == "--fraggenescars-path":
                arg_fraggenescanrs_path = sys.argv[i+1]
            elif sys.argv[i] == "-hp" or sys.argv[i] == "--hmmscan-path":
                arg_hmmscan_path = sys.argv[i+1]
            elif sys.argv[i] == "-hpp" or sys.argv[i] == "--hmmpress-path":
                arg_hmmpress_path = sys.argv[i+1]
            elif sys.argv[i] == "-hfp" or sys.argv[i] == "--hmmfetch-path":
                arg_hmmfetch_path = sys.argv[i+1]
            elif sys.argv[i] == "-dp" or sys.argv[i] == "--diamond-path":
                arg_diamond_path = sys.argv[i+1]
            elif sys.argv[i] == "-tkp" or sys.argv[i] == "--taxonkit-path":
                arg_taxonkit_path = sys.argv[i+1]
            elif sys.argv[i] == "-php" or sys.argv[i] == "--phobius-folder-path":
                arg_phobius_path = sys.argv[i+1]
            elif sys.argv[i] == "-bbp" or sys.argv[i] == "--bowtie-build-path":
                arg_bowtie_build_path = sys.argv[i+1]
            elif sys.argv[i] == "-bwp" or sys.argv[i] == "--bowtie-path":
                arg_bowtie_path = sys.argv[i+1]
            elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
                help_message_show.help_message()
                exit()
            arg_input_command = "{} {} {}".format(arg_input_command, sys.argv[i], sys.argv[i+1])
    proteoseek(arg_input_folder, arg_sra_code, arg_contigs, arg_protein_input, arg_adapters_path, arg_protein_db_path, arg_kraken_db_path, arg_profiles_path, arg_profiles_taxonomy_path, arg_profiles_broad_path, arg_swissprot_path, arg_motifs_path, arg_options_file_path, arg_output_path, arg_seek_mode, arg_taxonomy_mode, arg_seek_route, arg_taxonomy_route, arg_seek_family_code, arg_taxonomy_family_code, arg_seek_db_name, arg_taxonomy_db_name, arg_input_seek_protein_names_status, arg_input_seek_protein_names, arg_input_taxonomy_protein_names_status, arg_input_taxonomy_protein_names, arg_name_thr, arg_paired_end, arg_compressed, arg_prefetch_size, arg_adapters_status, arg_add_seek_info, arg_add_taxonomy_info, arg_skip_fastqc, arg_bbduk_max_ram, arg_clear_space, arg_k_list, arg_kraken_threshold, arg_kraken_memory_mapping, arg_bracken_length, arg_bracken_level, arg_bracken_threshold, arg_binning_tool, arg_bin_ram_ammount, arg_bin_num_contig_len, arg_bin_num_kmer, arg_comebin_batch_size, arg_cd_hit_t, arg_cd_hit_mem, arg_prs_source, arg_genetic_code, arg_val_type, arg_second_dom_search, arg_e_value_nodom_thr, arg_add_type, arg_add_info, arg_thread_num, arg_pdf_threads, arg_after_trimming, arg_after_assembly, arg_after_gene_pred, arg_after_binning, arg_after_mapping, arg_after_db, arg_after_tm, arg_after_ap, arg_up_to_sra, arg_up_to_databases, arg_up_to_trimming_com, arg_up_to_trimming_uncom, arg_up_to_alignment, arg_sra_env, arg_fastqc_env, arg_bbduk_env, arg_megahit_env, arg_kraken_env, arg_bracken_env, arg_metabinner_env, arg_comebin_env, arg_cdhit_env, arg_genepred_env, arg_hmmer_env, arg_diamond_env, arg_taxonkit_env, arg_phobius_env, arg_bowtie_env, arg_conda_bin, arg_conda_sh, arg_prefetch_path, arg_vdb_validate_path, arg_fastq_dump_path, arg_fastqc_path, arg_gzip_path, arg_cat_path, arg_bbduk_path, arg_megahit_path, arg_kraken_path, arg_bracken_path, arg_alpha_diversity_path, arg_metabinner_bin_path, arg_comebin_bin_path, arg_fraggenescanrs_path, arg_hmmscan_path, arg_hmmpress_path, arg_hmmfetch_path, arg_diamond_path, arg_cd_hit_path, arg_taxonkit_path, arg_phobius_path, arg_bowtie_build_path, arg_bowtie_path, arg_input_command)
    
