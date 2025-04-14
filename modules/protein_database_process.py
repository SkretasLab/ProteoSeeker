import os
import protein_database_filter
# ProteoSeeker modules
import command_process


def process_pr_db(fpr_db_folder, diamond_db_bash_name, diamond_env, diamond_path, fpr_db_fasta, fpr_db_name, diamond_makedb_version_path, diamond_makedb_stdoe_path, thread_num, input_log_file, output_log_file, conda_sh_path):
    # Create the database with "makeblstdb", else with "diamond".
    # Create the Bash script.
    # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
    diamond_db_bash_script = "{}/{}".format(fpr_db_folder, diamond_db_bash_name)
    new_file_bash = open(diamond_db_bash_script, "w")
    phrase = "#!/bin/bash"
    new_file_bash.write("{}\n".format(phrase))
    phrase_1 = "cat \"{}/\"*.fasta > \"{}\"".format(fpr_db_folder, fpr_db_fasta)
    if diamond_env:
        phrase_s = "source \"{}\"".format(conda_sh_path)
        phrase_a = "conda activate \"{}\"".format(diamond_env)
        new_file_bash.write("{}\n".format(phrase_s))
        new_file_bash.write("{}\n".format(phrase_a))
    if diamond_path:
        phrase_2 = "\"{}\" makedb --threads {} --in \"{}\" --db \"{}\" &> \"{}\"".format(diamond_path, thread_num, fpr_db_fasta, fpr_db_name, diamond_makedb_stdoe_path)
        phrase_3 = 'echo "No version available for diamond makedb." > {}'.format(diamond_makedb_version_path)
    else:
        phrase_2 = "diamond makedb --threads {} --in \"{}\" --db \"{}\" &> \"{}\"".format(thread_num, fpr_db_fasta, fpr_db_name, diamond_makedb_stdoe_path)
        phrase_3 = 'echo "No version available for diamond makedb." > \"{}\"'.format(diamond_makedb_version_path)
    new_file_bash.write("{}\n".format(phrase_1))
    new_file_bash.write("{}\n".format(phrase_2))
    new_file_bash.write("{}\n".format(phrase_3))
    if diamond_env:
        phrase = "conda deactivate"
        new_file_bash.write("{}\n".format(phrase))
    new_file_bash.close()
    # Making the Bash script executable.
    # Sending command to run.
    phrase_b1 = "chmod +x {}".format(diamond_db_bash_script)
    phrase_b2 = "chmod --version"
    title_1 = "Making a Bash script executable:"
    title_2 = "Version of chmod:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Sending command to run.
    phrase_b1 = "bash {}".format(diamond_db_bash_script)
    phrase_b2 = "bash --version"
    title_1 = "Running bash:"
    title_2 = "Version of bash:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Deleting part fasta files.
    pd_folder_files = os.listdir(fpr_db_folder)
    for pd_file in pd_folder_files:
        pd_file_path = "{}/{}".format(fpr_db_folder, pd_file)
        if "_part_" in pd_file:
            os.remove(pd_file_path)
    

def create_pr_db(pr_names_dict, protein_db_path, diamond_db_bash_name, diamond_env, diamond_path, thread_num, pdf_threads, input_log_file, output_log_file, conda_sh_path):
    pr_names_analysis = pr_names_dict[0][0]
    pr_names_taxonomy = pr_names_dict[1][0]

    # Duplicate names are removed.
    if (not pr_names_analysis) and (not pr_names_taxonomy):
        print("\nNo protein names were found. The protein database will not be filtered.")
    else:
        fpr_db_fasta = pr_names_dict[0][1]
        fpr_db_fasta_taxonomy = pr_names_dict[1][1]
        fpr_db_folder = pr_names_dict[0][2]
        fpr_db_folder_taxonomy = pr_names_dict[1][2]
        fpr_db_name = pr_names_dict[0][3]
        fpr_db_name_taxonomy = pr_names_dict[1][3]
        fpr_db_fasta_prefix = pr_names_dict[0][4]
        fpr_db_fasta_prefix_taxonomy  = pr_names_dict[1][4]
        # Create the folder for the database.
        if pr_names_analysis:
            os.mkdir(fpr_db_folder)
        if pr_names_taxonomy:
            os.mkdir(fpr_db_folder_taxonomy)
        print("\nNames collected to create the seek pnr: {}".format(pr_names_analysis))
        print("Names collected to create the phylo pnr: {}".format(pr_names_taxonomy))
        if pr_names_analysis:
            fpr_seeknames_path = "{}/seek_names.txt".format(fpr_db_folder)
            fpr_seeknames_file = open(fpr_seeknames_path, "w")
            for sn_item in pr_names_analysis:
                fpr_seeknames_file.write("{}\n".format(sn_item))
            fpr_seeknames_file.close()
        if pr_names_taxonomy:
            fpr_taxonomynames_path = "{}/phylo_names.txt".format(fpr_db_folder_taxonomy)
            fpr_taxonomynames_file = open(fpr_taxonomynames_path, "w")
            for pn_item in pr_names_taxonomy:
                fpr_taxonomynames_file.write("{}\n".format(pn_item))
            fpr_taxonomynames_file.close()
        info_file_path=""
        # Lines in nr: 2.706.460.600
        # Proteins in nr: 439.976.609
        print("\nCreating the fpr databases...")
        db_info_dict = {}
        temp_0 = [fpr_db_fasta_prefix, pr_names_analysis, False, None]
        temp_1 = [fpr_db_fasta_prefix_taxonomy, pr_names_taxonomy, False, None]
        # The log file is created in the folder of the seek database if both seek and phylo databases are to be created. Otherwise, the log file
        # is created in the seek or the phylo folder depending on which of the two functionalities has been selected.
        if pr_names_analysis and pr_names_taxonomy:
            log_file_path = "{}/pd_filter_log.txt".format(fpr_db_folder)
        elif pr_names_analysis:
            log_file_path = "{}/pd_filter_log.txt".format(fpr_db_folder)
        elif pr_names_taxonomy:
            log_file_path = "{}/pd_filter_log.txt".format(fpr_db_folder_taxonomy)
        # Database information
        if pr_names_analysis:
            db_info_dict[0] = temp_0
        if pr_names_taxonomy:
            db_info_dict[1] = temp_1
        protein_database_filter.prfilter(protein_db_path, info_file_path, db_info_dict, log_file_path, pdf_threads)
        
    # Processing and creating the databases.
    if pr_names_analysis:
        diamond_makedb_version_path = "{}/diamond_makedb_version.txt".format(fpr_db_folder)
        diamond_makedb_stdoe_path = "{}/diamond_makedb_stdoe.txt".format(fpr_db_folder)
        process_pr_db(fpr_db_folder, diamond_db_bash_name, diamond_env, diamond_path, fpr_db_fasta, fpr_db_name, diamond_makedb_version_path, diamond_makedb_stdoe_path, thread_num, input_log_file, output_log_file, conda_sh_path)
    if pr_names_taxonomy:
        diamond_makedb_version_path = "{}/diamond_makedb_version.txt".format(fpr_db_folder_taxonomy)
        diamond_makedb_stdoe_path = "{}/diamond_makedb_stdoe.txt".format(fpr_db_folder_taxonomy)
        process_pr_db(fpr_db_folder_taxonomy, diamond_db_bash_name, diamond_env, diamond_path, fpr_db_fasta_taxonomy, fpr_db_name_taxonomy, diamond_makedb_version_path, diamond_makedb_stdoe_path, thread_num, input_log_file, output_log_file, conda_sh_path)
