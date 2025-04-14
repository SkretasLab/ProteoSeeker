import re
import os
import shutil
# ProteoSeeker modules
import command_process
import supportive_functions


def diamond_first_round_fpd(diamond_env, output_path_hmmer, output_path_blastp, prs_with_enz_domains, blastp_results_no_doms_nr_file, blastp_info_no_doms_nr_file, blastp_no_doms_below_threshold, blastp_info_no_doms_below_threshold, cd_hit_results_path, file_prs_seq_no_enzs_name, fpr_db_fasta, fpr_db_name, e_value_nodom_thr, thread_num, input_log_file, output_log_file, diamond_path, blastp_nodoms_script, diamond_blastp_nr1_stdoe_path, conda_sh_path):
    print("\nRunning DIAMOND BLASTP for proteins without seek domains...")
    # This is the second time this folder is needed, so it has to be created only if it has not been already created,
    # but not overwrite the folder of this path, if it has already been created.
    if not os.path.exists(output_path_hmmer):
        os.mkdir(output_path_hmmer)
    # This is the first time this folder is needed, so it has to be created or overwrite this folder path.
    if os.path.exists(output_path_blastp):
        shutil.rmtree(output_path_blastp)
    os.mkdir(output_path_blastp)
    # Patterns
    pattern_blastp_info = re.compile(r'^(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)')
    # The list of proteins with the specified enzyme domains is either empty (if no analysi mode 1 occured) or it has
    # already been formed after analysis mode 1.
    # Write the proteins which did not have enzyme domains in a file.
    # If the analysis for pecified domains is not performed (or proteins with those domains have not been found) then
    # all of the proteins are written in a file and are further analyzed.
    prs_blast_thr = []
    new_file_prs_blast_thr = open(file_prs_seq_no_enzs_name, "w")
    with open(cd_hit_results_path) as cdhit_results_file:
        for line in cdhit_results_file:
            line = line.rstrip("\n")
            if line[0] == ">":
                enz_protein = False
                for i in prs_with_enz_domains:
                    if i in line:
                        enz_protein = True
            # The header of the protein sequence has already determined whether a enzyme domain has been found in the protein or not. Thus every line after the header will be written in the file
            # only if in the protein no enzyme domains had been found. Otherwise no line of the protein's sequence will be written in the file.
            if not enz_protein:
                if line[0] == ">":
                    prs_blast_thr.append(line)
                new_file_prs_blast_thr.write("{}\n".format(line))
    new_file_prs_blast_thr.close()
    # BLASTP for the proteins without enzyme domains.
    if os.path.exists(file_prs_seq_no_enzs_name):
        if os.path.exists(fpr_db_fasta):
            # Start time
            if diamond_path:
                phrase_1 = "\"{}\" blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(diamond_path, thread_num, fpr_db_name, blastp_results_no_doms_nr_file, file_prs_seq_no_enzs_name, diamond_blastp_nr1_stdoe_path)
                phrase_2 = ""
            else:
                phrase_1 = "diamond blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(thread_num, fpr_db_name, blastp_results_no_doms_nr_file, file_prs_seq_no_enzs_name, diamond_blastp_nr1_stdoe_path)
                phrase_2 = ""
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(blastp_nodoms_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if diamond_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(diamond_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            if phrase_2:
                new_file_bash.write("{}\n".format(phrase_2))
            if diamond_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(blastp_nodoms_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(blastp_nodoms_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        else:
            print("\nBLASTP of the proteins without any of the seek domains was not performed against the SFPD. The SFPD ({}) was not located.".format(fpr_db_fasta))
        # Collect results from BLASTP in the annotation folder.
        blastp_output_info = open(blastp_info_no_doms_nr_file, "w")
        blastp_output_info.write("Protein_name\tNR_ID\tPercentage_identity\tE_value\tBitscore\n")
        if os.path.exists(blastp_results_no_doms_nr_file):
            with open(blastp_results_no_doms_nr_file) as blastp_lines:
                for line in blastp_lines:
                    line = line.rstrip("\n")
                    result_blastp_info = pattern_blastp_info.search(line)
                    if result_blastp_info:
                        protein_name = result_blastp_info.group(1)
                        uniprot_ac = result_blastp_info.group(2)
                        # Insert check for name of entry with that Uniprot ID and possible GO annotations with it.
                        percentage_identity = result_blastp_info.group(3)
                        e_value = result_blastp_info.group(4)
                        bitscore = result_blastp_info.group(5)
                        blastp_output_info.write("{}\t{}\t{}\t{}\t{}\n".format(protein_name, uniprot_ac, percentage_identity, e_value, bitscore))
        blastp_output_info.close()
    else:
        print("\nBLASTP of the proteins without any of the seek domains was not performed against the SFPD. The file ({}) with the proteins without any of the seek domains, was not located.".format(file_prs_seq_no_enzs_name))
    # Information gathering, based on a given threshold, from the results of BLASTP of proteins without enzyme domains against the reductd nr.
    protein_ids_below_thr = []
    if os.path.exists(blastp_info_no_doms_nr_file):
        e_value_nodom_thr = float(e_value_nodom_thr)
        blastp_lines = supportive_functions.read_file(blastp_info_no_doms_nr_file)
        header = True
        header_write = True
        for line in blastp_lines:
            if header:
                header = False
                continue
            line_splited = line.split("\t")
            protein_acc = line_splited[0]
            e_value = float(line_splited[-2])
            # 0.00001 = 1e-5
            if e_value <= e_value_nodom_thr:
                if protein_acc not in protein_ids_below_thr:
                    protein_ids_below_thr.append(protein_acc)
                if header_write:
                    new_file = open(blastp_info_no_doms_below_threshold, "w")
                    new_file.write("Protein_name\tNR_ID\tPercentage_identity\tE_value\tBitscore\n")
                    header_write = False
                new_file.write("{}\n".format(line))
        # If the file opened then close it. This can be checked based on whether the header line for the new file was written in it or not.
        if not header_write:
            new_file.close()
        # Write the sequences of the proteins with no domains and e-value euqal to or smaller than the set threshold.
        if protein_ids_below_thr:
            new_file = open(blastp_no_doms_below_threshold, "w")
            with open(cd_hit_results_path) as lines_seqs:
                for line_pr_seq in lines_seqs:
                    line_pr_seq = line_pr_seq.rstrip("\n")
                    if line_pr_seq[0] == ">":
                        protein_boolean = False
                        protein_acc = line_pr_seq[1:]
                        # The protein_ids_below_thr list contain the IDs of all the proteins which were found not to have enzyme domains and whose best hit in NR had an e-value lower than the specified threshold.
                        if protein_acc in protein_ids_below_thr:
                            protein_boolean = True
                    if protein_boolean:
                        new_file.write("{}\n".format(line_pr_seq))
            new_file.close()
    return prs_blast_thr, protein_ids_below_thr


def diamond_second_round_swissprot_fpd(diamond_env, output_path_blastp, blastp_results_swissprot_file, blastp_doms_nr_file, blastp_info_swissprot_file, blastp_doms_info_nr_file, file_prs_seq_enz_domains_name, proteins_combined_file, swissprot_path, fpr_db_fasta, fpr_db_name, thread_num, input_log_file, output_log_file, diamond_path, blast_swissprot_bash_script, blast_nr_bash_script, phylo_analysis, conda_sh_path, diamond_blastp_swissprot_stdoe_path, diamond_blastp_nr2_stdoe_path):
    print("\nRunning DIAMOND for the combined dataset...")
    # Second time this folder is needed, so if it already exits, then it is not overwritten. It is creted only if it
    # has not already.
    if not os.path.exists(output_path_blastp):
        os.mkdir(output_path_blastp)
    # Initialiazing variables
    dict_sp = None
    # Pattern
    pattern_blastp_info = re.compile(r'^(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)')
    # If the file with the combined predictions for proteins exists.
    # Proteins from analysis mode 1 need to be blasted against the databases: SwissProt, nr
    # Proteins from analysis mode 2 need to be blasted against the databases: SwissProt
    # Therefore the combined proteines are blasted against SwissProt and only the proteins from analysis
    # mdoe 1 are blasted against nr and their results are added to the result of Blastp against nr of the proteins
    # from analysi mode 2.
    # ---------DIAMOND: combined proteins VS swissprot---------
    if os.path.exists(proteins_combined_file):
        # If the SwissProt database exists.
        if os.path.exists(swissprot_path):
            if diamond_path:
                phrase_1 = "\"{}\" blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(diamond_path, thread_num, swissprot_path, blastp_results_swissprot_file, proteins_combined_file, diamond_blastp_swissprot_stdoe_path)
                phrase_2 = ""
            else:
                phrase_1 = "diamond blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(thread_num, swissprot_path, blastp_results_swissprot_file, proteins_combined_file, diamond_blastp_swissprot_stdoe_path)
                phrase_2 = ""
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(blast_swissprot_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if diamond_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(diamond_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            if phrase_2:
                new_file_bash.write("{}\n".format(phrase_2))
            if diamond_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(blast_swissprot_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(blast_swissprot_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        else:
            print("\nBlastp was not performed against the Swiss-Prot database.")
        # Gathering information.
        # Collect results from DIAMOND in the annotation folder for the SwissProt database. Based on this database the protein codes are checked for the GO annotations.
        dict_sp = {}
        min_evalue = None
        new_file = open(blastp_info_swissprot_file, "w")
        new_file.write("Protein_name\tUniprot_AC\tPercentage_identity\tE_value\tBitscore\n")
        if os.path.exists(blastp_results_swissprot_file):
            with open(blastp_results_swissprot_file) as blastp_lines:
                for line in blastp_lines:
                    line = line.rstrip("\n")
                    result_blastp_info = pattern_blastp_info.search(line)
                    if result_blastp_info:
                        protein_name = result_blastp_info.group(1)
                        uniprot_ac = result_blastp_info.group(2)
                        # Insert check for name of entry with that Uniprot ID and possible GO annotations with it.
                        percentage_identity = result_blastp_info.group(3)
                        e_value = float(result_blastp_info.group(4))
                        bitscore = result_blastp_info.group(5)
                        new_line = "{}\t{}\t{}\t{}\t{}".format(protein_name, uniprot_ac, percentage_identity, e_value, bitscore)
                        new_file.write("{}\n".format(new_line))
                        # Check if the current e-value of the Uniprot AC for the protein name (if any) is smaller or higher. If higher continue. If smaller change the value of the
                        # dictionary for the protein_name with the current line of the Uniprot AC.
                        if protein_name not in dict_sp.keys():
                            dict_sp[protein_name] = new_line
                            min_evalue = e_value
                        else:
                            if e_value < min_evalue:
                                dict_sp[protein_name] = new_line
                                min_evalue = e_value
        new_file.close()
    # ---------DIAMOND: phylo domain proteins VS swissprot---------
    # DIAMOND of protein with phylo domains against the Swiss-Prot database.
    if phylo_analysis:
        # If the SwissProt database exists.
        if os.path.exists(swissprot_path) and os.path.exists(file_prs_seq_enz_domains_name):
            if diamond_path:
                phrase_1 = "\"{}\" blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(diamond_path, thread_num, swissprot_path, blastp_results_swissprot_file, file_prs_seq_enz_domains_name, diamond_blastp_swissprot_stdoe_path)
                phrase_2 = ""
            else:
                phrase_1 = "diamond blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(thread_num, swissprot_path, blastp_results_swissprot_file, file_prs_seq_enz_domains_name, diamond_blastp_swissprot_stdoe_path)
                phrase_2 = ""
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(blast_swissprot_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if diamond_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(diamond_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            if phrase_2:
                new_file_bash.write("{}\n".format(phrase_2))
            if diamond_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(blast_swissprot_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(blast_swissprot_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        else:
            print("\nBlastp was not performed against the Swiss-Prot database.")
        # Gathering information.
        # Collect results from DIAMOND in the annotation folder for the SwissProt database. Based on this database the protein codes are checked for the GO annotations.
        dict_sp = {}
        min_evalue = None
        new_file = open(blastp_info_swissprot_file, "w")
        new_file.write("Protein_name\tUniprot_AC\tPercentage_identity\tE_value\tBitscore\n")
        if os.path.exists(blastp_results_swissprot_file):
            with open(blastp_results_swissprot_file) as blastp_lines:
                for line in blastp_lines:
                    line = line.rstrip("\n")
                    result_blastp_info = pattern_blastp_info.search(line)
                    if result_blastp_info:
                        protein_name = result_blastp_info.group(1)
                        uniprot_ac = result_blastp_info.group(2)
                        # Insert check for name of entry with that Uniprot ID and possible GO annotations with it.
                        percentage_identity = result_blastp_info.group(3)
                        e_value = float(result_blastp_info.group(4))
                        bitscore = result_blastp_info.group(5)
                        new_line = "{}\t{}\t{}\t{}\t{}".format(protein_name, uniprot_ac, percentage_identity, e_value, bitscore)
                        new_file.write("{}\n".format(new_line))
                        # Check if the current e-value of the Uniprot AC for the protein name (if any) is smaller or higher. If higher continue. If smaller change the value of the
                        # dictionary for the protein_name with the current line of the Uniprot AC.
                        if protein_name not in dict_sp.keys():
                            dict_sp[protein_name] = new_line
                            min_evalue = e_value
                        else:
                            if e_value < min_evalue:
                                dict_sp[protein_name] = new_line
                                min_evalue = e_value
        new_file.close()
    # ---------DIAMOND: seek domain proteins VS swissprot---------
    # Specifically for the nr database, if the analysis mode if 2 or 3 only the proteins from mode analysis 1 are
    # used in BLAST.
    if os.path.exists(file_prs_seq_enz_domains_name):
        # If proteins from analysis mode 1 an the nr database exist, then run the proteins from analysis mode 1
        # against the nr database.
        if os.path.exists(fpr_db_fasta):
            if diamond_path:
                phrase_1 = "\"{}\" blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(diamond_path, thread_num, fpr_db_name, blastp_doms_nr_file, file_prs_seq_enz_domains_name, diamond_blastp_nr2_stdoe_path)
                phrase_2 = ""
            else:
                phrase_1 = "diamond blastp --threads {} --db \"{}\" --out \"{}\" --query \"{}\" --outfmt 6 &> \"{}\"".format(thread_num, fpr_db_name, blastp_doms_nr_file, file_prs_seq_enz_domains_name, diamond_blastp_nr2_stdoe_path)
                phrase_2 = ""
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(blast_nr_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if diamond_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(diamond_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            if phrase_2:
                new_file_bash.write("{}\n".format(phrase_2))
            if diamond_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(blast_nr_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(blast_nr_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Collect results from BLASTP in the annotation folder for the filtered nr database.
            new_file = open(blastp_doms_info_nr_file, "w")
            new_file.write("Protein_name\tNR_ID\tPercentage_identity\tE_value\tBitscore\n")
            if os.path.exists(blastp_doms_nr_file):
                with open(blastp_doms_nr_file) as blastp_lines:
                    for line in blastp_lines:
                        line = line.rstrip("\n")
                        result_blastp_info = pattern_blastp_info.search(line)
                        if result_blastp_info:
                            protein_name = result_blastp_info.group(1)
                            uniprot_ac = result_blastp_info.group(2)
                            # Insert check for name of entry with that Uniprot ID and possible GO annotations with it.
                            percentage_identity = result_blastp_info.group(3)
                            e_value = result_blastp_info.group(4)
                            bitscore = result_blastp_info.group(5)
                            new_file.write("{}\t{}\t{}\t{}\t{}\n".format(protein_name, uniprot_ac, percentage_identity, e_value, bitscore))
            new_file.close()
    return dict_sp


def combine_fpd_results(blastp_info_no_doms_below_threshold, blastp_doms_info_nr_file, blastp_info_comb_nr_file):
    print("\nCombining predictions from the screen against the seek filtered protein database (SFPD)...")
    if os.path.exists(blastp_info_no_doms_below_threshold) and os.path.exists(blastp_doms_info_nr_file):
        filenames = [blastp_info_no_doms_below_threshold, blastp_doms_info_nr_file]
    elif os.path.exists(blastp_info_no_doms_below_threshold):
        filenames = [blastp_info_no_doms_below_threshold]
    elif os.path.exists(blastp_doms_info_nr_file):
        filenames = [blastp_doms_info_nr_file]
    else:
        filenames = []
    # The first line (header) of the first file is the only one written in file from the first lines
    # of any of the files.
    dict_nr = {}
    min_evalue = None
    with open(blastp_info_comb_nr_file, 'w') as outfile:
        outfile.write("Protein_name\tFPD_ID\tPercentage_identity\tE_value\tBitscore\n")
        for i in range(0, len(filenames)):
            fname = filenames[i]
            header = True
            with open(fname) as infile:
                for line in infile:
                    if not header:
                        outfile.write(line)
                        line_splited = line.split("\t")
                        pr_name = line_splited[0]
                        e_value = float(line_splited[3])
                        if pr_name not in dict_nr.keys():
                            dict_nr[pr_name] = line
                            min_evalue = e_value
                        else:
                            if e_value < min_evalue:
                                dict_nr[pr_name] = line
                                min_evalue = e_value
                    header = False
    return dict_nr