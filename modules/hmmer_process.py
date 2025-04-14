import re
import os
import shutil
# ProteoSeeker modules
import command_process


def hmmer_spec(hmmer_env, output_path_hmmer, hmmer_dmbl_results, hmmer_simple_results, hmmer_enz_domains_all_proteins, cd_hit_results_path, hmmscan_path, profiles_path, val_type, file_prs_seq_enz_domains_name, thread_num, input_log_file, output_log_file, hmmer_spec_bash_script, phylo_analysis, hmmscan_spec_version_path, hmmscan_spec_stdoe_path, conda_sh_path):
    print("\nRunning HMMER...")
    if not phylo_analysis:
        if os.path.exists(output_path_hmmer):
            shutil.rmtree(output_path_hmmer)
        os.mkdir(output_path_hmmer)
    else:
        if not os.path.exists(output_path_hmmer):
            os.mkdir(output_path_hmmer)
    # If the proteins used for the analysis, after CD-HIT are not found by FragGeneScanRs then it will be needed to add a step where the headers of the proteins will not have empty characters because HMMER only reports the headers
    # of the protein to which it finds domains based only on their part of the name until it meets the first empty character. FragGeneScanRs produces names without empty characters so there is no problem currently. If another gene prediction program is used there might be.
    if os.path.exists(cd_hit_results_path) and os.path.exists(profiles_path):
        if hmmscan_path:
            phrase_1 = "\"{}\" -o \"{}\" --domtblout \"{}\" --notextw {}--cpu {} \"{}\" \"{}\" &> \"{}\"".format(hmmscan_path, hmmer_simple_results, hmmer_dmbl_results, val_type, thread_num, profiles_path, cd_hit_results_path, hmmscan_spec_stdoe_path)
            phrase_2 = "\"{}\" -h > \"{}\"".format(hmmscan_path, hmmscan_spec_version_path)
        else:
            phrase_1 = "hmmscan -o \"{}\" --domtblout \"{}\" --notextw {}--cpu {} \"{}\" \"{}\" &> \"{}\"".format(hmmer_simple_results, hmmer_dmbl_results, val_type, thread_num, profiles_path, cd_hit_results_path, hmmscan_spec_stdoe_path)
            phrase_2 = "hmmscan -h > \"{}\"".format(hmmscan_spec_version_path)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(hmmer_spec_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if hmmer_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(hmmer_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        if hmmer_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(hmmer_spec_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(hmmer_spec_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Processing the output from hmmscan.
    pattern_domains = re.compile(r'^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(.*)')
    # It is mandatory to write all the domains found for one protein in contineous lines, that is the reason for using a dictionary.
    domains_dict_spec = {}
    prs_with_enz_domains = []
    if os.path.exists(hmmer_dmbl_results):
        with open(hmmer_dmbl_results) as db_hmmer_file:
            for line in db_hmmer_file:
                line = line.rstrip("\n")
                if line[0] != "#":
                    result_domains = pattern_domains.search(line)
                    if result_domains:
                        pfam_name_short = result_domains.group(1)
                        pfam_accession = result_domains.group(2)
                        pfam_name_long = result_domains.group(12)
                        protein_name = result_domains.group(3)
                        e_value = result_domains.group(4)
                        score = result_domains.group(5)
                        c_e_value = result_domains.group(6)
                        i_e_value = result_domains.group(7)
                        hmm_from = result_domains.group(8)
                        hmm_to = result_domains.group(9)
                        env_from = result_domains.group(10)
                        env_to = result_domains.group(11)
                        line_w = ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(protein_name, e_value, score, pfam_name_short, pfam_accession, pfam_name_long, c_e_value, i_e_value, hmm_from, hmm_to, env_from, env_to))
                        if protein_name not in domains_dict_spec.keys():
                            domains_dict_spec[protein_name] = [line_w]
                        else:
                            domains_dict_spec[protein_name].append(line_w)
                        if protein_name not in prs_with_enz_domains:
                            prs_with_enz_domains.append(protein_name)
    # The information are written in a file.
    if domains_dict_spec.keys():
        new_file_hmmer_output = open(hmmer_enz_domains_all_proteins, "w")
        new_file_hmmer_output.write("Protein_name\tE_value\tScore\tPfam_domain_short_name\tPfam_domain_accession_name\tPfam_domain_long_name\tc-E-value\ti-E-value\thmm_from\thmm_to\tenv_from\tenv_to\n")
        for key_pr in domains_dict_spec.keys():
            for domains_line in domains_dict_spec[key_pr]:
                new_file_hmmer_output.write("{}\n".format(domains_line))
        new_file_hmmer_output.close()
    else:
        print("\nNo proteins found with enzyme domains.\n")
    # Write the proteins found with the specified enzyme domains in a file.
    if prs_with_enz_domains:
        new_file = open(file_prs_seq_enz_domains_name, "w")
        with open(cd_hit_results_path) as cdhit_results_file:
            for line in cdhit_results_file:
                line = line.rstrip("\n")
                if line[0] == ">":
                    # The status showing a protein with any of the specified domains for a new protein sequence is
                    # initialized to False.
                    enz_pr_status = False
                    # If any header/name of the proteins found to have such domains is present in the current header,
                    # then this sequences is a protein with at least one of the specified domains.
                    for prd in prs_with_enz_domains:
                        if prd in line:
                            enz_pr_status = True
                if enz_pr_status:
                    new_file.write("{}\n".format(line))
        new_file.close()
    return prs_with_enz_domains



def hmmer_broad(hmmer_env, prs_with_enz_domains, hmmer_dmbl_results, hmmer_simple_results, hmmer_enz_domains_all_proteins, proteins_combined_file, comb_dmbl_results, comb_simple_results, comb_all_domains_proteins, profiles_broad_path, hmmscan_path, val_type, second_dom_search, thread_num, input_log_file, output_log_file, hmmer_broad_bash_script, hmmscan_broad_version_path, hmmscan_broad_stdoe_path, conda_sh_path):
    print("\nRunning HMMER against Pfam...")
    # Example:
    #                                                                                    --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
    # target name        accession   tlen query name                   accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description
    # (Pro_CA)               (PF00484.21)   156 (contig270_439198_it3_meta_344131_344436_-) -            101   (3.5e-11)   (32.7)   0.1   1   1   4.7e-12   3.7e-11   32.7   0.1    38    85    12    55     1    86 0.74 (Carbonic anhydrase)
    pattern_domains = re.compile(r'^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(.*)')
    # Rerun HMMER against all HMMs only for the proteins identified with any of the known enzyme domains.
    # File with protein sequences which have enzyme domains.
    # The file with all the representative sequences from CD-HIT is parsed and every header of the file is checked
    # on whether it contains any of the headers of the proteins found with enzyme domains. If such a protein is found
    # then it is removed from the list with the headers of proteins found to have enzyme domains as not to be checked
    # again in any of the lines of the (CD-HIT) file.
    dict_hm = {}
    if os.path.exists(proteins_combined_file):
        # Run against all HMMs of the second database (if any).
        # In case someone nees to use all the domains available from Pfam at the first search then there is no need to make
        # the second search. Therefore, in that case either the second database is left empty or the option second_dom_search
        # is to false.
        if (os.path.exists(profiles_broad_path)) and (second_dom_search is not False):
            if hmmscan_path:
                phrase_1 = "\"{}\" -o \"{}\" --domtblout \"{}\" --notextw {}--cpu {} \"{}\" \"{}\" &> \"{}\"".format(hmmscan_path, comb_simple_results, comb_dmbl_results, val_type, thread_num, profiles_broad_path, proteins_combined_file, hmmscan_broad_stdoe_path)
                phrase_2 = "\"{}\" -h > \"{}\"".format(hmmscan_path, hmmscan_broad_version_path)
            else:
                phrase_1 = "hmmscan -o \"{}\" --domtblout \"{}\" --notextw {}--cpu {} \"{}\" \"{}\" &> \"{}\"".format(comb_simple_results, comb_dmbl_results, val_type, thread_num, profiles_broad_path, proteins_combined_file, hmmscan_broad_stdoe_path)
                phrase_2 = "hmmscan -h > \"{}\"".format(hmmscan_broad_version_path)
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(hmmer_broad_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if hmmer_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(hmmer_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            if hmmer_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(hmmer_broad_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(hmmer_broad_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Recollect information
            # The information to collect from the hits of all domains in the proteins with enzyme domain are set here.
            # It is mandatory to write all the domains found for one protein in contineous lines, that is the reason for using a dictionary.
            if os.path.exists(comb_dmbl_results):
                with open(comb_dmbl_results) as db_hmmer_file:
                    for line in db_hmmer_file:
                        line = line.rstrip("\n")
                        if line[0] != "#":
                            result_domains = pattern_domains.search(line)
                            if result_domains:
                                pfam_name_short = result_domains.group(1)
                                pfam_accession = result_domains.group(2)
                                pfam_name_long = result_domains.group(12)
                                protein_name = result_domains.group(3)
                                e_value = result_domains.group(4)
                                score = result_domains.group(5)
                                c_e_value = result_domains.group(6)
                                i_e_value = result_domains.group(7)
                                hmm_from = result_domains.group(8)
                                hmm_to = result_domains.group(9)
                                env_from = result_domains.group(10)
                                env_to = result_domains.group(11)
                                line_w = ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(protein_name, e_value, score, pfam_name_short, pfam_accession, pfam_name_long, c_e_value, i_e_value, hmm_from, hmm_to, env_from, env_to))
                                if protein_name not in dict_hm.keys():
                                    dict_hm[protein_name] = [line_w]
                                else:
                                    dict_hm[protein_name].append(line_w)
            # The information are written in a file.
            if dict_hm.keys():
                new_file = open(comb_all_domains_proteins, "w")
                new_file.write("Protein_name\tE_value\tScore\tPfam_domain_short_name\tPfam_domain_accession_name\tPfam_domain_long_name\tc-E-value\ti-E-value\thmm_from\thmm_to\tenv_from\tenv_to\n")
                for key_pr in dict_hm.keys():
                    for domains_line in dict_hm[key_pr]:
                        new_file.write("{}\n".format(domains_line))
                new_file.close()
        else:
            # If no such database exists, then just copy the results of the search based on the firt profile database.
            if os.path.exists(hmmer_dmbl_results):
                shutil.copyfile(hmmer_dmbl_results, comb_dmbl_results)
            if os.path.exists(hmmer_simple_results):
                shutil.copyfile(hmmer_simple_results, comb_simple_results)
            if os.path.exists(hmmer_enz_domains_all_proteins):
                shutil.copyfile(hmmer_enz_domains_all_proteins, comb_all_domains_proteins)
    return prs_with_enz_domains, dict_hm