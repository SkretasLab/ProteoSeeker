import os
import copy
import time
import shutil
import command_process
import supportive_functions


def bracken_filtering(kraken_threshold, krakentools_bash_script, alpha_diversity_path, alpha_diversity_version_path, alpha_diversity_stdoe_path, input_log_file, output_log_file, bracken_output_path, kraken_species_thr_path, bracken_filters_path):
    first_kt = True
    bracken_species_thr_dict = None
    bracken_filters_file = open(bracken_filters_path, "w")
    for kt in kraken_threshold:
        kraken_species_thr_dict_temp = {}
        if "." in kt:
            kt = float(kt)
        else:
            kt = int(kt)
        kraken_prop = None
        shannon_index = None
        bracken_filters_file.write("Filter: {}\n".format(kt))

        # Compute the threshold, if needed.
        if kt in [-1, -2]:
            # Create the Bash script.
            new_file_bash = open(krakentools_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if alpha_diversity_path:
                phrase_1 = "python \"{}\" -f \"{}\" -a Sh &> \"{}\"".format(alpha_diversity_path, bracken_output_path, alpha_diversity_stdoe_path)
                phrase_2 = "python \"{}\" -h > \"{}\"".format(alpha_diversity_path, alpha_diversity_version_path)
            else:
                phrase_1 = "alpha_diversity.py -f \"{}\" -a Sh &> \"{}\"".format(bracken_output_path, alpha_diversity_stdoe_path)
                phrase_2 = "alpha_diversity.py -h > \"{}\"".format(alpha_diversity_version_path)
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            new_file_bash.close()

            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(krakentools_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(krakentools_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

            # Get the Shannon Index.
            alpha_div_lines = supportive_functions.read_file(alpha_diversity_stdoe_path)
            shannon_line_label = "Shannon's diversity: "
            for line in alpha_div_lines:
                if shannon_line_label in line:
                    line_splited = line.split(shannon_line_label)
                    shannon_part = line_splited[1]
                    shannon_part = shannon_part.strip()
                    shannon_index = float(shannon_part)
                    shannon_index = round(shannon_index, 2)
            if kt == -1:
                if 0 <= shannon_index <= 2.5:
                    kraken_prop = 1
                elif 2.5 < shannon_index <= 4.5:
                    kraken_prop = 0.1
                elif 4.5 < shannon_index:
                    kraken_prop = 0
            elif kt == -2:
                if 0 <= shannon_index <= 2.5:
                    kraken_prop = 0.1
                elif 2.5 < shannon_index:
                    kraken_prop = 0
            if kraken_prop is not None:
                kraken_prop_float = float(kraken_prop)
            else:
                kraken_prop_float = "-"

            # Create a list with species above the threshold.
            if os.path.exists(bracken_output_path):
                with open(bracken_output_path) as info_lines:
                    for line in info_lines:
                        line = line.rstrip("\n")
                        line_splited = line.split("\t")
                        clade_type = line_splited[2]
                        if clade_type == "S":
                            taxonid = int(line_splited[1])
                            taxon_count = int(line_splited[5])
                            rel_abu = float(line_splited[6])
                            abu_perc = rel_abu * 100
                            if abu_perc >= kraken_prop:
                                kraken_species_thr_dict_temp[taxonid] = [taxon_count, abu_perc]
            bracken_filters_file.write("Shannon Index: {}\n".format(shannon_index))
            bracken_filters_file.write("Kraken percentage filter: {}%\n".format(kraken_prop_float))
        elif isinstance(kt, int) or isinstance(kt, float):
            if os.path.exists(bracken_output_path):
                with open(bracken_output_path) as info_lines:
                    for line in info_lines:
                        line = line.rstrip("\n")
                        line_splited = line.split("\t")
                        clade_type = line_splited[2]
                        if clade_type == "S":
                            taxonid = int(line_splited[1])
                            taxon_count = int(line_splited[5])
                            rel_abu = float(line_splited[6])
                            abu_perc = rel_abu * 100
                            if isinstance(kt, int):
                                if taxon_count >= kt:
                                    kraken_species_thr_dict_temp[taxonid] = [taxon_count, abu_perc]
                            elif isinstance(kt, float):
                                if abu_perc >= kt:
                                    kraken_species_thr_dict_temp[taxonid] = [taxon_count, abu_perc]
        else:
            print("\nAn error occured when filtering the species in the kraken report. Exiting.")
            exit()

        if kraken_prop is None:
            kraken_prop = "None"
            kraken_prop_temp = "None"
        else:
            kraken_prop_temp = round(kraken_prop, 5)
        abs_sum = 0
        perc_sum = 0
        kraken_species_thr_path_temp = "{}_{}_{}.tsv".format(kraken_species_thr_path, kt, kraken_prop_temp)
        new_file_species_thr = open(kraken_species_thr_path_temp, "w")
        for key in kraken_species_thr_dict_temp.keys():
            taxon_count = kraken_species_thr_dict_temp[key][0]
            abu_perc = kraken_species_thr_dict_temp[key][1]
            abs_sum += taxon_count
            perc_sum += abu_perc
            new_file_species_thr.write("{}\t{}\t{}\n".format(key, taxon_count, abu_perc))
        new_file_species_thr.write("Total\t{}\t{}\n".format(abs_sum, perc_sum))
        new_file_species_thr.close()

        # The kraken threshold used to create the bins based on the filtered species from kraken is the first threshold used in the kraken threshold list.
        if first_kt:
            bracken_species_thr_dict = copy.deepcopy(kraken_species_thr_dict_temp)
            first_kt = False
    bracken_filters_file.close()
    return bracken_species_thr_dict


def kraken(paired_end, tr_ex_file_paths_p, tr_ex_file_paths, output_path_kraken, kraken_db_path, kraken_results_path, kraken_report_path, kraken_threshold, kraken_species_path, kraken_species_thr_path, kraken_reads_path, kraken_taxname_path, conda_sh_path, kraken_env, kraken_path, kraken_bash_script, kraken_stde_path, kraken_version_path, kraken_memory_mapping, bracken_bash_script, bracken_path, bracken_env, bracken_output_path, bracken_report_path, bracken_length, bracken_level, bracken_threshold, bracken_stde_path, bracken_version_path, krakentools_bash_script, alpha_diversity_path, alpha_diversity_version_path, alpha_diversity_stdoe_path, thread_num, time_dict, bracken_filters_path, input_log_file, output_log_file, kraken_par, bracken_par):
    print("\nRunning Kraken2...")
    read_to_species_dict = {}
    taxid_to_species_dict = {}
    if kraken_memory_mapping:
        kraken_memory_mapping = "--memory-mapping "
    else:
        kraken_memory_mapping = ""
    start_time_kraken_spec = time.time()

    # Tool parameters. Add whitespaces to the user-defined parameters, if needed.
    if kraken_par is None:
        kraken_par = " --threads {} --use-names".format(thread_num)
    else:
        if kraken_par[0] != " ":
            kraken_par = " {}".format(kraken_par)

    if os.path.exists(output_path_kraken):
        shutil.rmtree(output_path_kraken)
    os.mkdir(output_path_kraken)
    # If the option "--output" is not provided to Kraken2, then the information that would be stored in the
    # file provided to that option will be printed in the standard output and should be caught by redirection
    # for stdout (">").
    # Create the Bash script.
    new_file_bash = open(kraken_bash_script, "w")
    phrase = "#!/bin/bash"
    new_file_bash.write("{}\n".format(phrase))
    if kraken_env:
        phrase_s = "source \"{}\"".format(conda_sh_path)
        phrase_a = "conda activate \"{}\"".format(kraken_env)
        new_file_bash.write("{}\n".format(phrase_s))
        new_file_bash.write("{}\n".format(phrase_a))
    if kraken_path:
        phrase_1 = "\"{}\" --db \"{}\" --output \"{}\" --report \"{}\" {}--paired{}".format(kraken_path, kraken_db_path, kraken_results_path, kraken_report_path, kraken_memory_mapping, kraken_par)
        phrase_2 = "\"{}\" --version > \"{}\"".format(kraken_path, kraken_version_path)
    else:
        phrase_1 = "kraken2 --db \"{}\" --output \"{}\" --report \"{}\" {}--paired{}".format(kraken_db_path, kraken_results_path, kraken_report_path, kraken_memory_mapping, kraken_par)
        phrase_2 = "kraken2 --version > \"{}\"".format(kraken_version_path)
    if paired_end:
        phrase_1 = "{} \"".format(phrase_1)
        tr_ex_file_paths_keys = list(tr_ex_file_paths_p.keys())
        for key_i in range(0, len(tr_ex_file_paths_keys)):
            key = tr_ex_file_paths_keys[key_i]
            if key_i + 1 == len(tr_ex_file_paths_keys):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths_p[key][0])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths_p[key][0])
        phrase_1 = "{} \"".format(phrase_1)
        for key_i in range(0, len(tr_ex_file_paths_p.keys())):
            key = tr_ex_file_paths_keys[key_i]
            if key_i + 1 == len(tr_ex_file_paths_keys):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths_p[key][1])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths_p[key][1])
    else:
        phrase_1 = "{} \"".format(phrase_1)
        for i in range(0, len(tr_ex_file_paths)):
            if i + 1 == len(tr_ex_file_paths):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths[i])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths[i])
    phrase_1 = "{} 2> \"{}\"".format(phrase_1, kraken_stde_path)
    new_file_bash.write("{}\n".format(phrase_1))
    new_file_bash.write("{}\n".format(phrase_2))
    if kraken_env:
        phrase = "conda deactivate"
        new_file_bash.write("{}\n".format(phrase))
    new_file_bash.close()

    # Making the Bash script executable.
    # Sending command to run.
    phrase_b1 = "chmod +x {}".format(kraken_bash_script)
    phrase_b2 = "chmod --version"
    title_1 = "Making a Bash script executable:"
    title_2 = "Version of chmod:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Sending command to run.
    phrase_b1 = "bash {}".format(kraken_bash_script)
    phrase_b2 = "bash --version"
    title_1 = "Running bash:"
    title_2 = "Version of bash:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    label = "Process of kraken2 process (specific):"
    elpased_time = supportive_functions.end_time_analysis(label, start_time_kraken_spec, output_log_file)
    time_dict["kraken_specific_time"] = elpased_time

    # Tool parameters. Add whitespaces to the user-defined parameters, if needed.
    if bracken_par is None:
        bracken_par = " -r \"{}\" -l \"{}\" -t \"{}\" ".format(bracken_length, bracken_level, bracken_threshold)
    else:
        if bracken_par[0] != " ":
            bracken_par = " {}".format(bracken_par)
        if bracken_par[-1] != " ":
            bracken_par = "{} ".format(bracken_par)

    # Run Bracken to get the abundancies of the species.
    # Create the Bash script.
    new_file_bash = open(bracken_bash_script, "w")
    phrase = "#!/bin/bash"
    new_file_bash.write("{}\n".format(phrase))
    if bracken_env:
        phrase_s = "source \"{}\"".format(conda_sh_path)
        phrase_a = "conda activate \"{}\"".format(bracken_env)
        new_file_bash.write("{}\n".format(phrase_s))
        new_file_bash.write("{}\n".format(phrase_a))
    if bracken_path:
        phrase_1 = "\"{}\" -d \"{}\" -i \"{}\" -o \"{}\" -w \"{}\"{}&> \"{}\"".format(bracken_path, kraken_db_path, kraken_report_path, bracken_output_path, bracken_report_path, bracken_par, bracken_stde_path)
        phrase_2 = "\"{}\" -v > \"{}\"".format(bracken_path, bracken_version_path)
    else:
        phrase_1 = "bracken -d \"{}\" -i \"{}\" -o \"{}\" -w \"{}\"{}&> \"{}\"".format(kraken_db_path, kraken_report_path, bracken_output_path, bracken_report_path, bracken_par, bracken_stde_path)
        phrase_2 = "bracken -v > \"{}\"".format(bracken_version_path)
    new_file_bash.write("{}\n".format(phrase_1))
    new_file_bash.write("{}\n".format(phrase_2))
    if bracken_env:
        phrase = "conda deactivate"
        new_file_bash.write("{}\n".format(phrase))
    new_file_bash.close()

    # Making the Bash script executable.
    # Sending command to run.
    phrase_b1 = "chmod +x {}".format(bracken_bash_script)
    phrase_b2 = "chmod --version"
    title_1 = "Making a Bash script executable:"
    title_2 = "Version of chmod:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Sending command to run.
    phrase_b1 = "bash {}".format(bracken_bash_script)
    phrase_b2 = "bash --version"
    title_1 = "Running bash:"
    title_2 = "Version of bash:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    # Store the Kraken2 species and their relative abundances.
    kraken_species_dict = {}
    if os.path.exists(bracken_output_path):
        with open(bracken_output_path) as report_lines:
            for line in report_lines:
                line = line.rstrip("\n")
                line_splited = line.split("\t")
                clade_type = line_splited[2]
                if clade_type == "S":
                    taxonid = int(line_splited[1])
                    taxon_count = int(line_splited[5])
                    rel_abu = float(line_splited[6])
                    kraken_species_dict[taxonid] = [taxon_count, rel_abu]

    # Kraken filtering thresholds are applied.
    bracken_species_thr_dict = bracken_filtering(kraken_threshold, krakentools_bash_script, alpha_diversity_path, alpha_diversity_version_path, alpha_diversity_stdoe_path, input_log_file, output_log_file, bracken_output_path, kraken_species_thr_path, bracken_filters_path)
    # Check for the results
    if os.path.exists(kraken_results_path):
        with open(kraken_results_path) as results_lines:
            for line in results_lines:
                line = line.rstrip("\n")
                line_splited = line.split("\t")
                result_type = line_splited[0]
                if result_type == "C":
                    read_id = line_splited[1]
                    taxonid_phrase = line_splited[2]
                    if "taxid " in taxonid_phrase:
                        taxonid_phrase_splited = taxonid_phrase.split("taxid ")
                        taxonid_phrase_part = taxonid_phrase_splited[1]
                        if taxonid_phrase_part[-1] == ")":
                            taxonid_phrase_part = taxonid_phrase_part[:-1]
                            taxonid = int(taxonid_phrase_part)
                            # The read_to_species_dict is based on the dictionary (bracken_species_thr_dict) returned by the
                            # filtering which is the one generated by the first filter applied to the species of the kraken report.
                            # The binning based on the taxonomy assignment of kraken is based on the read_to_species_dict dictionary.
                            # Thus the binning process is based on the information for the associations of reads and species which species
                            # were classified and were above the set threshold of the first filtering threshold used by ProteoSeeker.
                            # Each read has a list as a value which contains the taxid associated with the read and whether or not the
                            # taxid passed the filtering threshold (1) or not (0).
                            read_to_species_dict[read_id] = [taxonid, 0]
                            if taxonid in bracken_species_thr_dict.keys():
                                read_to_species_dict[read_id][1] = 1
                            # The taxid_to_species_dict is not depended on the bracken_species_thr_dict dictionary.
                            if taxonid not in taxid_to_species_dict.keys():
                                taxname = taxonid_phrase_splited[0]
                                if taxname[-1] == "(":
                                    taxname = taxname[:-1]
                                taxname = taxname.strip()
                                taxid_to_species_dict[taxonid] = taxname

    new_file_species = open(kraken_species_path, "w")
    for key in kraken_species_dict.keys():
        taxon_count = kraken_species_dict[key][0]
        abu_perc = kraken_species_dict[key][1]
        new_file_species.write("{}\t{}\t{}\n".format(key, taxon_count, abu_perc))
    new_file_species.close()

    new_file_reads = open(kraken_reads_path, "w")
    for key in read_to_species_dict.keys():
        taxonid = read_to_species_dict[key][0]
        thr_state = read_to_species_dict[key][1]
        new_file_reads.write("{}\t{}\t{}\n".format(key, taxonid, thr_state))
    new_file_reads.close()

    new_file_taxname = open(kraken_taxname_path, "w")
    for key_tax in taxid_to_species_dict.keys():
        tax_name = taxid_to_species_dict[key_tax]
        new_file_taxname.write("{}\t{}\n".format(key_tax, tax_name))
    new_file_taxname.close()
    return kraken_species_dict, bracken_species_thr_dict, read_to_species_dict, taxid_to_species_dict, time_dict