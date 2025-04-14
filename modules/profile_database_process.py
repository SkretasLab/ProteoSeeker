import os
# ProteoSeeker modules
import command_process
import supportive_functions


def craete_phmm_db(family_code, family_group_name, enzyme_domains_folder, fam_nums_file_name, fam_pfam_file_name, pfam_domains_names, hmmer_env, hmmfetch_path, hmmpress_path, profiles_broad_path, input_log_file, output_log_file, conda_sh_path):
    # Initialization of variables.
    d_unique = None
    profiles_name = None
    analysis_fam_names = []
    family_to_profile_dict = {}
    # Folder and file paths. Some are initialized only if at least one protein family was selected.
    phmm_database_folder_name = "{}/profile_dbs/phmm_db_{}".format(enzyme_domains_folder, family_group_name)
    fam_prof_info_path = "{}/family_profile_info.tsv".format(phmm_database_folder_name)
    if family_group_name is not None:
        new_file_dom_codes_name = "{}/profile_codes.txt".format(phmm_database_folder_name)
        new_file_dom_names_name = "{}/profile_names.txt".format(phmm_database_folder_name)
        profiles_name = "{}/profiles".format(phmm_database_folder_name)
        hmmer_db_bash_script = "{}/hmmer_db.sh".format(phmm_database_folder_name)
        hmmfetch_stdoe_path = "{}/hmmfetch_stdoe.txt".format(phmm_database_folder_name)
        hmmfetch_version_path = "{}/hmmfetch_version.txt".format(phmm_database_folder_name)
        hmmpress_stdoe_path = "{}/hmmpress_stdoe.txt".format(phmm_database_folder_name)
        hmmpress_version_path = "{}/hmmpress_version.txt".format(phmm_database_folder_name)
    # For each protein family selected...
    if family_code is not None:
        d_unique = []
        # Create the correspondances between the families and their numbering.
        fam_nums_lines = supportive_functions.read_file(fam_nums_file_name)
        fam_num_dict = {}
        for line in fam_nums_lines:
            line_splited = line.split("\t")
            fam_num = line_splited[0]
            fam_name = line_splited[1]
            fam_num_dict[fam_num] = fam_name
        for family_code_cur in family_code:
            # Find the name of the selected protein family.
            # This must be changed to a dictionary and then search the protein length against the all mean legnth of each selected protein family and give out the family with the closest length.
            sel_fam = fam_num_dict[family_code_cur]
            print("Selected faimly: {}".format(sel_fam))
            analysis_fam_names.append(sel_fam)
            if sel_fam not in family_to_profile_dict.keys():
                family_to_profile_dict[sel_fam] = {}
            # The domain sets for each protein family are already those of the highest requency for the protein family.
            # When selecting multiple protein families, all protein sets from each protein family are collected, because the point is to
            # search for the domains corresponding to each of the selected protein families as to later on determine proteins that may belong
            # to any of the selected protein families.
            # Collect the domain and their counts, corresponding to the name of the selected protein family.
            # Collect each set of domains related to the protein family with their counts.
            # Collect the unique domains related to the protein family.
            pfam_lines = supportive_functions.read_file(fam_pfam_file_name)
            for line in pfam_lines:
                line_splited = line.split("\t")
                line_fam = line_splited[0]
                if line_fam == sel_fam:
                    domain_sets = line_splited[1]
                    # Remove the first and last parenthesis of the outer list.
                    domain_sets = domain_sets[1:-1]
                    # Remove the first parenthesis of the first list and the last parenthesis of the last list.
                    # These lists might be the same.
                    domain_sets = domain_sets[1:-1]
                    # Now split the domains sets.
                    # If the domains sets are more than one then split them.
                    if "], [" in domain_sets:
                        # 'PF00541_1;PF00608_2'], ['PF00541_1'], ['PF00608_7'
                        domain_sets_splited = domain_sets.split("], [")
                        for ds in domain_sets_splited:
                            # Remove the single quotes.
                            ds = ds[1:-1]
                            # If the set of domains contains multiple domains split them.
                            if ";" in ds:
                                # PF00541_1;PF00608_2
                                ds_splited = ds.split(";")
                                # Extract the domain codes and their counters.
                                # ["PF00541_1", "PF00608_2"]
                                for item in ds_splited:
                                    item_splited = item.split("_")
                                    dname = item_splited[0]
                                    if dname not in d_unique:
                                        d_unique.append(dname)
                                    family_to_profile_dict[sel_fam][dname] = 0
                            else:
                                # PF00541_1
                                ds_splited = ds.split("_")
                                dname = ds_splited[0]
                                if dname not in d_unique:
                                    d_unique.append(dname)       
                                family_to_profile_dict[sel_fam][dname] = 0                 
                    else:
                        # 'PF01757_1;PF00541_1'
                        # Remove the single quotes.
                        domain_sets = domain_sets[1:-1]
                        # PF01757_1;PF00541_1
                        if ";" in domain_sets:
                            ds_splited = domain_sets.split(";")
                            # ["PF01757_1", "PF00541_1"]
                            for item in ds_splited:
                                item_splited = item.split("_")
                                dname = item_splited[0]
                                if dname not in d_unique:
                                    d_unique.append(dname)
                                family_to_profile_dict[sel_fam][dname] = 0
                        else:
                            # PF00541_1
                            ds_splited = domain_sets.split("_")
                            dname = ds_splited[0]
                            if dname not in d_unique:
                                d_unique.append(dname)
                            family_to_profile_dict[sel_fam][dname] = 0
        print("\nSelected domain codes: {}".format(d_unique))
        # If the path to the folder of the pHMM library does not exist, then create the folder and then the library.
        if not os.path.exists(phmm_database_folder_name):
            os.mkdir(phmm_database_folder_name)
            # Creating a file to store the families and their profiles.
            fam_prof_info_file = open(fam_prof_info_path, "w")
            for key_fam in family_to_profile_dict.keys():
                for key_prof in family_to_profile_dict[key_fam].keys():
                    fam_prof_info_file.write("{}\t{}\n".format(key_fam, key_prof))
            fam_prof_info_file.close()
            # Create dictionary of Pfam codes to names.
            # The Pfam code is correponded to the latest Pfam version of the Pfam code.
            pfam_names_doms_dict = {}
            pfam_doms_names_lines = supportive_functions.read_file(pfam_domains_names)
            for line in pfam_doms_names_lines:
                splited_line = line.split("\t")
                dom_code = splited_line[0]
                dom_name = splited_line[1]
                if "." in dom_code:
                    dom_code_splited = dom_code.split(".")
                    dom_code = dom_code_splited[0]
                pfam_names_doms_dict[dom_code] = dom_name
            # Find the current Pfam codes and create a file with these codes.
            new_file_dom_codes = open(new_file_dom_codes_name, "w")
            new_file_dom_names = open(new_file_dom_names_name, "w")
            domains_selected = []
            for dom_code in d_unique:
                if dom_code in pfam_names_doms_dict.keys():
                    dom_name = pfam_names_doms_dict[dom_code]
                    domains_selected.append(dom_name)
                    new_file_dom_codes.write("{}\n".format(dom_code))
                    new_file_dom_names.write("{}\n".format(dom_name))
            new_file_dom_codes.close()
            new_file_dom_names.close()
            # Createa a pHMM file with the corresponding profile names.
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(hmmer_db_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if hmmer_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(hmmer_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            if hmmfetch_path:
                phrase_1 = "\"{}\" -o \"{}\" -f \"{}\" \"{}\" &> \"{}\"".format(hmmfetch_path, profiles_name, profiles_broad_path, new_file_dom_names_name, hmmfetch_stdoe_path)
                phrase_2 = "\"{}\" \"{}\" &> \"{}\"".format(hmmpress_path, profiles_name, hmmpress_stdoe_path)
                phrase_3 = "\"{}\" -h > \"{}\"".format(hmmfetch_path, hmmfetch_version_path)
                phrase_4 = "\"{}\" -h > \"{}\"".format(hmmpress_path, hmmpress_version_path)
            else:
                phrase_1 = "hmmfetch -o \"{}\" -f \"{}\" \"{}\" &> \"{}\"".format(profiles_name, profiles_broad_path, new_file_dom_names_name, hmmfetch_stdoe_path)
                phrase_2 = "hmmpress \"{}\" &> \"{}\"".format(profiles_name, hmmpress_stdoe_path)
                phrase_3 = "hmmfetch -h > \"{}\"".format(hmmfetch_version_path)
                phrase_4 = "hmmpress -h > \"{}\"".format(hmmpress_version_path)
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            new_file_bash.write("{}\n".format(phrase_3))
            new_file_bash.write("{}\n".format(phrase_4))
            if hmmer_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(hmmer_db_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(hmmer_db_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        else:
            if os.path.exists(fam_prof_info_path):
                fam_prof_info_lines = supportive_functions.read_file(fam_prof_info_path)
                for line in fam_prof_info_lines:
                    line_splited = line.split("\t")
                    fam_temp = line_splited[0]
                    prof_temp = line_splited[1]
                    if fam_temp not in family_to_profile_dict.keys():
                        family_to_profile_dict[fam_temp] = {}
                    family_to_profile_dict[fam_temp][prof_temp] = 0
            print("\nThe pHMM database for this protein family already exists. Skipping creation of the pHMM database.")
    else:
        if os.path.exists(fam_prof_info_path):
            fam_prof_info_lines = supportive_functions.read_file(fam_prof_info_path)
            for line in fam_prof_info_lines:
                line_splited = line.split("\t")
                fam_temp = line_splited[0]
                prof_temp = line_splited[1]
                if fam_temp not in family_to_profile_dict.keys():
                    family_to_profile_dict[fam_temp] = {}
                family_to_profile_dict[fam_temp][prof_temp] = 0
    return profiles_name, analysis_fam_names, family_to_profile_dict