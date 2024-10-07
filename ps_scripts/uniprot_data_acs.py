import re
import csv
import sys
import copy
from statistics import mean, median


def help_par(text, step):
    words = text.split(" ")
    lines = [words[0]]
    for word in words[1:]:
        if len(lines[-1]) + len(word) < step:
            lines[-1] += (" " + word)
        else:
            lines.append(word)
    return lines


def help_message():
    print("uniprot_data_acs.py Version 1.0.0")
    print()
    print("Usage:")
    print("python uniprot_data_acs.py -i <swiss_prot_file>")
    print()
    print("This tool can be used to process information from the Swiss-Prot database downloaded as a file from UniprotKB.")
    print("It outputs a series of files each containing different types of information based on the Swiss-Prot file.")
    print("The output files are created in the same directory where this Python script is run.")
    print()
    print("Option description:")
    print("1. Parameter type")
    print("2. Req: Required, Opt: Optional")
    print("3. Default value (shown if not empty or none)")
    print("4. Description")
    print()
    print("Options:")
    print("---------Input and output options---------")
    help_notes_dict = {
        "-i/--input": "Str -Req: uniprot_sprot.dat- The path to the Swiss-Prot file.",
        "-c/--combination": "Int -Opt: 1- The combination type for the protein families identified from the Swiss-Prot file."
    }
    split_step = 60
    for key_hn in help_notes_dict.keys():
        description = help_notes_dict[key_hn]
        des_length = len(description)
        if split_step >= des_length:
            print('{:2} {:30} {}\n'.format("", key_hn, description))
        else:
            pieces = help_par(description, split_step)
            for pmi in range(0, len(pieces)):
                piece_mod = pieces[pmi]
                if pmi == 0:
                    print('{:2} {:30} {}'.format("", key_hn, piece_mod))
                elif pmi + 1 == len(pieces):
                    print('{:33} {}\n'.format("", piece_mod))
                else:
                    print('{:33} {}'.format("", piece_mod))


def read_file(file_path):
    file_handle = open(file_path, "r")
    pre_lines = file_handle.readlines()
    lines = []
    for i in pre_lines:
        i = i.rstrip("\n")
        lines.append(i)
    file_handle.close()
    return lines


def unidat(uniprot_falt_file="uniprot_sprot.dat", family_combination_type=1):
    pattern_id = re.compile(r'^ID\s+(\S+)')
    pattern_length = re.compile(r'^ID\s+\S+\s+\S+\s+(\S+)')
    pattern_ac = re.compile(r'^AC\s+(.*)')
    pattern_sim = re.compile(r'^CC\s+-!-\s+SIMILARITY:(.*)')
    pattern_sim_cont = re.compile(r'^CC\s+(.*)')
    pattern_interpro = re.compile(r'^DR\s+InterPro;\s+(\S+)\s+(\S+)')
    pattern_pfam = re.compile(r'^DR\s+Pfam;\s+(\S+)\s+(\S+)\s+(\S+)')
    # The "RecName" is followed always by "Full".
    pattern_de_recname_full = re.compile(r'^DE\s+RecName:\s+Full=([^{]*)')
    pattern_de_altname_full = re.compile(r'^DE\s+AltName:\s+Full=([^{]*)')
    pattern_de_empty_ec = re.compile(r'^DE\s+EC=([^{]*)')
    pattern_de_empty_short = re.compile(r'^DE\s+Short=([^{]*)')

    protein_information = {}
    sim_found = False

    pr_counter = 0

    with open(uniprot_falt_file) as unf:
        for line in unf:
            # Skip all the other lines not starting from the desirable cahracters to save time.
            if line[0:2] in ["AC", "ID", "DE", "CC", "DR"]:
                line = line.rstrip("\n")
                # Find the ID of the protein. Each protein has a unique ID.
                if line[0:2] == "ID":
                    result_id = pattern_id.search(line)
                    if result_id:
                        protein_id = result_id.group(1)
                        if protein_id not in protein_information.keys():
                            # A counter for the different proteins of the file. It is based on their unique protein IDs.
                            pr_counter += 1
                            if pr_counter % 10000 == 0:
                                print("Parsed {} IDs out of 569793.".format(pr_counter))
                            # The boolean that restarts the search for similarity lines.
                            sim_found = False
                            first_ac = True
                            protein_first_ac = None
                            protein_information[protein_id] = {
                                "id": protein_id,
                                "acs": [],
                                "names": []
                            }
                        else:
                            print("Duplicate IDs found for different proteins. Exiting...")
                            exit()
                    # Determine the length of the protein.
                    result_length = pattern_length.search(line)
                    if result_length:
                        protein_length = int(result_length.group(1))
                        protein_information[protein_id]["length"] = protein_length

                if line[0:2] == "AC":
                    # Find the latest AC of the protein. One protein may have had several different ACs.
                    result_ac = pattern_ac.search(line)
                    if result_ac:
                        protein_ac = result_ac.group(1)
                        # Split based on ";".
                        protein_ac_splited = protein_ac.split(";")
                        # Remove the last empty element, if any.
                        if protein_ac_splited[-1] == "":
                            protein_ac_splited = protein_ac_splited[:-1]
                        # Remove spaces from the start and end of the each string.
                        # If the first line of ACs has already been found then search for the entry with the key of the first AC of the protein, otherwie the protein ID is used.
                        for item in protein_ac_splited:
                            item = item.lstrip()
                            item = item.rstrip()
                            if first_ac:
                                protein_information[protein_id]["acs"].append(item)
                            else:
                                protein_information[protein_first_ac]["acs"].append(item)
                    # When the first line with AC(s) for a protein is found. Change the key of an entry (ID) with the first AC found for the protein.
                    if first_ac:
                        protein_first_ac = protein_information[protein_id]["acs"][0]
                        protein_information[protein_first_ac] = copy.deepcopy(protein_information[protein_id])
                        del protein_information[protein_id]
                        first_ac = False

                if line[0:2] == "DE":
                    # Find the names of the protein.
                    # Recname Full
                    result_de_recname_full = pattern_de_recname_full.search(line)
                    if result_de_recname_full:
                        recname_full = result_de_recname_full.group(1)
                        # Removing whitespaces before and after the matched phrase.
                        recname_full = recname_full.strip()
                        # Remove the semicolon after the phrase, if any.
                        if recname_full[-1] == ";":
                            recname_full = recname_full[:-1]
                        if recname_full not in protein_information[protein_first_ac]["names"]:
                            protein_information[protein_first_ac]["names"].append(recname_full)
                    # Altname Full
                    result_de_altname_full = pattern_de_altname_full.search(line)
                    if result_de_altname_full:
                        altname_full = result_de_altname_full.group(1)
                        # Removing whitespaces before and after the matched phrase.
                        altname_full = altname_full.strip()
                        # Remove the semicolon after the phrase, if any.
                        if altname_full[-1] == ";":
                            altname_full = altname_full[:-1]
                        if altname_full not in protein_information[protein_first_ac]["names"]:
                            protein_information[protein_first_ac]["names"].append(altname_full)
                    # Empty EC
                    result_de_empty_ec = pattern_de_empty_ec.search(line)
                    if result_de_empty_ec:
                        empty_ec = result_de_empty_ec.group(1)
                        # Removing whitespaces before and after the matched phrase.
                        empty_ec = empty_ec.strip()
                        # Remove the semicolon after the phrase, if any.
                        if empty_ec[-1] == ";":
                            empty_ec = empty_ec[:-1]
                        if empty_ec not in protein_information[protein_first_ac]["names"]:
                            protein_information[protein_first_ac]["names"].append(empty_ec)
                    # Empty Short
                    result_de_empty_short = pattern_de_empty_short.search(line)
                    if result_de_empty_short:
                        empty_short = result_de_empty_short.group(1)
                        # Removing whitespaces before and after the matched phrase.
                        empty_short = empty_short.strip()
                        # Remove the semicolon after the phrase, if any.
                        if empty_short[-1] == ";":
                            empty_short = empty_short[:-1]
                        if empty_short not in protein_information[protein_first_ac]["names"]:
                            protein_information[protein_first_ac]["names"].append(empty_short)

                if line[0:2] == "CC":
                    # Find the first line for the protein family/subfamily.
                    # In the case the line is not the header of the similarity section check for the next lines of this section.
                    result_sim = pattern_sim.search(line)
                    if result_sim:
                        protein_sim_firstline = result_sim.group(1)
                        protein_sim_firstline = protein_sim_firstline.lstrip()
                        protein_sim_firstline = protein_sim_firstline.rstrip()
                        protein_sim_firstline = "{}".format(protein_sim_firstline)
                        protein_information[protein_first_ac]["family"] = protein_sim_firstline
                        sim_found = True
                    else:
                        # Find the next lines for the protein family.subfamily, if any.
                        if sim_found:
                            result_sim_cont = pattern_sim_cont.search(line)
                            if result_sim_cont:
                                # Strip any spaces from the start of the string.
                                # Examine whether the leftmost characters include the symbols "-!-" or "---".
                                protein_sim_nextline = result_sim_cont.group(1)
                                protein_sim_nextline = protein_sim_nextline.lstrip()
                                if protein_sim_nextline[:3] == "-!-" or protein_sim_nextline[:3] == "---":
                                    sim_found = False
                                else:
                                    protein_information[protein_first_ac]["family"] = "{} {}".format(protein_information[protein_first_ac]["family"], protein_sim_nextline)
                    if "family" in protein_information[protein_first_ac].keys():
                        if protein_information[protein_first_ac]["family"] == "":
                            print(line)
                            exit()

                if line[0:2] == "DR":
                    # Find the InterPro domains.
                    result_interpro = pattern_interpro.search(line)
                    if result_interpro:
                        protein_interpro_code = result_interpro.group(1)
                        protein_interpro_name = result_interpro.group(2)
                        if protein_interpro_code[-1] == ";":
                            protein_interpro_code = protein_interpro_code[:-1]
                        if protein_interpro_name[-1] == ".":
                            protein_interpro_name = protein_interpro_name[:-1]
                        if "InterPro" not in protein_information[protein_first_ac].keys():
                            protein_information[protein_first_ac]["InterPro"] = [protein_interpro_code]
                        else:
                            protein_information[protein_first_ac]["InterPro"].append(protein_interpro_code)
                    # Find the Pfam domains.
                    result_pfam = pattern_pfam.search(line)
                    if result_pfam:
                        protein_pfam_code = result_pfam.group(1)
                        protein_pfam_name = result_pfam.group(2)
                        protein_pfam_count = result_pfam.group(3)
                        if protein_pfam_code[-1] == ";":
                            protein_pfam_code = protein_pfam_code[:-1]
                        if protein_pfam_name[-1] == ".":
                            protein_pfam_name = protein_pfam_name[:-1]
                        if protein_pfam_count[-1] == ".":
                            protein_pfam_count = protein_pfam_count[:-1]
                        protein_pfam_comb = "{}_{}".format(protein_pfam_code, protein_pfam_count)
                        if "Pfam" not in protein_information[protein_first_ac].keys():
                            protein_information[protein_first_ac]["Pfam"] = [protein_pfam_comb]
                        else:
                            protein_information[protein_first_ac]["Pfam"].append(protein_pfam_comb)
    print("The number of proteins IDs found in the file is: {}".format(pr_counter))
    
    # Remove anything unrelated to the families for the family information of each protein.
    for protein_first_ac in protein_information.keys():
        # If a string for the family has been formed. Filter out anything unrelated to the families.
        if "family" in protein_information[protein_first_ac].keys():
            family_info_pre = protein_information[protein_first_ac]["family"]
            phrase_1 = "Belongs to the "
            phrase_2 = "belongs to the "
            non_family_part = family_info_pre[:15]
            if family_combination_type == 1:
                if (non_family_part == phrase_1) or (non_family_part == phrase_2):
                    protein_information[protein_first_ac]["family"] = family_info_pre[15:]
            elif family_combination_type == 2:
                if (phrase_1 in family_info_pre) or (phrase_2 in family_info_pre):
                    if (phrase_1 in family_info_pre):
                        family_info_pre_splited = family_info_pre.split(phrase_1)
                        family_info_pre_part = family_info_pre_splited[1]
                        protein_information[protein_first_ac]["family"] = family_info_pre_part
                    elif (phrase_2 in family_info_pre):
                        family_info_pre_splited = family_info_pre.split(phrase_2)
                        family_info_pre_part = family_info_pre_splited[1]
                        protein_information[protein_first_ac]["family"] = family_info_pre_part
            else:
                if (non_family_part == phrase_1) or (non_family_part == phrase_2):
                    protein_information[protein_first_ac]["family"] = family_info_pre[15:]
            # If ECO terms are present in the family string.
            family_info = protein_information[protein_first_ac]["family"]
            if "ECO" in family_info:
                # Remove anything after the next to last period.
                family_str = family_info
                family_str_splited = family_str.split(".")
                family_processed = family_str_splited[:-2]
                family_processed = ".".join(family_processed)
                protein_information[protein_first_ac]["family"] = family_processed
                if family_processed == "":
                    print("A family wihtout information (a name) was found. Exiting...")
                    exit()

    # Create a dictionary with keys all the different families.
    unique_families = {}
    for protein_first_ac in protein_information.keys():
        if "family" in protein_information[protein_first_ac].keys():
            pr_family = protein_information[protein_first_ac]["family"]
            if pr_family not in unique_families.keys():
                unique_families[pr_family] = []
    # Add a key for proteins without a protein family.
    unique_families["Unresolved"] = []

    # Group the proteins based on their families.
    # Parse each family name and find all the proteins of that family.
    for pr_fi_ac in protein_information.keys():
        if "family" in protein_information[pr_fi_ac].keys():
            pr_family = protein_information[pr_fi_ac]["family"]
            unique_families[pr_family].append(protein_information[pr_fi_ac])
        else:
            unique_families["Unresolved"].append(protein_information[pr_fi_ac])

    # Sort the protein families alphabetically.
    unique_family_names = list(unique_families.keys())
    unique_family_names.sort(key=str.lower)
    unique_families_sorted = {}
    for ufn in unique_family_names:
        for key in unique_families:
            if ufn == key:
                unique_families_sorted[key] = unique_families[key]

    # For each protein family, find the frequencies of each set of InterPro and Pfam domains related to it based on the
    # proteins assigned to the protein family.
    # unique_families_sorted: {
    # 'GDSL' lipolytic enzyme family: [
    # {'id': 'AAE_RAUSE', 'acs': ['Q3MKY2'], 'family': "'GDSL' lipolytic enzyme family", 'InterPro': ['IPR001087', 'IPR036514', 'IPR035669'], 'Pfam': ['PF00657_1']},
    # {'id': 'AAE_RAUVE', 'acs': ['P86830'], 'family': "'GDSL' lipolytic enzyme family"},
    # {'id': 'ACHE_MAIZE', 'acs': ['B4FZ87', 'A0A3L6E326', 'Q5FC14'], 'family': "'GDSL' lipolytic enzyme family", 'InterPro': ['IPR036915', 'IPR001087', 'IPR036514', 'IPR035669'], 'Pfam': ['PF00657_1']},
    # {'id': 'APG_BRANA', 'acs': ['P40603'], 'family': "'GDSL' lipolytic enzyme family", 'InterPro': ['IPR001087', 'IPR008265', 'IPR036514', 'IPR035669'], 'Pfam': ['PF00657_1']}
    # ...
    # {'id': 'ZWINT_HUMAN', 'acs': ['O95229', 'A6NNV6', 'Q0D2I3', 'Q9BWD0'], 'length': 277, 'InterPro': ['IPR029092'], 'Pfam': ['PF15556_1']}, 
    # {'id': 'ZWINT_RAT', 'acs': ['Q8VIL3', 'Q546Y5'], 'length': 266, 'InterPro': ['IPR029092'], 'Pfam': ['PF15556_1']}, 
    # {'id': 'ZYG1_DICMU', 'acs': ['Q76N59'], 'length': 268}, 
    # {'id': 'ZZZ3_MOUSE', 'acs': ['Q6KAQ7', 'Q3TMK6', 'Q3V189'], 'length': 910, 'InterPro': ['IPR009057', 'IPR017930', 'IPR001005', 'IPR000433', 'IPR043145', 'IPR037830', 'IPR041981'], 'Pfam': ['PF00249_1', 'PF00569_1']}, {'id': 'Z_BPPHM', 'acs': ['Q9G058'], 'length': 64}, 
    # {'id': 'Z_SHEEP', 'acs': ['P08105'], 'length': 79}
    # ]
    # }
    # Connect each family with the set(s) of Pfam domains of the highest count.
    prfam_pfamdoms = {}
    for fam_key in unique_families_sorted.keys():
        count_dom_sets = {}
        for item in unique_families_sorted[fam_key]:
            if "id" not in item.keys():
                print("Protein without ID found in a protein family. Exiting...")
                exit()
            if "acs" not in item.keys():
                print("Protein without ACs found in a protein family. Exiting...")
                exit()
            pr_ac = item["acs"][0]
            if "Pfam" in item.keys():
                pfam_doms = item["Pfam"]
                # Sort the Pfam domains alphabetically.
                pfam_doms.sort(key=str.lower)
                comb_info = ";".join(pfam_doms)
                if comb_info not in count_dom_sets.keys():
                    count_dom_sets[comb_info] = [1, [pr_ac]]
                else:
                    count_dom_sets[comb_info][0] += 1
                    # May want to add a check for duplicate protein IDs (although not such duplicates should exists at
                    # this point).
                    count_dom_sets[comb_info][1].append(pr_ac)
        # If any Pfam domains were found for any protein of the protein family then continue processing.
        if count_dom_sets:
            # Find the set(s) of Pfam domains the highest count for the protein family.
            dom_sets_values = []
            for key_cds in count_dom_sets.keys():
                cds_val = count_dom_sets[key_cds][0]
                dom_sets_values.append(cds_val)
            dom_sets_values = sorted(dom_sets_values, key=int)
            dom_sets_values = dom_sets_values[::-1]
            highest_count = dom_sets_values[0]
            highest_scoring_sets = []
            highest_domscoring_prs = []
            for key_dom in count_dom_sets.keys():
                if highest_count == count_dom_sets[key_dom][0]:
                    highest_scoring_sets.append([key_dom])
                    for pr_ac in count_dom_sets[key_dom][1]:
                        highest_domscoring_prs.append(pr_ac)
            # Assign the set(s) of Pfam domains of the highest count to the protein family.
            prfam_pfamdoms[fam_key] = [highest_count, highest_scoring_sets, highest_domscoring_prs]
            
    # Connect each family with the set(s) of InterPro domains of the highest count.
    prfam_interprodoms = {}
    for fam_key in unique_families_sorted.keys():
        count_dom_sets = {}
        for item in unique_families_sorted[fam_key]:
            if "id" not in item.keys():
                print("Protein without ID found in a protein family. Exiting...")
                exit()
            if "acs" not in item.keys():
                print("Protein without ACs found in a protein family. Exiting...")
                exit()
            pr_ac = item["acs"][0]
            if "InterPro" in item.keys():
                interpro_doms = item["InterPro"]
                # Sort the InterPro domains alphabetically.
                interpro_doms.sort(key=str.lower)
                comb_info = ";".join(interpro_doms)
                if comb_info not in count_dom_sets.keys():
                    count_dom_sets[comb_info] = [1, [pr_ac]]
                else:
                    count_dom_sets[comb_info][0] += 1
                    # May want to add a check for duplicate protein IDs (although not such duplicates should exists at
                    # this point).
                    count_dom_sets[comb_info][1].append(pr_ac)
        # If any InterPro domains were found for any protein of the protein family then continue processing.
        if count_dom_sets:
            # Find the set(s) of InterPro domains the highest count for the protein family.
            dom_sets_values = []
            for key_cds in count_dom_sets.keys():
                cds_val = count_dom_sets[key_cds][0]
                dom_sets_values.append(cds_val)
            dom_sets_values = sorted(dom_sets_values, key=int)
            dom_sets_values = dom_sets_values[::-1]
            highest_count = dom_sets_values[0]
            highest_scoring_sets = []
            highest_domscoring_prs = []
            for key_dom in count_dom_sets.keys():
                if highest_count == count_dom_sets[key_dom][0]:
                    highest_scoring_sets.append([key_dom])
                    for pr_ac in count_dom_sets[key_dom][1]:
                        highest_domscoring_prs.append(pr_ac)
            # Assign the set(s) of InterPro domains of the highest count to the protein family.
            prfam_interprodoms[fam_key] = [highest_count, highest_scoring_sets, highest_domscoring_prs]

    # Find the mean and median length of each protein family.
    prfam_length = {}
    for fam_key in unique_families_sorted.keys():
        pr_set_lengths = []
        fam_unique_acs = []
        for item in unique_families_sorted[fam_key]:
            if "length" in item.keys():
                pr_length = item["length"]
                pr_set_lengths.append(pr_length)
                pr_ac_set = item["acs"]
                for pr_ac in pr_ac_set:
                    # Collect the unique protein IDs for each family.
                    if pr_ac not in fam_unique_acs:
                        fam_unique_acs.append(pr_ac)
        if pr_set_lengths:
            pr_fam_mean_length = mean(pr_set_lengths)
            pr_fam_mean_length = round(pr_fam_mean_length, 2)
            pr_fam_median_length = median(pr_set_lengths)
            pr_fam_median_length = round(pr_fam_median_length, 2)
            prfam_length[fam_key] = [pr_fam_mean_length, pr_fam_median_length, fam_unique_acs]

    # Connect each family with the highest occuring names.
    fam_names_dict = {}
    # For each family...
    for fam_key in unique_families_sorted.keys():
        fam_names_dict[fam_key] = {}
        # For each protein of the fmaily...
        for pr_dict in unique_families_sorted[fam_key]:
            if "names" in pr_dict.keys():
                for pr_name in pr_dict["names"]:
                    if pr_name not in fam_names_dict[fam_key].keys():
                        fam_names_dict[fam_key][pr_name] = 1
                    else:
                        fam_names_dict[fam_key][pr_name] += 1
    # Sorting the names based on their frequency.
    # The consecutive symbols "::" do not exist in any of the protein names of any protein for any of the protein families.
    fam_names_dict_sorted = {}
    for key in fam_names_dict.keys():
        fam_names_dict_sorted[key] = {}
        name_freqs = list(fam_names_dict[key].values())
        name_freqs = list(dict.fromkeys(name_freqs))
        name_freqs = sorted(name_freqs, key=int)
        name_freqs = name_freqs[::-1]
        for nf in name_freqs:
            for key_name in fam_names_dict[key].keys():
                if fam_names_dict[key][key_name] == nf:
                    fam_names_dict_sorted[key][key_name] = copy.deepcopy(fam_names_dict[key][key_name])
        fam_names_dict[key] = copy.deepcopy(fam_names_dict_sorted[key])

    # Write the information to TXT files.
    # Protein families, numbered.
    info_file_families = open("prfamilies_numbered.tsv", "w")
    info_file_families.write("Family Number\tProtein Family\tProtein Names\n")
    prf_num = 0
    for key in prfam_pfamdoms.keys():
        info_file_families.write("{}\t{}".format(prf_num, key))
        if key in fam_names_dict_sorted.keys():
            for key_name in fam_names_dict_sorted[key]:
                info_file_families.write("\t{}::{}".format(key_name, fam_names_dict_sorted[key][key_name]))
        else:
            info_file_families.write("\tNone")
        info_file_families.write("\n")
        prf_num += 1
    info_file_families.close()
    # Pfam domains
    info_file_pfam_domains = open("prfamilies_pfamdomains.tsv", "w")
    info_file_pfam_domains.write("Protein Family\tPfam Domain Set\tFrequency\n")
    for key in prfam_pfamdoms.keys():
        dom_sets = []
        for cur_dom_set in prfam_pfamdoms[key][1]:
            dom_sets.append(cur_dom_set)
        pr_acs = prfam_pfamdoms[key][2]
        info_file_pfam_domains.write("{}\t{}\t{}\t{}\n".format(key, dom_sets, prfam_pfamdoms[key][0], pr_acs))
    info_file_pfam_domains.close()
    # InterPro domains
    info_file_interpro_domains = open("prfamilies_interprodomains.tsv", "w")
    info_file_interpro_domains.write("Protein Family\tInterPro Domain Set\tFrequency\n")
    for key in prfam_interprodoms.keys():
        dom_sets = []
        for cur_dom_set in prfam_interprodoms[key][1]:
            dom_sets.append(cur_dom_set)
        pr_acs = prfam_interprodoms[key][2]
        info_file_interpro_domains.write("{}\t{}\t{}\t{}\n".format(key, dom_sets, prfam_interprodoms[key][0], pr_acs))
    info_file_interpro_domains.close()
    # Lengths based on the families by Uniprot.
    info_file_lengths = open("prfamilies_length.tsv", "w")
    info_file_lengths.write("Protein Family\tMean Length\tMedian Length\tProtein IDs\n")
    for key in prfam_length.keys():
        info_file_lengths.write("{}\t{}\t{}\t{}\n".format(key, prfam_length[key][0], prfam_length[key][1], prfam_length[key][2]))
    info_file_lengths.close()

    # Write the information to CSV files.
    # Pfam domains
    with open("prfamilies_pfamdomains.csv", 'w', newline='') as info_file_pfam_domains_csv:
        writer = csv.writer(info_file_pfam_domains_csv)
        dom_row = ["Protein Family", "Pfam Domain Set(s)", "Frequency", "Protein ACs"]
        writer.writerow(dom_row)
        for key in prfam_pfamdoms.keys():
            dom_sets = []
            for cur_dom_set in prfam_pfamdoms[key][1]:
                dom_sets.append(cur_dom_set)
            pr_acs = prfam_pfamdoms[key][2]
            dom_row = [key, dom_sets, prfam_pfamdoms[key][0], pr_acs]
            writer.writerow(dom_row)
    # InterPro domains
    with open("prfamilies_interprodomains.csv", 'w', newline='') as info_file_interpro_domains_csv:
        writer = csv.writer(info_file_interpro_domains_csv)
        dom_row = ["Protein Family", "InterPro Domain Set(s)", "Frequency", "Protein ACs"]
        writer.writerow(dom_row)
        for key in prfam_interprodoms.keys():
            dom_sets = []
            for cur_dom_set in prfam_interprodoms[key][1]:
                dom_sets.append(cur_dom_set)
            pr_acs = prfam_interprodoms[key][2]
            dom_row = [key, dom_sets, prfam_interprodoms[key][0], pr_acs]
            writer.writerow(dom_row)
    # Lengths
    with open("prfamilies_length.csv", 'w', newline='') as info_file_lengths_csv:
        writer = csv.writer(info_file_lengths_csv)
        len_row = ["Protein Family", "Mean Length", "Median Length"]
        writer.writerow(len_row)
        for key in prfam_length.keys():
            len_row = [key, prfam_length[key][0], prfam_length[key][1], prfam_length[key][2]]
            writer.writerow(len_row)


if __name__ == "__main__":
    arg_uniprot_falt_file = "uniprot_sprot.dat"
    arg_family_combination_type = 1
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input":
                arg_uniprot_falt_file = sys.argv[i+1]
            elif sys.argv[i] == "-c" or sys.argv[i] == "--combination":
                arg_family_combination_type = int(sys.argv[i+1])
            elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
                help_message()
                exit()
    unidat(arg_uniprot_falt_file, arg_family_combination_type)

