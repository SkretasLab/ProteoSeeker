import os
# ProteoSeeker modules
import supportive_functions


def pr_fam_pred(analysis_fam_names, blastp_info_swissprot_file, proteins_combined_file, pr_fams_path, family_info_path):
    print("\nPredicting the protein family of each protein based on its hit of lowest e-value against the Swiss-Prot protein database...")
    # Keeping only the protein with the lowest e-value from the hits against the Swiss-Prot database.
    dict_sp_protac = {}
    protein_acc = None
    min_evalue = None
    header = True
    with open(blastp_info_swissprot_file) as lines_sp:
        for line_sp in lines_sp:
            if header:
                header = False
                continue
            line_sp = line_sp.rstrip("\n")
            splited_sp = line_sp.split("\t")
            protein_acc = splited_sp[0]
            swissprot_prot = splited_sp[1]
            e_value = float(splited_sp[3])
            if "." in swissprot_prot:
                swissprot_prot_splited = swissprot_prot.split(".")
                swissprot_prot = swissprot_prot_splited[0]
            if protein_acc not in dict_sp_protac.keys():
                dict_sp_protac[protein_acc] = swissprot_prot
                min_evalue = e_value
            else:
                if e_value < min_evalue:
                    dict_sp_protac[protein_acc] = swissprot_prot
                    min_evalue = e_value
    # Protein sequences
    dict_seqs = {}
    if os.path.exists(proteins_combined_file):
        with open(proteins_combined_file) as lines_seqs:
            for line_pr_seq in lines_seqs:
                line_pr_seq = line_pr_seq.rstrip("\n")
                if line_pr_seq[0] == ">":
                    protein_acc = line_pr_seq[1:]
                    if protein_acc in dict_seqs.keys():
                        print("Error. A duplicate header was found. Exiting.")
                        exit()
                    else:
                        dict_seqs[protein_acc] = ""
                else:
                    dict_seqs[protein_acc] += line_pr_seq
    # Protein families from the Swiss-Prot/UniprotKB database.
    fams_lenth_prs_dict = {}
    fam_lines = supportive_functions.read_file(pr_fams_path)
    header = True
    for line in fam_lines:
        if header:
            header = False
            continue
        else:
            line_splited = line.split("\t")
            pr_fam_name = line_splited[0]
            swissprot_proteins = line_splited[-1]
            swissprot_proteins = swissprot_proteins[2:-2]
            swissprot_proteins_splited = swissprot_proteins.split("', '")
            temp_list = line_splited[1:-1] + [swissprot_proteins_splited]
            fams_lenth_prs_dict[pr_fam_name] = temp_list
    # For each putative protein being alazyed (key):
    # Based on its best hit (swissprot_pr) against the Swiss-Prot database:
    # For each protein family (key_fam) collected from the Swiss-Prot database, take its proteins (fam_prs) and:
    # If the best hit protein (swissprot_pr) is part of fam_prs:
    # And the length of the protein family (key_fam) is not None:
    # Add the putative protein (key) in a dictionary as a key and the family of its best hit as it value in a list
    # Add the length information of the faimly of the best hit (predicted protein family)
    # In the end, if no protein families had been initially selected to seek add "-", otherwise
    # After checking if any of the selected proteins family to seek matches the predicted one, if at least one matches add "1", otherwise if none matches
    # add "0". The same procedure is applied for any protein family the swissprot_pr is part of.
    swiss_fams_len_comp_dict = {}
    for key in dict_sp_protac.keys():
        swissprot_pr = dict_sp_protac[key]
        pr_len = None
        if key in dict_seqs.keys():
            pr_len = float(len(dict_seqs[key]))
        for key_fam in fams_lenth_prs_dict.keys():
            if fams_lenth_prs_dict[key_fam][2]:
                fam_prs = fams_lenth_prs_dict[key_fam][2]
                if swissprot_pr in fam_prs:
                    if pr_len is not None:
                        if key not in swiss_fams_len_comp_dict.keys():
                            swiss_fams_len_comp_dict[key] = [[key_fam]]
                        else:
                            swiss_fams_len_comp_dict[key].append([key_fam])
                        fam_mean_len = float(fams_lenth_prs_dict[key_fam][0])
                        mean_len_dist = abs(pr_len - fam_mean_len)
                        mean_len_dist = round(mean_len_dist, 2)
                        relative_change = (mean_len_dist / fam_mean_len) * 100
                        relative_change = round(relative_change, 2)
                        swiss_fams_len_comp_dict[key][-1].append(fam_mean_len)
                        swiss_fams_len_comp_dict[key][-1].append(mean_len_dist)
                        swiss_fams_len_comp_dict[key][-1].append(relative_change)
                        if not analysis_fam_names:
                            swiss_fams_len_comp_dict[key][-1].append("-")
                        else:
                            swiss_fams_len_comp_dict[key][-1].append(0)
                            for sel_fam in analysis_fam_names:
                                if sel_fam == key_fam:
                                    swiss_fams_len_comp_dict[key][-1].append(1)
                                    break
    # Write the information in a file.
    new_file = open(family_info_path, "w")
    for key_comp in swiss_fams_len_comp_dict.keys():
        for set_info in swiss_fams_len_comp_dict[key_comp]:
            new_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(key_comp, set_info[0], set_info[1], set_info[2], set_info[3], set_info[4]))
    new_file.close()
    return swiss_fams_len_comp_dict, dict_seqs