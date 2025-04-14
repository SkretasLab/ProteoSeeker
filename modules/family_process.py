import os
import re
import copy
# ProteoSeeker modules
import supportive_functions


def find_fam_names(family_code, family_group_name, protein_db_path, fpd_gen_folder, fam_nums_file_name, name_thr, input_protein_names_status, input_protein_names, pr_names_dict, conda_sh_path):
    # Creation of the pnr database if the user selected the process to be performed and also provided a path for the nr database.
    pr_names = []
    fpd_name = ""
    fpd_fasta = ""
    fpd_folder = None
    fpd_fasta_prefix = ""

    # The name of the pnr database.
    # To create the pnr database the path to the nr database is needed.
    if os.path.exists(protein_db_path):
        fpd_folder = "{}/fpd_{}".format(fpd_gen_folder, family_group_name)
        fpd_fasta = "{}/fpd_{}.fasta".format(fpd_folder, family_group_name)
        fpd_fasta_prefix = "{}/fpd_{}".format(fpd_folder, family_group_name)
    else:
        print("\nThe path for the nr database has not been set or is wrong. Exiting.")
        exit()
    fpd_name = "{}/fpd_{}_database".format(fpd_folder, family_group_name)

    # If the pnr database already exists based on the selected family, do not create it again. It is already checked that, if no family code was
    # selected the pnr database does not exist based on the random counter.
    if not os.path.exists(fpd_folder):
        # Dictionary for protein family codes and protein family names.
        pr_code_famname_dict = {}
        pr_code_prsname_dict = {}
        fam_nums_lines = supportive_functions.read_file(fam_nums_file_name)
        for line in fam_nums_lines:
            line_splited = line.split("\t")
            fam_num = line_splited[0]
            fam_name = line_splited[1]
            pr_names = line_splited[2:]
            pr_code_famname_dict[fam_num] = fam_name
            pr_code_prsname_dict[fam_num] = pr_names
        # Create the list of initial enzyme names.
        if family_code is not None:
            pre_pr_names = []
            if pr_code_famname_dict.keys():
                for family_code_cur in family_code:
                    pre_pr_names_ori = pr_code_famname_dict[family_code_cur]
                    # For each family or subfamily case keep a name same as the the original text and a name without the "family" or "subfamily" word.
                    # AAA ATPase family. Katanin p60 subunit T.A1 subfamily. A-like 1 sub-subfamily
                    # pre_pr_names_ori_splited = ["AAA ATPase ", " Katanin p60 subunit T.A1 ", " A-like 1 sub-subfamily"]
                    # pre_pr_names_fam_words = ["AAA ATPase", "Katanin p60 subunit T.A1", "A-like 1 sub-subfamily"]
                    # pre_pr_names = pre_pr_names_fam_words
                    if ("family." in pre_pr_names_ori) or "subfamily." in pre_pr_names_ori:
                        pre_pr_names_fam_words = []
                        pre_pr_names_ori_splited = re.split('family.|subfamily.', pre_pr_names_ori)
                        for item_ori in pre_pr_names_ori_splited:
                            item_ori = item_ori.strip()
                            pre_pr_names_fam_words.append(item_ori)
                        pre_pr_names = copy.deepcopy(pre_pr_names_fam_words)
                        for item in pre_pr_names_fam_words:
                            if item[-10:] == " subfamily":
                                item = item[:-10]
                            if item[-7:] == " family":
                                item = item[:-7]
                            item = item.strip()
                            pre_pr_names.append(item)
                    else:
                        # Golgi pH regulator (TC 1.A.38) family
                        # pre_pr_names_ori = "Golgi pH regulator (TC 1.A.38) family"
                        # pre_pr_names = ["Golgi pH regulator (TC 1.A.38) family"]
                        # if pre_pr_names_ori[-10:] == " subfamily": False
                        # pre_pr_names_ori[-7:] == " family": True
                        # item = "Golgi pH regulator (TC 1.A.38)"
                        # pre_pr_names = ["Golgi pH regulator (TC 1.A.38) family", "Golgi pH regulator (TC 1.A.38)"]
                        pre_pr_names = [pre_pr_names_ori]
                        if pre_pr_names_ori[-10:] == " subfamily":
                            item = pre_pr_names_ori[:-10]
                        if pre_pr_names_ori[-7:] == " family":
                            item = pre_pr_names_ori[:-7]
                        item = item.strip()
                        pre_pr_names.append(item)
            for family_code_cur in family_code:
                # Add each protein name based on a threhsold of frequency.
                # Find the highest occuring name.
                if name_thr <= 0:
                    print("\nThe threshold for the naming selection of the protein family can not be 0 or negative. Exiting.")
                    exit()
                if isinstance(name_thr, float):
                    name_max_splited =  pr_code_prsname_dict[family_code_cur][0].split("::")
                    name_max = int(name_max_splited[1])
                    sel_thr = name_max * name_thr
                elif isinstance(name_thr, int):
                    sel_thr = float(name_thr)
                else:
                    print("\nAn error occured while computing the protein naming threshold. Exiting.")
                    exit()
                print("Protein family: {} - Protein name frequency threshold: {}".format(family_code_cur, sel_thr))
                for prn_cur in  pr_code_prsname_dict[family_code_cur]:
                    prn_cur_splited = prn_cur.split("::")
                    prn_cur_name = prn_cur_splited[0]
                    prn_cur_num = float(prn_cur_splited[1])
                    if prn_cur_num >= sel_thr:
                        pre_pr_names.append(prn_cur_name)
            # If selected the input enzyme names overwritten the ones found automatically (if any were found).
            # If any enzyme name have been given as input, they are added in the list of names for the selected enzyme
            # code.
            if not input_protein_names_status:
                if input_protein_names:
                    pr_names = pre_pr_names + input_protein_names
                else:
                    pr_names = copy.deepcopy(pre_pr_names)
            else:
                pr_names = copy.deepcopy(input_protein_names)
            # Deduplicate list of enzyme names.
            pr_names_dupl = copy.deepcopy(pr_names)
            pr_names = []
            for item in pr_names_dupl:
                if item not in pr_names:
                    pr_names.append(item)
    else:
        print("\nThe path to the fpr database ({}) for the selected protein family(ies) already exists. Skipping creation of the fpr database.".format(fpd_folder))

    temp_list = [pr_names, fpd_fasta, fpd_folder, fpd_name, fpd_fasta_prefix]
    if pr_names_dict:
        pr_names_dict[1] = copy.deepcopy(temp_list)
    else:
        pr_names_dict[0] = copy.deepcopy(temp_list)
    return fpd_fasta, fpd_name, fpd_folder, pr_names_dict