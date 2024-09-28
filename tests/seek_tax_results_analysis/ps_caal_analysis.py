import os
import sys


def check_i_value(i_value, option_str):
    if i_value in ["False", "false", "f", "0"]:
        i_value = False
    elif i_value in ["True", "true", "t", "1"]:
        i_value = True
    else:
        print("Wrong value given for the '{}' option. Exiting.".format(option_str))
        exit()
    return i_value


def read_file(filename):
    file_handle = open(filename, "r")
    lines_raw = file_handle.readlines()
    lines = []
    for line in lines_raw:
        line = line.rstrip("\n")
        lines.append(line)
    file_handle.close()
    return lines


def caalstats(ps_results, blastp_info_path):
    print("Analyzing the results by the CA and ALs runs...\n")

    if (not ps_results) or (ps_results == ""):
        print("The path to the input directory of the results by the CAs and ALs runs was not found. Exiting.")
        exit()

    if (not blastp_info_path) or (blastp_info_path == ""):
        print("The path to the input BLASTP information file was not found. Exiting.")
        exit()

    # Collect the information from the BLASTP file.
    # Structure:
    # drr163688:
    #   k141_9354_2378_2884_+: Pyrobaculum aerophilum
    #   k127_3944_212_385_+: Sphingobium yanoikuyae
    # srr3961740
    #   k141_7521_20443_21138_-: Novosphingobium sp.
    # ...
    pr_blastp_nr_tax_dict = {}
    pr_seq_label_dict = {}
    protein_sequence_set = set()
    blastp_lines = read_file(blastp_info_path)
    first_line = True
    for line in blastp_lines:
        if first_line:
            first_line = False
            continue
        line_splited = line.split("\t")
        protein_sequence = line_splited[0]
        blastp_nr_taxon = line_splited[1]
        if (len(line_splited) >= 3):
            pr_label = line_splited[2]
            pr_blastp_nr_tax_dict[pr_label] = blastp_nr_taxon
        else:
            pr_label = None
            pr_blastp_nr_tax_dict[protein_sequence] = blastp_nr_taxon
        protein_sequence_set.add(protein_sequence)
        if pr_label is not None:
            pr_seq_label_dict[protein_sequence] = pr_label

    # Collect the taxonomy predicted by ProteoSeeker for each protein of each method of each sample.
    # Each method for a sample contains the same protein IDs for the same proteins.
    sample_method_dict = {}
    st_file_names = os.listdir(ps_results)
    for sfn in st_file_names:
        print("Analyzing file: {}".format(sfn))
        sfn_splited = sfn.split("_")
        sample_name = sfn_splited[2]
        sfn_len = len(sfn_splited)
        if sfn_len == 4:
            method = sfn_splited[3]
        elif sfn_len == 5:
            method = "{}_{}".format(sfn_splited[3], sfn_splited[4])
        else:
            print("Unknown sample/method combination was found. Exiting.")
            exit()
        # Add the sample name as a key to the dictionary.
        if sample_name not in sample_method_dict.keys():
            sample_method_dict[sample_name] = {}
        # Ass the method as a key to the level of the sample.
        sample_method_dict[sample_name][method] = {}
        # Get the full path of the annotation file.
        sfn_full_path = "{}/{}".format(ps_results, sfn)
        sfn_annotation_path = "{}/annotation_results/annotation_info.txt".format(sfn_full_path)
        found_protein = False
        pr_seq_line = ""
        sfn_annotation_lines = read_file(sfn_annotation_path)
        for i in range(0, len(sfn_annotation_lines)):
            line = sfn_annotation_lines[i]
            if line == "Protein sequence:":
                for j in range(i+1, len(sfn_annotation_lines)):
                    next_line = sfn_annotation_lines[j]
                    if (not next_line) or (next_line == ""):
                        # Check if the protein sequence matches any of the given blastp information.
                        if pr_seq_line in protein_sequence_set:
                            taxon_protein = pr_seq_line
                            found_protein = True
                        # Empty the sequence.
                        pr_seq_line = ""
                        break
                    else:
                        pr_seq_line = "{}{}".format(pr_seq_line, next_line)
            if line[0:2] == "//":
                found_protein = False
            if found_protein:
                if "Protein taxonomy" in line:
                    tax_line = sfn_annotation_lines[i+1]
                    taxl_line_splited = tax_line.split("||")
                    pred_taxon = taxl_line_splited[0]
                    if taxon_protein in pr_seq_label_dict:
                        pr_label = pr_seq_label_dict[taxon_protein]
                        sample_method_dict[sample_name][method][pr_label] = pred_taxon
                    else:
                        sample_method_dict[sample_name][method][pr_seq_line] = pred_taxon

    # Creating a sorted dictionary.
    method_order = ["k8", "k16", "k72", "cnr", "mnr"]
    sample_method_dict_sorted = {}
    for key_sample in sample_method_dict.keys():
        sample_method_dict_sorted[key_sample] = {}
        for mo_item in method_order:
            for key_method in sample_method_dict[key_sample].keys():
                if mo_item == key_method:
                    method_value = sample_method_dict[key_sample][key_method]
                    if mo_item == "k72":
                        mo_item_mod = "k77"
                    else:
                        mo_item_mod = mo_item
                    sample_method_dict_sorted[key_sample][mo_item_mod] = method_value

    # Writing the results in a file.
    # If no taxonomy information was found for a protein in the expected sample that means that no taxon was associated with the protein
    # or the taxon was filtered out after the filtering threshold for the Kraken2 taxonomy route.
    blastp_taxa_comp_path = "results_analysis.tsv"
    blastp_taxa_comp_file = open(blastp_taxa_comp_path, "w")
    blastp_taxa_comp_file.write("sample\tmethod\tprotein\tblastp_nr_species\tpredicted_species\n")
    for key_sample in sample_method_dict_sorted.keys():
        for key_method in sample_method_dict_sorted[key_sample].keys():
            if sample_method_dict_sorted[key_sample][key_method]:
                for pr_lid in sample_method_dict_sorted[key_sample][key_method].keys():
                    blastp_nr_taxa = pr_blastp_nr_tax_dict[pr_lid]
                    pred_taxa_temp = sample_method_dict_sorted[key_sample][key_method][pr_lid]
                    pred_taxa_temp = pred_taxa_temp.split("\t")
                    pred_taxa = ",".join(pred_taxa_temp)
                    blastp_taxa_comp_file.write("{}\t{}\t{}\t{}\t{}\n".format(key_sample, key_method, pr_lid, blastp_nr_taxa, pred_taxa))
    blastp_taxa_comp_file.close()


if __name__ == "__main__":
    arg_ps_results = ""
    arg_blastp_info_path = ""
    if len(sys.argv) > 1:
        arg_input_command = "{}".format(sys.argv[0])
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-p" or sys.argv[i] == "--proteoseeker-results-path":
                arg_ps_results = sys.argv[i+1]
            if sys.argv[i] == "-b" or sys.argv[i] == "--blastp-path":
                arg_blastp_info_path = sys.argv[i+1]
    caalstats(arg_ps_results, arg_blastp_info_path)
