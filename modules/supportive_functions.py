import os
import time
import shutil
# ProteoSeeker modules
import command_process


def check_i_value(i_value, option_str):
    if i_value in ["False", "false", "f", "0"]:
        i_value = False
    elif i_value in ["True", "true", "t", "1"]:
        i_value = True
    else:
        print("Wrong value given for the '{}' option. Exiting.".format(option_str))
        exit()
    return i_value


def read_file(file_path):
    file_handle = open(file_path, "r")
    pre_lines = file_handle.readlines()
    lines = []
    for i in pre_lines:
        i = i.rstrip("\n")
        lines.append(i)
    file_handle.close()
    return lines


def end_time_analysis(label, start_time, output_log_file):
    # End time
    end_time = time.time()
    elpased_time = end_time - start_time
    elpased_time_min = elpased_time / 60
    elpased_time_hour = elpased_time_min / 60
    elpased_time = round(elpased_time, 3)
    elpased_time_min = round(elpased_time_min, 3)
    elpased_time_hour = round(elpased_time_hour, 3)
    print("\n{}".format(label))
    print("Time elapsed: {} seconds / {} minutes / {} hours".format(elpased_time, elpased_time_min, elpased_time_hour))
    output_log_file.write("{}\n".format(label))
    output_log_file.write("Time elapsed: {} seconds / {} minutes / {} hours\n".format(elpased_time, elpased_time_min, elpased_time_hour))
    output_log_file.write("\n{}\n\n\n".format(100*"-"))
    return elpased_time


def write_time_dict(time_dict, time_analyis_path):
    time_analysis_file = open(time_analyis_path, "w")
    for key in time_dict.keys():
        temp_line = "{}\t{}".format(key, time_dict[key])
        time_analysis_file.write("{}\n".format(temp_line))
    time_analysis_file.close()


def reduce_volume(clear_space, input_folder, output_path_trimmed, input_log_file, del_pick):
    if clear_space:
        input_log_file.write("Deleting files:\n")
        if 1 in del_pick:
            input_file_paths = os.listdir(input_folder)
            for i in input_file_paths:
                input_file_rel_path = "{}/{}".format(input_folder, i)
                print("Deleting file: {}".format(input_file_rel_path))
                os.remove(input_file_rel_path)
                input_log_file.write("File deleted: {}\n".format(input_file_rel_path))
        if 2 in del_pick:
            # The file paths in "ca_file_paths" contain the file paths from the "trimmed_results" folder which are compressed, which are the files that
            # contain the "fastq" name in their names. For being sure we also check that the suffix is ".gz" before deleting the compressed files.
            trimmed_file_paths = os.listdir(output_path_trimmed)
            for i in trimmed_file_paths:
                if "fastq" in i and i[-3:] == ".gz":
                    local_file_path = "{}/{}".format(output_path_trimmed, i)
                    print("Deleting file: {}".format(local_file_path))
                    os.remove(local_file_path)
                    input_log_file.write("File deleted: {}\n".format(local_file_path))
        input_log_file.write("\n")


def unzip_files(file_paths, unzip_path, clear, input_log_file, output_log_file, gzip_path):
    if unzip_path:
        if clear:
            if os.path.exists(unzip_path):
                shutil.rmtree(unzip_path)
        os.mkdir(unzip_path)
    title_1 = "Uncompressing files - gzip:"
    title_2 = "Version of gzip:"
    capture_status = True
    shell_status = True
    pr_status = False
    for i in file_paths:
        file_name = i.split("/")[-1]
        file_name_parts = file_name.split(".")
        file_name = ".".join(file_name_parts[:-1])
        if unzip_path is not None:
            if gzip_path:
                phrase_1 = "\"{}\" -dk -c \"{}\"".format(gzip_path, i)
                phrase_1 = "{} > \"{}{}\"".format(phrase_1, unzip_path, file_name)
                phrase_2 = "{} --version".format(gzip_path)
            else:
                phrase_1 = "gzip -dk -c \"{}\"".format(i)
                phrase_1 = "{} > \"{}{}\"".format(phrase_1, unzip_path, file_name)
                phrase_2 = "gzip --version"
        else:
            if gzip_path:
                phrase_1 = "\"{}\" -dk \"{}\"".format(gzip_path, i)
                phrase_2 = "\"{}\" --version".format(gzip_path)
            else:
                phrase_1 = "gzip -dk \"{}\"".format(i)
                phrase_2 = "gzip --version"
        command_process.command_run(phrase_1, phrase_2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)


def combine_predictions(hmmer_enz_domains_all_proteins, blastp_no_doms_below_threshold, proteins_combined_file):
    print("\nCombining predictions...")
    if os.path.exists(hmmer_enz_domains_all_proteins) and os.path.exists(blastp_no_doms_below_threshold):
        filenames = [hmmer_enz_domains_all_proteins, blastp_no_doms_below_threshold]
    elif os.path.exists(hmmer_enz_domains_all_proteins):
        filenames = [hmmer_enz_domains_all_proteins]
    elif os.path.exists(blastp_no_doms_below_threshold):
        filenames = [blastp_no_doms_below_threshold]
    else:
        filenames = []
    with open(proteins_combined_file, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def tax_info_convert(bin_group_tax_dict):
    print("\nReforming the information related to the binned proteins...")
    # The process is faster if the key is the protein and the value is its bin, species and caterogization method.
    pr_tax_info_dcit = {}
    if bin_group_tax_dict:
        for bin_id in bin_group_tax_dict.keys():
            cat_type = "U"
            taxonomies = bin_group_tax_dict[bin_id]["species"]
            for key_info in bin_group_tax_dict[bin_id].keys():
                if key_info == "protein_ids_d":
                    cat_type = "D"
                if key_info == "protein_ids_i":
                    cat_type = "I"
                if key_info == "protein_ids_k":
                    cat_type = "K"
                if key_info in ["protein_ids_d", "protein_ids_i", "protein_ids_k"]:
                    for key_pr_ac in bin_group_tax_dict[bin_id][key_info]:
                        pr_tax_info_dcit[key_pr_ac] = [bin_id, taxonomies, cat_type]
    return pr_tax_info_dcit