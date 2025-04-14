import os
# ProteoSeeker modules
import supportive_functions


def file_reads(output_path, paired_end, file_paths, file_paths_p, file_read_info_path):
    print("\nAnalyzing FastQC results...")
    # Number of total sequences from all input files
    total_sequences_num = 0
    paired_sequences_num = 0
    single_sequences_num = 0
    # The paired sequence refer to the sum of the sequences of only 1 of the 2 files of each pair of files.
    paired_to_seqnum_dict = {}
    single_to_sequences_dict = {}
    fastqc_files = os.listdir(output_path)
    for fqc_file in fastqc_files:
        if fqc_file[-5:] == ".html":
            # Counting the total number of sequences based on all input files.
            if paired_end:
                for key_p_temp in file_paths_p.keys():
                    file_path_p_temp_1 = file_paths_p[key_p_temp][0]
                    file_path_p_temp_2 = file_paths_p[key_p_temp][1]
                    # Matching the files from the fastqc folder to the ones provided.
                    # Catch the name after the last slash, if any exists.
                    if "/" in file_path_p_temp_1:
                        file_path_p_temp_1_splited = file_path_p_temp_1.split("/")
                        file_path_p_temp_1 = file_path_p_temp_1_splited[-1]
                    if "/" in file_path_p_temp_2:
                        file_path_p_temp_2_splited = file_path_p_temp_2.split("/")
                        file_path_p_temp_2 = file_path_p_temp_2_splited[-1]
                    # Remove the suffix of ".fastq" or ".fast.gz", if any exists.
                    if file_path_p_temp_1[-6:] == ".fastq":
                        file_path_p_temp_1 = file_path_p_temp_1[:-6]
                    elif file_path_p_temp_1[-9:] == ".fastq.gz":
                        file_path_p_temp_1 = file_path_p_temp_1[:-9]
                    if file_path_p_temp_2[-6:] == ".fastq":
                        file_path_p_temp_2 = file_path_p_temp_2[:-6]
                    elif file_path_p_temp_2[-9:] == ".fastq.gz":
                        file_path_p_temp_2 = file_path_p_temp_2[:-9]
                    # The path of the first parsed from the provided files.
                    fqc_file_base_name = fqc_file.split("_fastqc.html")[0]
                    if (fqc_file_base_name == file_path_p_temp_1) or (fqc_file_base_name == file_path_p_temp_2):
                        fqc_temp_file_path = "{}/{}".format(output_path, fqc_file)
                        fqc_lines_temp = supportive_functions.read_file(fqc_temp_file_path)
                        for line_fqc_temp in fqc_lines_temp:
                            if "Total Sequences</td><td>" in line_fqc_temp:
                                splited_fqc_temp_1 = line_fqc_temp.split("Total Sequences</td><td>")
                                splited_fqc_temp_2 = splited_fqc_temp_1[1].split("</td></tr>")
                                sequences_num = int(splited_fqc_temp_2[0])
                                total_sequences_num += sequences_num
                                if key_p_temp not in paired_to_seqnum_dict.keys():
                                    paired_to_seqnum_dict[key_p_temp] = sequences_num
            else:
                for temp_i in range(0, len(file_paths)):
                    if temp_i not in single_to_sequences_dict.keys():
                        file_path_temp = file_paths[temp_i]
                        fqc_file_base_name = fqc_file.split("_fastqc.html")[0]
                        if (fqc_file_base_name in file_path_temp):
                            fqc_temp_file_path = "{}/{}".format(output_path, fqc_file)
                            fqc_lines_temp = supportive_functions.read_file(fqc_temp_file_path)
                            for line_fqc_temp in fqc_lines_temp:
                                if "Total Sequences</td><td>" in line_fqc_temp:
                                    splited_fqc_temp_1 = line_fqc_temp.split("Total Sequences</td><td>")
                                    splited_fqc_temp_2 = splited_fqc_temp_1[1].split("</td></tr>")
                                    sequences_num = int(splited_fqc_temp_2[0])
                                    total_sequences_num += sequences_num
                                    single_to_sequences_dict[temp_i] = sequences_num
    if paired_end:
        paired_sequences_num = 0
        for key in paired_to_seqnum_dict.keys():
            paired_sequences_num += paired_to_seqnum_dict[key]
    else:
        single_sequences_num = total_sequences_num
    # Pritn information
    print("\nThe input number of reads is: {}".format(total_sequences_num))
    if paired_end:
        print("The input number of paired-end reads, of either direction (forward or reverse) is: {}".format(paired_sequences_num))
    else:
        print("The input number of single-end reads is: {}".format(single_sequences_num))
    # Write information
    if not os.path.exists(file_read_info_path):
        file_read_info_file = open(file_read_info_path, "w")
        file_read_info_file.write("Input reads\t{}\n".format(total_sequences_num))
        if paired_end:
            file_read_info_file.write("Input paired-end reads\t{}\n".format(paired_sequences_num))
        else:
            file_read_info_file.write("Input single-end reads\t{}\n".format(single_sequences_num))
        file_read_info_file.close()
    else:
        file_read_info_file = open(file_read_info_path, "a+")
        file_read_info_file.write("Trimmed reads\t{}\n".format(total_sequences_num))
        if paired_end:
            file_read_info_file.write("Trimmed paired-end reads\t{}\n".format(paired_sequences_num))
        else:
            file_read_info_file.write("Trimmed single-end reads\t{}\n".format(single_sequences_num))
        file_read_info_file.close()