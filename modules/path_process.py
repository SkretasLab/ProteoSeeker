import os


def paired_end_file_paths(folder, ignore_comp=True):
    file_paths_p = {}
    passed = []
    r = 0
    passed_counter = 0
    pre_file_paths = os.listdir(folder)
    for i in range(0, len(pre_file_paths)):
        if "Summary" in pre_file_paths[i] or "Error" in pre_file_paths[i]:
            continue
        if ignore_comp:
            if (pre_file_paths[i][-3:] == ".gz") or (pre_file_paths[i][-3:] == ".sh") or (pre_file_paths[i][-4:] == ".txt"):
                continue
        passed_counter += 1
        if i in passed:
            continue
        # The last file checked should always have been found earlier as the pair of another file. Thus it should always be in the list passed. But in case a mistake has been made in the names of the files
        # then it may not be in the list passed and the counter will be increased thus in the end an error will occur based on the difference of the length of the list passed and of the counter.
        if i + 1 == len(pre_file_paths):
            break
        for j in range(i + 1, len(pre_file_paths)):
            if len(pre_file_paths[i]) != len(pre_file_paths[j]):
                continue
            dif = 0
            dif_chars = []
            for char_index in range(0, len(pre_file_paths[i])):
                if pre_file_paths[i][char_index] != pre_file_paths[j][char_index]:
                    dif += 1
                    dif_chars.append(pre_file_paths[i][char_index])
                    dif_chars.append(pre_file_paths[j][char_index])
            if dif == 1:
                if len(dif_chars) == 2:
                    if "1" in dif_chars and "2" in dif_chars:
                        passed.append(i)
                        passed.append(j)
                        splited_1_a = pre_file_paths[i].split("_")
                        splited_1_b = splited_1_a[-1].split(".")
                        temp_path_1 = "{}/{}".format(folder, pre_file_paths[i])
                        temp_path_2 = "{}/{}".format(folder, pre_file_paths[j])
                        if splited_1_b[0] == "1":
                            file_paths_p[r] = [temp_path_1, temp_path_2]
                        else:
                            file_paths_p[r] = [temp_path_2, temp_path_1]
                        r += 1
    if len(passed) != passed_counter:
        print("Error. Not all files are paired-end reads or not all files could be mapped as paired-end reads due to their names being different in more than one character, which is expected to be number 1 and number 2.")
        exit()
    return file_paths_p


def file_paths_creation(input_folder, paired_end):
    file_paths = []
    file_paths_p = {}
    if input_folder[-1] == "/":
        input_folder = input_folder[:-1]
    # If the files are paired-end reads then except of the list file_paths which has the names of the files there is also a dictionary
    # which has as values lists with two paired end reads each.
    pre_file_paths = os.listdir(input_folder)
    for i in pre_file_paths:
        local_file_path = "{}/{}".format(input_folder, i)
        file_paths.append(local_file_path)
    if paired_end:
        file_paths_p = paired_end_file_paths(input_folder, False)
    return file_paths, file_paths_p
                

def locate_nc_files(paired_end, output_path_trimmed):
    tr_ex_file_paths_p = {}
    tr_ex_file_paths = []
    if paired_end:
        tr_ex_file_paths_p = paired_end_file_paths(output_path_trimmed, True)
    else:
        pre_tr_ex_file_paths = os.listdir(output_path_trimmed)
        for i in pre_tr_ex_file_paths:
            if "fastq" in i and i[-2:] != "gz":
                local_file_path = "{}/{}".format(output_path_trimmed, i)
                tr_ex_file_paths.append(local_file_path)
    return tr_ex_file_paths_p, tr_ex_file_paths