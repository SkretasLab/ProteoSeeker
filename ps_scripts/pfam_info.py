import re
import sys


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
    print("pfam_info.py Version 1.0.0")
    print()
    print("Usage:")
    print("python pfam_info.py -i <pfam_hmm_file_path>")
    print()
    print("This tool can be used to process information from the Pfam-A.hmm file.")
    print("It outputs a series of files each containing different types of information")
    print("based on the Pfam-A.hmm file. The output files are created in the same")
    print("directory where this Python script is run.")
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
        "-i/--input": "Str -Req: pfam.dat- The path to the Pfam-A.hmm file."
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


def pfamdata(pfam_path):
    pattern_name = re.compile(r'^NAME\s+(.*)')
    pattern_acc = re.compile(r'^ACC\s+(.*)')
    pattern_len = re.compile(r'^LENG\s+(.*)')
    found_name = False
    found_acc = False
    name_acc_dict = {}
    counter_lines = 0
    with open(pfam_path) as pfam_lines:
        for line in pfam_lines:
            counter_lines += 1
            if(counter_lines % 100000 == 0):
                print("Line {} of 10332549".format(counter_lines))
            line = line.rstrip()
            result_name = pattern_name.search(line)
            result_acc = pattern_acc.search(line)
            result_len = pattern_len.search(line)
            if result_name:
                found_name = True
                temp_name = result_name.group(1)
            elif result_acc:
                if not found_name:
                    print("Found accession number without finding a new name prior to the accession number. Exiting.")
                    exit()
                found_name = False
                found_acc = True
                temp_acc = result_acc.group(1)
                if temp_acc in name_acc_dict.keys():
                    print("Duplicate accession number. Exiting.")
                    exit()
                name_acc_dict[temp_acc] = [temp_name]
            elif result_len:
                if not found_acc:
                    print("Found length without finding a new accession number prior to the length. Exiting.")
                    exit()
                found_acc = False
                temp_len = result_len.group(1)
                name_acc_dict[temp_acc].append(temp_len)
    print("Line {} of 10332549".format(counter_lines))
    print()

    new_file_name = "pfam_accs_names.tsv"
    new_file_1 = open(new_file_name, "w")
    for key in name_acc_dict.keys():
        domain_acc = key
        domain_name = name_acc_dict[key][0]
        new_file_1.write("{}\t{}\n".format(domain_acc, domain_name))
    new_file_1.close()

    new_file_name = "profiles_lengths.tsv"
    new_file_2 = open(new_file_name, "w")
    for key in name_acc_dict.keys():
        domain_acc = key
        domain_name = name_acc_dict[key][0]
        domain_length = name_acc_dict[key][1]
        new_file_2.write("{}\t{}\t{}\n".format(domain_name, domain_acc, domain_length))
    new_file_2.close()


if __name__ == "__main__":
    pfam_path = "Pfam-A.hmm"
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input":
                pfam_path = sys.argv[i+1]
            elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
                help_message()
                exit()
    pfamdata(pfam_path)

