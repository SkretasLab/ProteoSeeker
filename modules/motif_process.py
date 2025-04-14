import re
import os
import shutil
# ProteoSeeker modules
import supportive_functions


def input_motif_search(motifs_path, output_path_motifs, input_motifs_results, proteins_combined_file, input_log_file):
    print("\nSearching for the input motifs...")
    if os.path.exists(output_path_motifs):
        shutil.rmtree(output_path_motifs)
    os.mkdir(output_path_motifs)
    protein_acc = None
    dict_input_motifs = {}
    if motifs_path is not None:
        # If the file with the CD-HIT filtered proteins has a non-zero size.
        if os.path.getsize(proteins_combined_file) != 0:
            motif_lines = supportive_functions.read_file(motifs_path)
            motif_sqs = []
            for line in motif_lines:
                # If the line is not empty and contains a motif then add it to the list.
                if line:
                    motif_sqs.append(line)
            print("\nInput motifs to be searched:")
            print(motif_sqs)
            if motif_sqs:
                motif_patterns = []
                for mof in motif_sqs:
                    mof_pat = "({})".format(mof)
                    pattern = re.compile(mof_pat)
                    motif_patterns.append(pattern)
                first_line = True
                new_file = open(input_motifs_results, "w")
                with open(proteins_combined_file) as lines_cd_hit:
                    for line in lines_cd_hit:
                        line = line.rstrip("\n")
                        if line[0] == ">":
                            if not first_line:
                                for pi in range(0, len(motif_patterns)):
                                    pattern = motif_patterns[pi]
                                    mof = motif_sqs[pi]
                                    result = re.finditer(pattern, protein_seq)
                                    if result:
                                        for item in result:
                                            # Numbering starts from 0 and the ending of the pattern found is one position
                                            # higher. So, for a numbering starting from 1, the start of the pattern is
                                            # one higher than the one given while the end remains the same.
                                            result_start = item.start() + 1
                                            result_end = item.end()
                                            temp_mof_list = [mof, result_start, result_end]
                                            if protein_acc not in dict_input_motifs.keys():
                                                dict_input_motifs[protein_acc] = [temp_mof_list]
                                            else:
                                                dict_input_motifs[protein_acc].append(temp_mof_list)
                                            new_file.write("{}\t{}\t{}\t{}\n".format(protein_acc, mof, result_start, result_end))
                            protein_acc = line[1:]
                            protein_seq = ""
                            first_line = False
                        else:
                            protein_seq += line
                # This is the same motif analysi as in the loop above but it is also set here because the last protein sequence must also be examined.
                # The whole sequence will be colleceted at the last iteration so that sequence should also be examined.
                for pi in range(0, len(motif_patterns)):
                    pattern = motif_patterns[pi]
                    mof = motif_sqs[pi]
                    result = re.finditer(pattern, protein_seq)
                    if result:
                        for item in result:
                            result_start = item.start() + 1
                            result_end = item.end()
                            temp_mof_list = [mof, result_start, result_end]
                            if protein_acc not in dict_input_motifs.keys():
                                dict_input_motifs[protein_acc] = [temp_mof_list]
                            else:
                                dict_input_motifs[protein_acc].append(temp_mof_list)
                            new_file.write("{}\t{}\t{}\t{}\n".format(protein_acc, mof, result_start, result_end))
                new_file.close()
            else:
                print("\nFile for motifs is empty. No input motifs are searched in the proteins.")
        else:
            print("\nNo proteins to be searched for motifs.")
            input_log_file.write("No proteins to be searched for motifs.\n")
    else:
        print("\nNo input motifs will be searched.\n")
        input_log_file.write("No input motifs were searched.\n")
    return dict_input_motifs