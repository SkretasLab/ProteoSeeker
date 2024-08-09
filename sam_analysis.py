import os
import sys


def analyse_file(sam_input_path):
    # All reads are present in the SAM file. Aligned (Mapped) and unaligned (unmapped) reads are present in the SAM file.
    # For paired-end reads, a read is considered unaligned if both reads of the pair are unaligned. Even if one is aligned
    # or ib both are aligned and discordant or both are aligned, the pair is considered aligned.
    # The unaligned reads of the SAM file can be found if the contig name is "*" and the 4th field is 0.
    # If both reads of a pair are aligned to different contigs, the contig of the read of the pair with the higher Phred score
    # (5th field) is retained.
    # Create a dictionary to correspond read IDs (headers) with contig IDs (headers).
    print("\nAnalyzing the SAM file...")
    reads_to_contigs_dict = {}
    contigs_to_reads_dict = {}
    with open(sam_input_path) as mapped_lines:
        for line in mapped_lines:
            line = line.rstrip()
            if line[0] == "@":
                continue
            else:
                line_splited = line.split("\t")
                read_id = line_splited[0]
                contig_id = line_splited[2]
                phred_score = line_splited[3]
                if (contig_id != "*") and (phred_score != 0):
                    if contig_id not in contigs_to_reads_dict.keys():
                        contigs_to_reads_dict[contig_id] = [0]
                    if read_id not in reads_to_contigs_dict.keys():
                        reads_to_contigs_dict[read_id] = [contig_id, phred_score]
                        contigs_to_reads_dict[contig_id].append(read_id)
                    else:
                        if contig_id != reads_to_contigs_dict[read_id]:
                            if phred_score > reads_to_contigs_dict[read_id][1]:
                                reads_to_contigs_dict[read_id][0] = contig_id
    read_sum = 0
    for key_contig in contigs_to_reads_dict.keys():
        read_temp_sum = len(contigs_to_reads_dict[key_contig])
        read_sum += read_temp_sum
    mapped_contig_num = len(list(contigs_to_reads_dict.keys()))
    print("\nContigs mapped to reads: {}".format(mapped_contig_num))
    print("Reads mapped to contigs: {}".format(read_sum))
    return reads_to_contigs_dict, contigs_to_reads_dict


def readmapal(sam_input_path):
    reads_to_contigs_dict = {}
    contigs_to_reads_dict = {}
    if not os.path.exists(sam_input_path):
        print("\nThe path to the input SAM file does not exists.")
        exit()
    if os.stat(sam_input_path).st_size == 0:
        print("\nThe input SAM file is empty.")
        exit()
    reads_to_contigs_dict, contigs_to_reads_dict = analyse_file(sam_input_path)
    return reads_to_contigs_dict, contigs_to_reads_dict


if __name__ == "__main__":
    arg_sam_input_path = ""
    if len(sys.argv) > 1:
        arg_input_command = "{}".format(sys.argv[0])
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input":
                sam_input_path = sys.argv[i+1]
    readmapal(sam_input_path)
