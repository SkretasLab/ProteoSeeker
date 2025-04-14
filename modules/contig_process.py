import os
import shutil
# ProteoSeeker modules
import command_process


def contig_formation(contigs, file_paths, output_final_contigs, output_megahit_contigs, output_final_contigs_formated, input_log_file, output_log_file, cat_path):
    print("\nLocating the file with the contigs...")
    # Copy contigs to contigs folder
    # If the input files contain files with contigs then their combination is the final contigs file.
    # If the input files contain only FASTQ files then the file with the contigs is produced by Megahit and it is found in the corresponding results folder (of Megahit).
    if contigs:
        if cat_path:
            phrase_1 = "\"{}\"".format(cat_path)
            phrase_2 = "\"{}\" --version".format(cat_path)
        else:
            phrase_1 = "cat"
            phrase_2 = "cat --version"
        # Making sure only Fasta files are copied to the folder with the contigs and not FastQ files.
        # We consider, based on the information of usage given for the program, that Fasta files include only contigs.
        for i in file_paths:
            if i[-5:] != "fastq":
                phrase_1 = "{} \"{}\"".format(phrase_1, i)
        phrase_1 = "\n{} > \"{}\"".format(phrase_1, output_final_contigs)
        capture_status = True
        shell_status = True
        pr_status = False
        title_1 = "Concatenate files - cat:"
        title_2 = "Version of cat:"
        command_process.command_run(phrase_1, phrase_2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    else:
        if os.path.exists(output_megahit_contigs):
            shutil.copy(output_megahit_contigs, output_final_contigs)
        else:
            print("\nNo files with contigs were generated from the assembly. Exiting.")
            exit()
    # Writing the final contigs in FASTA format.
    split_size = 70
    new_file_final_contigs_formated = open(output_final_contigs_formated, "w")
    with open(output_final_contigs, "r") as final_contigs_lines:
        for c_line in final_contigs_lines:
            c_line = c_line.rstrip("\n")
            if c_line:
                if c_line[0] == ">":
                    new_file_final_contigs_formated.write("{}\n".format(c_line))
                else:
                    c_fasta_splited = [c_line[ci:ci + split_size] for ci in range(0, len(c_line), split_size)]
                    for c_fasta_line in c_fasta_splited:
                        new_file_final_contigs_formated.write("{}\n".format(c_fasta_line))
    new_file_final_contigs_formated.close()