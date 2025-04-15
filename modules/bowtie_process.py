import os
import shutil
# ProteoSeeker modules
import command_process
import sam_file_process


def bowtie(output_path_bowtie, paired_end, tr_ex_file_paths, tr_ex_file_paths_p, conda_sh_path, bowtie_build_path, bowtie_build_version_path, bowtie_build_stdoe_path, bowtie_env, bowtie_bash_script, bowtie_path, bowtie_version_path, bowtie_stdoe_path, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, mapped_reads_path, thread_num, output_final_contigs_formated, bowtie_contigs_basename, bowtie_status, input_log_file, output_log_file):
    print("\nRunning bowtie2...")
    reads_to_contigs_dict = {}
    contigs_to_reads_dict = {}
    if bowtie_status:
        if os.path.exists(output_path_bowtie):
            shutil.rmtree(output_path_bowtie)
        os.mkdir(output_path_bowtie)
    if os.path.exists(output_final_contigs_formated):
        if bowtie_status:
            # bowtie_single_unaligned_path: Unpaired reads (meaning single-end reads) that did not align.
            # bowtie_paired_unaligned_con_path: Paired reads that did not align concordantly.
            # where are the paired reads that did not align at all?
            if bowtie_build_path:
                phrase_1 = "\"{}\" --threads {} \"{}\" \"{}\" &> \"{}\"".format(bowtie_build_path, thread_num, output_final_contigs_formated, bowtie_contigs_basename, bowtie_build_stdoe_path)
                phrase_2 = "\"{}\" --version > \"{}\"".format(bowtie_build_path, bowtie_build_version_path)
            else:
                phrase_1 = "bowtie2-build --threads {} \"{}\" \"{}\" &> \"{}\"".format(thread_num, output_final_contigs_formated, bowtie_contigs_basename, bowtie_build_stdoe_path)
                phrase_2 = "bowtie2-build --version > \"{}\"".format(bowtie_build_version_path)
            if paired_end:
                paired_one_phrase = ""
                tr_ex_file_paths_keys = list(tr_ex_file_paths_p.keys())
                for key_file in range(0, len(tr_ex_file_paths_keys)):
                    key = tr_ex_file_paths_keys[key_file]
                    file_path_1 = tr_ex_file_paths_p[key][0]
                    if key_file + 1 == len(tr_ex_file_paths_keys):
                        paired_one_phrase = "{}{}".format(paired_one_phrase, file_path_1)
                    else:
                        paired_one_phrase = "{}{},".format(paired_one_phrase, file_path_1)
                paired_two_phrase = ""
                tr_ex_file_paths_keys = list(tr_ex_file_paths_p.keys())
                for key_file in range(0, len(tr_ex_file_paths_keys)):
                    key = tr_ex_file_paths_keys[key_file]
                    file_path_2 = tr_ex_file_paths_p[key][1]
                    if key_file + 1 == len(tr_ex_file_paths_keys):
                        paired_two_phrase = "{}{}".format(paired_two_phrase, file_path_2)
                    else:
                        paired_two_phrase = "{}{},".format(paired_two_phrase, file_path_2)
                if bowtie_path:
                    phrase_3 = "\"{}\" -x \"{}\" -1 \"{}\" -2 \"{}\" --very-sensitive --un \"{}\" --un-conc \"{}\" --met-file \"{}\" --omit-sec-seq -p {} -S \"{}\" &> \"{}\"".format(bowtie_path, bowtie_contigs_basename, paired_one_phrase, paired_two_phrase, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, thread_num, mapped_reads_path, bowtie_stdoe_path)
                    phrase_4 = "\"{}\" --version > \"{}\"".format(bowtie_path, bowtie_version_path)
                else:
                    phrase_3 = "bowtie2 -x \"{}\" -1 \"{}\" -2 \"{}\" --very-sensitive --un \"{}\" --un-conc \"{}\" --met-file \"{}\" --omit-sec-seq -p {} -S \"{}\" &> \"{}\"".format(bowtie_contigs_basename, paired_one_phrase, paired_two_phrase, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, thread_num, mapped_reads_path, bowtie_stdoe_path)
                    phrase_4 = "bowtie2 --version > \"{}\"".format(bowtie_version_path)
            else:
                single_phrase = ""
                for i in range(0, len(tr_ex_file_paths)):
                    file_path = tr_ex_file_paths[i]
                    if i + 1 == len(tr_ex_file_paths):
                        single_phrase = "{}{}".format(single_phrase, file_path)
                    else:
                        single_phrase = "{}{},".format(single_phrase, file_path)
                if bowtie_path:
                    phrase_3 = "\"{}\" -x \"{}\" -U \"{}\" --very-sensitive --un \"{}\" --un-conc \"{}\" --met-file \"{}\" --omit-sec-seq -S \"{}\" -p {} &> \"{}\"".format(bowtie_path, bowtie_contigs_basename, single_phrase, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, mapped_reads_path, thread_num, bowtie_stdoe_path)
                    phrase_4 = "\"{}\" --version > \"{}\"".format(bowtie_path, bowtie_version_path)
                else:
                    phrase_3 = "bowtie2 -x \"{}\" -U \"{}\" --very-sensitive --un \"{}\" --un-conc \"{}\" --met-file \"{}\" --omit-sec-seq -S \"{}\" -p {} &> \"{}\"".format(bowtie_contigs_basename, single_phrase, bowtie_single_unaligned_path, bowtie_paired_unaligned_con_path, bowtie_stats_path, mapped_reads_path, thread_num, bowtie_stdoe_path)
                    phrase_4 = "bowtie2 --version > \"{}\"".format(bowtie_version_path)
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(bowtie_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if bowtie_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(bowtie_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            new_file_bash.write("{}\n".format(phrase_3))
            new_file_bash.write("{}\n".format(phrase_4))
            if bowtie_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(bowtie_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(bowtie_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Find the aligned reads.
        reads_to_contigs_dict, contigs_to_reads_dict = sam_file_process.readmapal(mapped_reads_path)
    else:
        print("\nOne of the files needed for the Bowtie analysis is missing.")
    return reads_to_contigs_dict, contigs_to_reads_dict
