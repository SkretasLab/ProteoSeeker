import os
import shutil
# ProteoSeeker modules
import command_process


def bbduk(bbduk_env, output_path_trimmed, file_paths, adapters_ap_path, paired_end, file_paths_p, bbduk_path, thread_num, bbduk_max_ram, output_path_summaries_errors, input_log_file, output_log_file, bbduk_bash_script, bbduk_version_path, bbduk_stdoe_path, conda_sh_path):
    print("\nPreprocessing reads with BBDuk...")
    if os.path.exists(output_path_trimmed):
        shutil.rmtree(output_path_trimmed)
    os.mkdir(output_path_trimmed)
    # Run bbduk.
    stats_file_counter = 1
    phrase_1 = ""
    if paired_end:
        for key in file_paths_p.keys():
            seq_fastq_1 = file_paths_p[key][0]
            seq_fastq_2 = file_paths_p[key][1]
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(bbduk_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if bbduk_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(bbduk_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            if bbduk_path:
                phrase_1 = "\"{}\" \"in={}\" \"in2={}\" \"out={}/trimmed_results_{}\" \"out2={}/trimmed_results_{}\" interleaved=f \"stats={}/trimmed_stats_file_{}\" \"ref={}\" ktrim=r k=21 mink=11 hammingdistance=2 qtrim=rl trimq=20 minavgquality=20 minlength=25 tpe=t tbo=t threads={} -Xmx{}g &> \"{}\"".format(bbduk_path, seq_fastq_1, seq_fastq_2, output_path_trimmed, seq_fastq_1.split("/")[-1], output_path_trimmed, seq_fastq_2.split("/")[-1], output_path_summaries_errors, stats_file_counter, adapters_ap_path, thread_num, bbduk_max_ram, bbduk_stdoe_path)
                phrase_2 = "\"{}\" > {}".format(bbduk_path, bbduk_version_path)
            else:
                phrase_1 = "bbduk.sh \"in={}\" \"in2={}\" \"out={}/trimmed_results_{}\" \"out2={}/trimmed_results_{}\" interleaved=f \"stats={}/trimmed_stats_file_{}\" \"ref={}\" ktrim=r k=21 mink=11 hammingdistance=2 qtrim=rl trimq=20 minavgquality=20 minlength=25 tpe=t tbo=t threads={} -Xmx{}g &> \"{}\"".format(seq_fastq_1, seq_fastq_2, output_path_trimmed, seq_fastq_1.split("/")[-1], output_path_trimmed, seq_fastq_2.split("/")[-1], output_path_summaries_errors, stats_file_counter, adapters_ap_path, thread_num, bbduk_max_ram, bbduk_stdoe_path)
                phrase_2 = "bbduk.sh > \"{}\"".format(bbduk_version_path)
            stats_file_counter += 1
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            if bbduk_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(bbduk_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(bbduk_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    else:
        for i in file_paths:
            # Create the Bash script.
            # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
            new_file_bash = open(bbduk_bash_script, "w")
            phrase = "#!/bin/bash"
            new_file_bash.write("{}\n".format(phrase))
            if bbduk_env:
                phrase_s = "source \"{}\"".format(conda_sh_path)
                phrase_a = "conda activate \"{}\"".format(bbduk_env)
                new_file_bash.write("{}\n".format(phrase_s))
                new_file_bash.write("{}\n".format(phrase_a))
            if bbduk_path:
                phrase_1 = "\"{}\" \"in={}\" \"out={}/trimmed_results_{}\" \"interleaved=f stats={}/trimmed_stats_file_{}\" \"ref={}\" ktrim=r k=21 mink=11 hammingdistance=2 qtrim=rl trimq=20 minavgquality=20 minlength=25 tpe=t tbo=t threads={} -Xmx{}g &> \"{}\"".format(bbduk_path, i, output_path_trimmed, i.split("/")[-1], output_path_summaries_errors, stats_file_counter, adapters_ap_path, thread_num, bbduk_max_ram, bbduk_stdoe_path)
                phrase_2 = "\"{}\" --version > {}".format(bbduk_path, bbduk_version_path)
            else:
                phrase_1 = "bbduk.sh \"in={}\" \"out={}/trimmed_results_{}\" \"interleaved=f stats={}/trimmed_stats_file_{}\" \"ref={}\" ktrim=r k=21 mink=11 hammingdistance=2 qtrim=rl trimq=20 minavgquality=20 minlength=25 tpe=t tbo=t threads={} -Xmx{}g &> \"{}\"".format(i, output_path_trimmed, i.split("/")[-1], output_path_summaries_errors, stats_file_counter, adapters_ap_path, thread_num, bbduk_max_ram, bbduk_stdoe_path)
                phrase_2 = "bbduk.sh --version > {}".format(bbduk_version_path)
            stats_file_counter += 1
            new_file_bash.write("{}\n".format(phrase_1))
            new_file_bash.write("{}\n".format(phrase_2))
            if bbduk_env:
                phrase = "conda deactivate"
                new_file_bash.write("{}\n".format(phrase))
            new_file_bash.close()
            # Making the Bash script executable.
            # Sending command to run.
            phrase_b1 = "chmod +x {}".format(bbduk_bash_script)
            phrase_b2 = "chmod --version"
            title_1 = "Making a Bash script executable:"
            title_2 = "Version of chmod:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
            # Sending command to run.
            phrase_b1 = "bash {}".format(bbduk_bash_script)
            phrase_b2 = "bash --version"
            title_1 = "Running bash:"
            title_2 = "Version of bash:"
            capture_status = True
            shell_status = True
            pr_status = False
            command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)