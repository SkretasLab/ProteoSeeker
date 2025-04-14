import os
import shutil
# ProteoSeeker modules
import command_process


def megahit(megahit_env, output_path_megahit, tr_ex_file_paths, paired_end, tr_ex_file_paths_p, megahit_path, k_list, thread_num, input_log_file, output_log_file, megahit_bash_script, megahit_version_path, megahit_stdoe_path, conda_sh_path):
    print("\nRunning Megahit...")
    # Megahit (with the option "-o" as used in this case) demands that the output folder does not already exists. In this case the output folder its not "output_path_megahit", its "output_path_megahit/Megahit_Contigs/".
    if os.path.exists(output_path_megahit):
        shutil.rmtree(output_path_megahit)
    os.mkdir(output_path_megahit)
    # --k-list 15,17,21,29,39,59,79,99,119,141 --merge-level=20,0.7
    # --k-list 15,17,21,29,39,59,79,99,119,141
    # --presets meta-large: --min-count 1 --k-list 21,29,39,49,...,129,141
    # --presets meta-sensitive: --k-min 27 --k-max 127 --k-step 10
    if paired_end:
        if megahit_path:
            phrase_1 = "\"{}\" -t {}".format(megahit_path, thread_num)
            phrase_2 = "\"{}\" --version > \"{}\"".format(megahit_path, megahit_version_path)
        else:
            phrase_1 = "megahit -t {}".format(thread_num)
            phrase_2 = "megahit --version > \"{}\"".format(megahit_version_path)
        if k_list is not None:
            phrase_1 = "{} \"--k-list={}\"".format(phrase_1, k_list)
        phrase_1 = "{} -1 \"".format(phrase_1)
        tr_ex_file_paths_keys = list(tr_ex_file_paths_p.keys())
        for key_i in range(0, len(tr_ex_file_paths_keys)):
            key = tr_ex_file_paths_keys[key_i]
            if key_i + 1 == len(tr_ex_file_paths_keys):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths_p[key][0])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths_p[key][0])
        phrase_1 = "{} -2 \"".format(phrase_1)
        for key_i in range(0, len(tr_ex_file_paths_p.keys())):
            key = tr_ex_file_paths_keys[key_i]
            if key_i + 1 == len(tr_ex_file_paths_keys):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths_p[key][1])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths_p[key][1])
        phrase_1 = "{} --presets meta-large -o {}/megahit_contigs &> {}".format(phrase_1, output_path_megahit, megahit_stdoe_path)
    else:
        if megahit_path:
            phrase_1 = "\"{}\" -t {}".format(megahit_path, thread_num)
            phrase_2 = "\"{}\" --version > \"{}\"".format(megahit_path, megahit_version_path)
        else:
            phrase_1 = "megahit -t {}".format(thread_num)
            phrase_2 = "megahit --version > \"{}\"".format(megahit_version_path)
        if k_list is not None:
            phrase_1 = "{} \"--k-list={}\"".format(phrase_1, k_list)
        phrase_1 = "{} -r \"".format(phrase_1)
        for i in range(0, len(tr_ex_file_paths)):
            if i + 1 == len(tr_ex_file_paths):
                phrase_1 = "{}{}\"".format(phrase_1, tr_ex_file_paths[i])
            else:
                phrase_1 = "{}{},".format(phrase_1, tr_ex_file_paths[i])
        phrase_1 = "{} --presets meta-large -o \"{}/megahit_contigs\" &> \"{}\"".format(phrase_1, output_path_megahit, megahit_stdoe_path)
    # Create the Bash script.
    # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
    new_file_bash = open(megahit_bash_script, "w")
    phrase = "#!/bin/bash"
    new_file_bash.write("{}\n".format(phrase))
    if megahit_env:
        phrase_s = "source \"{}\"".format(conda_sh_path)
        phrase_a = "conda activate \"{}\"".format(megahit_env)
        new_file_bash.write("{}\n".format(phrase_s))
        new_file_bash.write("{}\n".format(phrase_a))
    new_file_bash.write("{}\n".format(phrase_1))
    new_file_bash.write("{}\n".format(phrase_2))
    if megahit_env:
        phrase = "conda deactivate"
        new_file_bash.write("{}\n".format(phrase))
    new_file_bash.close()
    # Making the Bash script executable.
    # Sending command to run.
    phrase_b1 = "chmod +x {}".format(megahit_bash_script)
    phrase_b2 = "chmod --version"
    title_1 = "Making a Bash script executable:"
    title_2 = "Version of chmod:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Sending command to run.
    phrase_b1 = "bash {}".format(megahit_bash_script)
    phrase_b2 = "bash --version"
    title_1 = "Running bash:"
    title_2 = "Version of bash:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)