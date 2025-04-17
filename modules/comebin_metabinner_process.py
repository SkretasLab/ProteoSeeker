import os
import shutil
import command_process
import supportive_functions


def binning(binning_tool, metabinner_env, comebin_env, contigs, protein_input, output_path_bin, metabinner_bin_path, comebin_bin_path, fullpath_contigs_formated, bin_bam_folder, enz_dir, output_path_conitgs, full_coverage_folder, fullpath_tr_fastq_files, full_bin_results_folder, full_coverage_profile_path, thread_num, input_log_file, output_log_file, metabinner_bash_script, comebin_bash_script, bin_ram_ammount, bin_num_contig_len, bin_num_kmer, bin_gen_coverage_stdoe_path, bin_gen_kmer_stdoe_path, bin_filter_tooshort_stdoe_path, binning_stdoe_path, comebin_batch_size, conda_sh_path, comebin_par, metabinner_par):
    print("\nPerforming binning with MetaBinner or COMEBin...")
    # Delete other files from the contigs directory, related to previous binning processes.
    if os.path.exists(output_path_conitgs):
        fn_group = os.listdir(output_path_conitgs)
        for fn in fn_group:
            if fn not in ["contigs.fa", "contigs_formated.fna"]:
                fn_path = "{}/{}".format(output_path_conitgs, fn)
                os.remove(fn_path)

    if contigs or protein_input:
        return 0
    if os.path.exists(output_path_bin):
        shutil.rmtree(output_path_bin)
    os.mkdir(output_path_bin)
    if os.path.exists(full_bin_results_folder):
        shutil.rmtree(full_bin_results_folder)

    if binning_tool == 1:
        # Tool parameters. Add whitespaces to the user-defined parameters, if needed.
        if metabinner_par is None:
            metabinner_par = " -t {} ".format(thread_num)
        else:
            if metabinner_par[0] != " ":
                metabinner_par = " {}".format(metabinner_par)
            if metabinner_par[-1] != " ":
                metabinner_par = "{} ".format(metabinner_par)
                
        # Change directory to the scripts of the metabinner environment.
        gen_coverage_file_path = "{}/scripts/gen_coverage_file.sh".format(metabinner_bin_path)
        gen_kmer_path = "{}/scripts/gen_kmer.py".format(metabinner_bin_path)
        filter_tooshort_path = "{}/scripts/Filter_tooshort.py".format(metabinner_bin_path)
        relpath_filtered_contigs_path = "{}/contigs_formated_kmer_4_f{}.csv".format(output_path_conitgs, bin_num_contig_len)
        full_filtered_contigs_path = os.path.abspath(relpath_filtered_contigs_path)
        run_metabinner_path = "{}/run_metabinner.sh".format(metabinner_bin_path)
        contnigs_formated_splited = fullpath_contigs_formated.split(".")
        contnigs_formated_splited = contnigs_formated_splited[:-1]
        contnigs_filtered_path = ".".join(contnigs_formated_splited)
        contnigs_filtered_path = "{}_{}.fa".format(contnigs_filtered_path, bin_num_contig_len)

        # Find the first line of the script for filtering contigs.
        # This indicates whether metabinner was installed from source from github or through anaconda.
        # If place_python is True then metabinner was installed from github, python is needed in front of the command to use the python of the environment but the 
        # full_coverage_profile_path path should be based on the selected contig length. If the place_python is False then metabinner was installed through anaconda,
        # no python is needed in front the command but the full_coverage_profile_path path should not change.
        place_python = False
        fts_lines = supportive_functions.read_file(filter_tooshort_path)
        first_line = fts_lines[0]
        if first_line == "#!/usr/bin/python":
            place_python = True
        else:
            full_coverage_profile_path = "{}/coverage_profile_f1k.tsv".format(full_coverage_folder)

        # Get the coverage profiles.
        phrase_1 = "\"{}\" -a \"{}\" -o \"{}\" -b \"{}\" -t {} -m {} -l {} {} &> \"{}\"".format(gen_coverage_file_path, fullpath_contigs_formated, full_coverage_folder, bin_bam_folder, thread_num, bin_ram_ammount, bin_num_contig_len, fullpath_tr_fastq_files, bin_gen_coverage_stdoe_path)
        # Create kmers.
        phrase_2 = "\"{}\" \"{}\" {} {} &> \"{}\"".format(gen_kmer_path, fullpath_contigs_formated, bin_num_contig_len, bin_num_kmer, bin_gen_kmer_stdoe_path)
        # Filter contigs based on their lengths.
        phrase_3 = "\"{}\" \"{}\" {} &> \"{}\"".format(filter_tooshort_path, fullpath_contigs_formated, bin_num_contig_len, bin_filter_tooshort_stdoe_path)
        if place_python:
            phrase_3 = "python {}".format(phrase_3)
        # Bin the contigs.
        phrase_4 = "\"{}\" -a \"{}\" -o \"{}\" -d \"{}\" -k \"{}\" -p \"{}\"{}&> \"{}\"".format(run_metabinner_path, contnigs_filtered_path, full_bin_results_folder, full_coverage_profile_path, full_filtered_contigs_path, metabinner_bin_path, metabinner_par, binning_stdoe_path)
        
        # Create the Bash script.
        new_file_bash = open(metabinner_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if metabinner_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(metabinner_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        new_file_bash.write("{}\n".format(phrase_3))
        new_file_bash.write("{}\n".format(phrase_4))
        if metabinner_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()

        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(metabinner_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(metabinner_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    else:
        # Tool parameters. Add whitespaces to the user-defined parameters, if needed.
        if comebin_par is None:
            comebin_par = " -t {} -b {} ".format(thread_num, comebin_batch_size)
        else:
            if comebin_par[0] != " ":
                comebin_par = " {}".format(comebin_par)
            if comebin_par[-1] != " ":
                comebin_par = "{} ".format(comebin_par)

        gen_coverage_file_path = "{}/scripts/gen_cov_file.sh".format(comebin_bin_path)
        filter_tooshort_path = "{}/scripts/Filter_tooshort.py".format(comebin_bin_path)
        run_comebin_path = "{}/run_comebin.sh".format(comebin_bin_path)
        contnigs_formated_splited = fullpath_contigs_formated.split(".")
        contnigs_formated_splited = contnigs_formated_splited[:-1]
        contnigs_filtered_path = ".".join(contnigs_formated_splited)
        contnigs_filtered_path = "{}_{}.fa".format(contnigs_filtered_path, bin_num_contig_len)

        # Change path to the comebin folder.
        phrase_1 = "cd \"{}\"".format(comebin_bin_path)
        # Get the coverage profiles.
        phrase_2 = "\"{}\" -a \"{}\" -o \"{}\" -b \"{}\" -t {} -m {} -l {} {} &> \"{}\"".format(gen_coverage_file_path, fullpath_contigs_formated, full_coverage_folder, bin_bam_folder, thread_num, bin_ram_ammount, bin_num_contig_len, fullpath_tr_fastq_files, bin_gen_coverage_stdoe_path)
        # Filter contigs based on their lengths.
        phrase_3 = "\"{}\" \"{}\" {} &> \"{}\"".format(filter_tooshort_path, fullpath_contigs_formated, bin_num_contig_len, bin_filter_tooshort_stdoe_path)
        # Bin the contigs.
        phrase_4 = "\"{}\" -a \"{}\" -o \"{}\" -p \"{}\"{}&> \"{}\"".format(run_comebin_path, contnigs_filtered_path, full_bin_results_folder, bin_bam_folder, comebin_par, binning_stdoe_path)
        # Change path to the proteoseeker folder
        phrase_5 = "cd \"{}\"".format(enz_dir)

        # Create the Bash script.
        new_file_bash = open(comebin_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if comebin_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(comebin_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        new_file_bash.write("{}\n".format(phrase_3))
        new_file_bash.write("{}\n".format(phrase_4))
        new_file_bash.write("{}\n".format(phrase_5))
        if comebin_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()

        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(comebin_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(comebin_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)