import os
import shutil
# ProteoSeeker modules
import command_process


def collect_sra(sra_env, sra_folder, prefetch_size, sra_run_folder, sra_file, fastq_folder, fastq_single_end, sra_code, input_log_file, output_log_file, prefetch_path, vdb_validate_path, fastq_dump_path, sra_bash_script, conda_sh_path):
    print("\nSRA code selected: {}".format(sra_code))
    if not os.path.exists(sra_folder):
        os.mkdir(sra_folder)
    # Create a folder for the specific run.
    if os.path.exists(fastq_folder):
        print("\nThe folder with the FASTQ files of the specified run already exists. Analysis is continued based on the FASTQ files present in that folder ('{}').".format(fastq_folder))
    else:
        os.mkdir(sra_run_folder)
        prefetch_version_path = "{}/prefetch_version.txt".format(sra_run_folder)
        vdb_vaildate_version_path = "{}/vdb_vaildate_version.txt".format(sra_run_folder)
        fastq_dump_version_path = "{}/fastq_dump_version.txt".format(sra_run_folder)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(sra_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if sra_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(sra_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        if prefetch_path:
            phrase_1 = "\"{}\" -X {}G -O \"{}\" \"{}\"".format(prefetch_path, prefetch_size, sra_run_folder, sra_code)
            phrase_2 = "\"{}\" --version > \"{}\"".format(prefetch_path, prefetch_version_path)
        else:
            phrase_1 = "prefetch -X {}G -O \"{}\" \"{}\"".format(prefetch_size, sra_run_folder, sra_code)
            phrase_2 = "prefetch --version > \"{}\"".format(prefetch_version_path)
        if vdb_validate_path:
            phrase_3 = "\"{}\" \"{}/\"".format(vdb_validate_path, sra_file)
            phrase_4 = "\"{}\" --version > \"{}\"".format(vdb_validate_path, vdb_vaildate_version_path)
        else:
            phrase_3 = "vdb-validate \"{}/\"".format(sra_file)
            phrase_4 = "vdb-validate --version > \"{}\"".format(vdb_vaildate_version_path)
        if fastq_dump_path:
            phrase_5 = "\"{}\" --skip-technical --gzip --split-3 -O \"{}\" \"{}\"".format(fastq_dump_path, fastq_folder, sra_file)
            phrase_6 = "\"{}\" --version > \"{}\"".format(fastq_dump_path, fastq_dump_version_path)
        else:
            phrase_5 = "fastq-dump --skip-technical --gzip --split-3 -O \"{}\" \"{}\"".format(fastq_folder, sra_file)
            phrase_6 = "fastq-dump --version > {}".format(fastq_dump_version_path)
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        new_file_bash.write("{}\n".format(phrase_3))
        new_file_bash.write("{}\n".format(phrase_4))
        new_file_bash.write("{}\n".format(phrase_5))
        new_file_bash.write("{}\n".format(phrase_6))
        if sra_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(sra_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(sra_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Determine whether the files are SINGLE or PAIRED-END reads.
    fastq_files = os.listdir(fastq_folder)
    paired_end = None
    paired_end_size = 0
    single_end_size = 0
    # If there are 3 files in the folder, then 2 of them are paired-end and one single-end.
    if len(fastq_files) == 3:
        if os.path.exists(fastq_single_end):
            shutil.rmtree(fastq_single_end)
        os.mkdir(fastq_single_end)
        for ffi in fastq_files:
            ffi_path = "{}/{}".format(fastq_folder, ffi)
            if (ffi[-11:] != "_1.fastq.gz") and (ffi[-11:] != "_2.fastq.gz"):
                new_ffi_path = "{}/{}".format(fastq_single_end, ffi)
                shutil.move(ffi_path, new_ffi_path)
                pes = os.path.getsize(new_ffi_path)
                single_end_size += pes
            else:
                pes = os.path.getsize(ffi_path)
                paired_end_size += pes
        paired_end = True
        if single_end_size >= (0.5 * paired_end_size):
            print("\nThe size of the single-end file (forward or reverse reads whose pair was not found) is equal or greater than 50% of the size of the paired-end file. Exiting.")
            exit()
    elif len(fastq_files) == 2:
        paired_end = True
    elif len(fastq_files) == 1:
        paired_end = False
    else:
        print("\nError. Could not determine whether the files are single-end or paired-end.")
        exit()
    if paired_end:
        print("\nThe files are paired-end.")
    else:
        print("\nThe file is single-end.")
    return paired_end, fastq_folder