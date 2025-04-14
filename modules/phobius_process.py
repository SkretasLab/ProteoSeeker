import re
import os
import shutil
# ProteoSeeker modules
import command_process


def phobius(phobius_env, proteins_combined_file, output_path_topology, phobius_input_file_name, phobius_output_file_name, topology_info_path, phobius_path, input_log_file, output_log_file, phobius_bash_script, phobius_version_path, phobius_stde_path, conda_sh_path):
    print("\nRunning Phobius...")
    if os.path.exists(output_path_topology):
        shutil.rmtree(output_path_topology)
    os.mkdir(output_path_topology)
    # If there is a path to Phobius and a file with combined proteins, then use Phobius for transmembrane topology
    # predictions.
    dict_top = {}
    new_file = open(topology_info_path, "w")
    if os.path.exists(proteins_combined_file):
        shutil.copy(proteins_combined_file, phobius_input_file_name)
        # Command
        if os.path.exists(phobius_path):
            # In this case, at first export the path to the phobius folder.
            phrase_1 = "export PATH=\"$PATH:{}\"".format(phobius_path)
            phrase_2 = "phobius.pl \"{}\" > \"{}\" 2> \"{}\"".format(phobius_input_file_name, phobius_output_file_name, phobius_stde_path)
            phrase_3 = "phobius.pl -h > \"{}\"".format(phobius_version_path)
        else:
            phrase_1 = "perl phobius.pl \"{}\" > \"{}\" 2> \"{}\"".format(phobius_input_file_name, phobius_output_file_name, phobius_stde_path)
            phrase_2 = "perl phobius.pl -h > \"{}\"".format(phobius_version_path)
            phrase_3 = ""
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(phobius_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if phobius_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(phobius_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        if phrase_3:
            new_file_bash.write("{}\n".format(phrase_3))
        if phobius_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(phobius_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(phobius_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Gathering information from the results of Phobius
        pattern_1 = re.compile(r'^ID\s+(.*)')
        pattern_2 = re.compile(r'^FT\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)')
        pattern_3 = re.compile(r'^FT\s+(\S+)\s+(\S+)\s+(\S+)')
        with open(phobius_output_file_name) as phobius_output_lines:
            for line in phobius_output_lines:
                line = line.rstrip("\n")
                result_dom = pattern_1.search(line)
                if result_dom:
                    phobius_protein_id = result_dom.group(1)
                else:
                    region_name = None
                    result_2 = pattern_2.search(line)
                    result_3 = pattern_3.search(line)
                    if result_2:
                        region_name_1 = result_2.group(1)
                        left_bounder = result_2.group(2)
                        right_bounder = result_2.group(3)
                        region_name_2 = result_2.group(4)
                        if region_name_1:
                            region_name_1 = region_name_1.strip()
                        if region_name_2:
                            region_name_2 = region_name_2.strip()
                            if region_name_2[-1] == ".":
                                region_name_2 = region_name_2[:-1]
                            region_name = "{}||{}".format(region_name_1, region_name_2)
                        else:
                            region_name = region_name_1
                    elif result_3:
                            region_name = result_3.group(1)
                            left_bounder = result_3.group(2)
                            right_bounder = result_3.group(3)
                    if region_name is not None:
                        # Writing the information in the file.
                        new_file.write("{}\t{}\t{}\t{}\n".format(phobius_protein_id, left_bounder, right_bounder, region_name))
                        temp_phrase = "{}\t{}\t{}".format(left_bounder, right_bounder, region_name)
                        # Adding the information in the dictionary.
                        if phobius_protein_id not in dict_top.keys():
                            dict_top[phobius_protein_id] = []
                        dict_top[phobius_protein_id].append(temp_phrase)
    if not dict_top:
        new_file.write("No sequences with transmembrane or signal peptides found.\n")
    new_file.close()
    return dict_top