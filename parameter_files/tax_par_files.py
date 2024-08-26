import os
import stat
import shutil
import pathlib
import subprocess


def read_file(file_path):
    file_handle = open(file_path, "r")
    pre_lines = file_handle.readlines()
    lines = []
    for cur_line in pre_lines:
        cur_line = cur_line.rstrip("\n")
        lines.append(cur_line)
    file_handle.close()
    return lines


def cparf():
    prefix = "par"
    sra_codes_dict = {
        1: "SRR12829168",
        2: "SRR12829161",
        3: "SRR12829164",
        4: "SRR12829167",
        5: "SRR12829157",
        6: "SRR12829156",
        7: "SRR12829155",
        8: "SRR12829159",
        9: "SRR12829160",
        10: "SRR12829152",
        11: "SRR12829158",
        12: "SRR12829165",
        13: "SRR12829154",
        14: "SRR12829166",
        15: "SRR12829153",
        16: "SRR12829163",
        17: "SRR12829162",
        18: "SRR12829169",
        19: "SRR12829170"
    }
    basename_dict = {
        0: "phylo_k8",
        1: "phylo_k16",
        2: "phylo_k72",
        3: "phylo_m_nr",
        4: "phylo_c_nr"
    }

    # Getting the path to the parameter, the main proteseeker and the installation directory.
    parameter_path = pathlib.Path(__file__).parent.resolve()
    ps_path = parameter_path.parent.resolve()
    installation_path = "{}/installation".format(ps_path)

    # Getting the path for conda.
    conda_find_script = "{}/find_conda.sh".format(installation_path)
    command_phrase = "source {} && echo $CONDA_INST_DIR && echo $CONDA_SH_PATH".format(conda_find_script)
    result = subprocess.run(command_phrase, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    bash_output = result.stdout.strip().split('\n')
    conda_inst_dir = bash_output[0]
    conda_sh_path = bash_output[1]

    # Getting the paths for the COMEBin and MetaBinner bins.
    mb_cb_find_script = "{}/find_mb_cb.sh".format(installation_path)
    command_phrase = "source {} && echo $CTBPATH && echo $CBINPATH && echo $MTBPATH && echo $MBINPATH".format(mb_cb_find_script)
    result = subprocess.run(command_phrase, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    bash_output = result.stdout.strip().split('\n')
    ctbpath = bash_output[0]
    cbinpath = bash_output[1]
    mtbpath = bash_output[2]
    mbinpath = bash_output[3]

    # Getting the paths for the databases from the "setup.sh" script.
    setup_script = "{}/setup.sh".format(parameter_path)
    command_phrase = "source {} && echo $KRAKEN_8_DB_PATH && echo $KRAKEN_16_DB_PATH && echo $KRAKEN_72_DB_PATH && echo $PROTEIN_DB_PATH".format(setup_script)
    result = subprocess.run(command_phrase, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    bash_output = result.stdout.strip().split('\n')
    kraken_8_db_path = bash_output[0]
    kraken_16_db_path = bash_output[1]
    kraken_72_db_path = bash_output[2]
    protein_db_path = bash_output[3]

    # Setting the other paths.
    adapters_path = "{}/adapters.fa".format(ps_path)
    pfam_path = "{}/pfam_database/Pfam-A.hmm".format(ps_path)
    swissprot_path = "{}/swissprot_database/swissprot".format(ps_path)
    motifs_path = "{}/motifs.txt".format(ps_path)
    base_output_path = "{}/results".format(ps_path)
    fraggenescanrs_path = "{}/ps_tools/fgsrs/FragGeneScanRs".format(ps_path)
    phobius_path = "{}/ps_tools/phobius_files/phobius".format(ps_path)

    # The paths for the parameter files.
    par_demo_path = "{}/par_demo.txt".format(parameter_path)
    kraken_demo_path = "{}/run_kraken_demo.txt".format(parameter_path)
    comebin_demo_path = "{}/run_comebin_demo.txt".format(parameter_path)
    metabinner_demo_path = "{}/run_metabinner_demo.txt".format(parameter_path)
    all_demo_path = "{}/run_all_demo.txt".format(parameter_path)

    # Parameter lines for the taxonomy analysis
    par_lines = read_file(par_demo_path)
    # Parsing samples from 1 to 19.
    for sample_id in range(1, 20):
        print("Creating par files for index {}".format(sample_id))
        # Creating the directory for the sample, overwriting it if it already exists.
        dirpath = "{}".format(sample_id)
        if os.path.exists(dirpath):
            shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        # Parsing the parameter files for a sample.
        # Creating the parameter's file name
        for bs_key in basename_dict.keys():
            par_file_name = "{}_s{}_{}".format(prefix, sample_id, basename_dict[bs_key])
            par_file_path = "{}/{}.txt".format(dirpath, par_file_name)
            # 1. The first commands run for each method are about their smallest databases.
            # 2. After these commands are run, the output directory of the following command
            # is created.
            # 3. The initial and common directories from the first command are copied to the second one.
            # These directories include everything output up to the application of the binning method.
            # 4. The parameter file of the next command includes the name of the output directory and
            # is set to continue after the assembly for kraken and after the binning for COMEBin and
            # MetaBinner.
            # Creating the parameter files.
            new_file = open(par_file_path, "w")
            for line in par_lines:
                if line == "sra_code=\"\"":
                    sra_code = sra_codes_dict[sample_id]
                    line = "sra_code=\"{}\"".format(sra_code)
                elif line == "adapters_path=\"\"":
                    line = "adapters_path=\"{}\"".format(adapters_path)
                elif line == "profiles_broad_path=\"\"":
                    line = "profiles_broad_path=\"{}\"".format(pfam_path)
                elif line == "swissprot_path=\"\"":
                    line = "swissprot_path=\"{}\"".format(swissprot_path)
                elif line == "motifs_path=\"\"":
                    line = "motifs_path=\"{}\"".format(motifs_path)
                elif line == "conda_bin=\"\"":
                    line = "conda_bin=\"{}\"".format(conda_inst_dir)
                elif line == "conda_sh=\"\"":
                    line = "conda_sh=\"{}\"".format(conda_sh_path)
                elif line == "metabinner_bin_path=\"\"":
                    line = "metabinner_bin_path=\"{}\"".format(mbinpath)
                elif line == "comebin_bin_path=\"\"":
                    line = "comebin_bin_path=\"{}\"".format(cbinpath)
                elif line == "fraggenescanrs_path=\"\"":
                    line = "fraggenescanrs_path=\"{}\"".format(fraggenescanrs_path)
                elif line == "phobius_path=\"\"":
                    line = "phobius_path=\"{}\"".format(phobius_path)
                elif line == "protein_db_path=\"\"":
                    if bs_key in [3, 4]:
                        line = "protein_db_path=\"{}\"".format(protein_db_path)
                elif line == "kraken_db_path=\"\"":
                    if bs_key == 0:
                        line = "kraken_db_path=\"{}\"".format(kraken_8_db_path)
                    elif bs_key == 1:
                        line = "kraken_db_path=\"{}\"".format(kraken_16_db_path)
                    elif bs_key == 2:
                        line = "kraken_db_path=\"{}\"".format(kraken_72_db_path)
                elif line == "output_path=\"\"":
                    output_dir_path = "{}/sample_{}".format(base_output_path, sample_id)
                    if bs_key == 0:
                        output_dir_path = "{}/sample_{}_kraken_8".format(output_dir_path, sample_id)
                    elif bs_key == 1:
                        output_dir_path = "{}/sample_{}_kraken_16".format(output_dir_path, sample_id)
                    elif bs_key == 2:
                        output_dir_path = "{}/sample_{}_kraken_72".format(output_dir_path, sample_id)
                    elif bs_key == 3:
                        output_dir_path = "{}/sample_{}_metabinner_nr".format(output_dir_path, sample_id)
                    elif bs_key == 4:
                        output_dir_path = "{}/sample_{}_comebin_nr".format(output_dir_path, sample_id)
                    line =  "output_path=\"{}\"".format(output_dir_path)
                elif line == "family_code_phylo=\"\"":
                    if bs_key in [3, 4]:
                         line = "family_code_phylo=\"9050,9051,9052,9053,9054,9055,9056,9057,9058,9496\""
                elif line == "db_name_phylo=\"\"":
                    if bs_key in [3, 4]:
                         line = "db_name_phylo=\"nr_rna_pol\""
                elif line == "create_nr_db_status=\"\"":
                    if bs_key in [3, 4]:
                         line = "create_nr_db_status=\"True\""
                elif line == "kraken_mode=\"\"":
                    if bs_key in [3, 4]:
                         line = "kraken_mode=\"False\""
                elif line == "binning_tool=\"\"":
                    if bs_key == 4:
                         line = "binning_tool=\"2\""
                elif line == "after_gene_pred=\"\"":
                    if bs_key in [2, 4]:
                         line = "after_gene_pred=\"True\""
                new_file.write("{}\n".format(line))
            new_file.close()
        # Replacements
        i_str = str(sample_id)
        sample_str = "sample_{}".format(i_str)
        s_str = "s{}".format(i_str)
        slash_str = "/{}/".format(i_str)
        # Creating the Bash scripts.
        # Creating the kraken2 Bash script.
        bash_script_kraken_path = "run_{}_kraken.sh".format(sample_id)
        bash_script_kraken_file = open(bash_script_kraken_path, "w")
        kraken_demo_lines = read_file(kraken_demo_path)
        for line in kraken_demo_lines:
            line_r = line.replace("sample_X", sample_str)
            line_f = line_r.replace("sX", s_str)
            line_n = line_f.replace("/X/", slash_str)
            bash_script_kraken_file.write("{}\n".format(line_n))
        bash_script_kraken_file.close()
        file_stats = os.stat(bash_script_kraken_path)
        os.chmod(bash_script_kraken_path, file_stats.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        # Creating the COMEBin Bash script.
        bash_script_comebin_path = "run_{}_comebin.sh".format(sample_id)
        bash_script_comebin_file = open(bash_script_comebin_path, "w")
        comebin_demo_lines = read_file(comebin_demo_path)
        for line in comebin_demo_lines:
            line_r = line.replace("sample_X", sample_str)
            line_f = line_r.replace("sX", s_str)
            line_n = line_f.replace("/X/", slash_str)
            bash_script_comebin_file.write("{}\n".format(line_n))
        bash_script_comebin_file.close()
        file_stats = os.stat(bash_script_comebin_path)
        os.chmod(bash_script_comebin_path, file_stats.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        # Creating the MetaBinner Bash script.
        bash_script_metabinner_path = "run_{}_metabinner.sh".format(sample_id)
        bash_script_metabinner_file = open(bash_script_metabinner_path, "w")
        metabinner_demo_lines = read_file(metabinner_demo_path)
        for line in metabinner_demo_lines:
            line_r = line.replace("sample_X", sample_str)
            line_f = line_r.replace("sX", s_str)
            line_n = line_f.replace("/X/", slash_str)
            bash_script_metabinner_file.write("{}\n".format(line_n))
        bash_script_metabinner_file.close()
        file_stats = os.stat(bash_script_metabinner_path)
        os.chmod(bash_script_metabinner_path, file_stats.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        # Creating the total Bash script.
        bash_script_all_path = "run_{}_all.sh".format(sample_id)
        bash_script_all_file = open(bash_script_all_path, "w")
        all_demo_lines = read_file(all_demo_path)
        for line in all_demo_lines:
            line_n = line.replace("X", i_str)
            bash_script_all_file.write("{}\n".format(line_n))
        bash_script_all_file.close()
        file_stats = os.stat(bash_script_all_path)
        os.chmod(bash_script_all_path, file_stats.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


if __name__ == "__main__":
    cparf()
