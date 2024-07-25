import os
import stat
import shutil


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
    par_demo_path = "par_demo.txt"
    kraken_demo_path = "run_kraken_demo.txt"
    comebin_demo_path = "run_comebin_demo.txt"
    metabinner_demo_path = "run_metabinner_demo.txt"
    all_demo_path = "run_all_demo.txt"

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
    dbs_path = "/mnt/sda2/Databases"
    kraken_8_db_path = "{}/kraken_standard_8_db".format(dbs_path)
    kraken_16_db_path = "{}/kraken_standard_16_db".format(dbs_path)
    kraken_72_db_path = "{}/kraken_standard_72_db".format(dbs_path)
    nr_db_path = "{}/nr_db/nr".format(dbs_path)

    base_output_path = "/home/compteam/ProteoSeeker/results"

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
                elif line == "protein_db_path=\"\"":
                    if bs_key in [3, 4]:
                        line = "protein_db_path=\"{}\"".format(nr_db_path)
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
                    line =  "output_path=\"{}\"".format(output_dir_path)
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
