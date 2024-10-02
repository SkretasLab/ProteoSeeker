import os
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


def modify_par_files(target_dir, output_dict, ps_path, adapters_path, protein_db_path, pfam_path, swissprot_path, motifs_path, conda_inst_dir, conda_sh_path, mbinpath, cbinpath, fraggenescanrs_path, phobius_path, aldiv_path, kraken_8_db_path, kraken_16_db_path, kraken_72_db_path):
    target_dir_files = os.listdir(target_dir)
    for sfn in target_dir_files:
        sfn_full_path = "{}/{}".format(target_dir, sfn)
        sfn_lines = read_file(sfn_full_path)
        sfn_file = open(sfn_full_path, "w")
        output_type = None
        if "sra_dbs" in sfn:
            output_type = 0
        elif "k8_0c01" in sfn:
            output_type = 2
        elif "k16_0c01" in sfn:
            output_type = 4
        elif "k72_0c01" in sfn:
            output_type = 6
        elif "k8" in sfn:
            output_type = 1
        elif "k16" in sfn:
            output_type = 3
        elif "k72" in sfn:
            output_type = 5
        elif "cnr" in sfn:
            output_type = 7
        elif "mnr" in sfn:
            output_type = 8
        for line in sfn_lines:
            if "adapters_path=" in line:
                line = "adapters_path=\"{}\"".format(adapters_path)
            elif "protein_db_path=" in line:
                line = "protein_db_path=\"{}\"".format(protein_db_path)
            elif "profiles_broad_path=" in line:
                line = "profiles_broad_path=\"{}\"".format(pfam_path)
            elif "swissprot_path=" in line:
                line = "swissprot_path=\"{}\"".format(swissprot_path)
            elif "motifs_path=" in line:
                line = "motifs_path=\"{}\"".format(motifs_path)
            elif "conda_bin=" in line:
                line = "conda_bin=\"{}\"".format(conda_inst_dir)
            elif "conda_sh=" in line:
                line = "conda_sh=\"{}\"".format(conda_sh_path)
            elif "alpha_diversity_path=" in line:
                line = "alpha_diversity_path=\"{}\"".format(aldiv_path)
            elif "metabinner_bin_path=" in line:
                line = "metabinner_bin_path=\"{}\"".format(mbinpath)
            elif "comebin_bin_path=" in line:
                line = "comebin_bin_path=\"{}\"".format(cbinpath)
            elif "fraggenescanrs_path=" in line:
                line = "fraggenescanrs_path=\"{}\"".format(fraggenescanrs_path)
            elif "phobius_path=" in line:
                line = "phobius_path=\"{}\"".format(phobius_path)
            elif "kraken_db_path=" in line:
                if "k8" in sfn:
                    line = "kraken_db_path=\"{}\"".format(kraken_8_db_path)
                elif "k16" in sfn:
                    line = "kraken_db_path=\"{}\"".format(kraken_16_db_path)
                elif "k72" in sfn:
                    line = "kraken_db_path=\"{}\"".format(kraken_72_db_path)
            elif "output_path=" in line:
                output_part = output_dict[output_type]
                output_dir_path = "{}/results/{}".format(ps_path, output_part)
                line =  "output_path=\"{}\"".format(output_dir_path)
            sfn_file.write("{}\n".format(line))
        sfn_file.close()


def cparf():
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
    command_phrase = "source {} && echo $CBINPATH && echo $MBINPATH".format(mb_cb_find_script)
    result = subprocess.run(command_phrase, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    bash_output = result.stdout.strip().split('\n')
    cbinpath = bash_output[0]
    mbinpath = bash_output[1]

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
    fraggenescanrs_path = "{}/ps_tools/fgsrs/FragGeneScanRs".format(ps_path)
    phobius_path = "{}/ps_tools/phobius_files/phobius".format(ps_path)
    aldiv_path = "{}/ps_tools/KrakenTools/DiversityTools/alpha_diversity.py".format(ps_path)

    # Parameter lines for the seek/taxonomy analysis and the analysis run.
    # SRR3961740
    srr3961740_dir = "{}/cas_als/ca_run/SRR3961740".format(parameter_path)
    output_dict_srr3961740 = {
        0: "results_ca_srr3961740_sra_dbs",
        1: "results_ca_srr3961740_k8",
        2: "results_ca_srr3961740_k8_0c01",
        3: "results_ca_srr3961740_k16",
        4: "results_ca_srr3961740_k16_0c01",
        5: "results_ca_srr3961740_k72",
        6: "results_ca_srr3961740_k72_0c01",
        7: "results_ca_srr3961740_cnr",
        8: "results_ca_srr3961740_mnr"
    }
    modify_par_files(srr3961740_dir, output_dict_srr3961740, ps_path, adapters_path, protein_db_path, pfam_path, swissprot_path, motifs_path, conda_inst_dir, conda_sh_path, mbinpath, cbinpath, fraggenescanrs_path, phobius_path, aldiv_path, kraken_8_db_path, kraken_16_db_path, kraken_72_db_path)
    # Drr163688
    drr163688_dir = "{}/cas_als/ca_run/DRR163688".format(parameter_path)
    output_dict_drr163688 = {
        0: "results_ca_drr163688_sra_dbs",
        1: "results_ca_drr163688_k8",
        2: "results_ca_drr163688_k8_0c01",
        3: "results_ca_drr163688_k16",
        4: "results_ca_drr163688_k16_0c01",
        5: "results_ca_drr163688_k72",
        6: "results_ca_drr163688_k72_0c01",
        7: "results_ca_drr163688_cnr",
        8: "results_ca_drr163688_mnr"
    }
    modify_par_files(drr163688_dir, output_dict_drr163688, ps_path, adapters_path, protein_db_path, pfam_path, swissprot_path, motifs_path, conda_inst_dir, conda_sh_path, mbinpath, cbinpath, fraggenescanrs_path, phobius_path, aldiv_path, kraken_8_db_path, kraken_16_db_path, kraken_72_db_path)
    # SRR17771278
    srr17771278_dir = "{}/cas_als/al_run/SRR17771278".format(parameter_path)
    output_dict_srr17771278 = {
        0: "results_al_srr17771278_sra_dbs",
        1: "results_al_srr17771278_k8",
        2: "results_al_srr17771278_k8_0c01",
        3: "results_al_srr17771278_k16",
        4: "results_al_srr17771278_k16_0c01",
        5: "results_al_srr17771278_k72",
        6: "results_al_srr17771278_k72_0c01",
        7: "results_al_srr17771278_cnr",
        8: "results_al_srr17771278_mnr"
    }
    modify_par_files(srr17771278_dir, output_dict_srr17771278, ps_path, adapters_path, protein_db_path, pfam_path, swissprot_path, motifs_path, conda_inst_dir, conda_sh_path, mbinpath, cbinpath, fraggenescanrs_path, phobius_path, aldiv_path, kraken_8_db_path, kraken_16_db_path, kraken_72_db_path)


if __name__ == "__main__":
    cparf()
