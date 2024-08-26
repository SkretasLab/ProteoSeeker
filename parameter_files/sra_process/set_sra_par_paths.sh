#!/bin/bash
set -e

# Get the current directory
SRA_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARAMETER_PATH=$(dirname "${SRA_PATH}")
PS_PATH=$(dirname "${PARAMETER_PATH}")
GEN_RESULTS_PATH="${PS_PATH}/results"
SRA_RESULTS_PATH="${GEN_RESULTS_PATH}/samples_sra"
INSTALLATION_PATH="${PS_PATH}/installation"
FIND_COND_SCRIPT="${INSTALLATION_PATH}/find_conda.sh"

# Source the script that located the conda paths.
source "${FIND_COND_SCRIPT}"

# Sample 1
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_1_sra\"|" "${PARAMETER_PATH}/sra_process/par_s1_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s1_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s1_sra.txt"

# Sample 2
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_2_sra\"|" "${PARAMETER_PATH}/sra_process/par_s2_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s2_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s2_sra.txt"

# Sample 3
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_3_sra\"|" "${PARAMETER_PATH}/sra_process/par_s3_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s3_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s3_sra.txt"

# Sample 4
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_4_sra\"|" "${PARAMETER_PATH}/sra_process/par_s4_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s4_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s4_sra.txt"

# Sample 5
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_5_sra\"|" "${PARAMETER_PATH}/sra_process/par_s5_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s5_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s5_sra.txt"

# Sample 6
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_6_sra\"|" "${PARAMETER_PATH}/sra_process/par_s6_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s6_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s6_sra.txt"

# Sample 7
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_7_sra\"|" "${PARAMETER_PATH}/sra_process/par_s7_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s7_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s7_sra.txt"

# Sample 8
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_8_sra\"|" "${PARAMETER_PATH}/sra_process/par_s8_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s8_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s8_sra.txt"

# Sample 9
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_9_sra\"|" "${PARAMETER_PATH}/sra_process/par_s9_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s9_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s9_sra.txt"

# Sample 10
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_10_sra\"|" "${PARAMETER_PATH}/sra_process/par_s10_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s10_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s10_sra.txt"

# Sample 11
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_11_sra\"|" "${PARAMETER_PATH}/sra_process/par_s11_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s11_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s11_sra.txt"

# Sample 12
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_12_sra\"|" "${PARAMETER_PATH}/sra_process/par_s12_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s12_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s12_sra.txt"

# Sample 13
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_13_sra\"|" "${PARAMETER_PATH}/sra_process/par_s13_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s13_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s13_sra.txt"

# Sample 14
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_14_sra\"|" "${PARAMETER_PATH}/sra_process/par_s14_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s14_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s14_sra.txt"

# Sample 15
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_15_sra\"|" "${PARAMETER_PATH}/sra_process/par_s15_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s15_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s15_sra.txt"

# Sample 16
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_16_sra\"|" "${PARAMETER_PATH}/sra_process/par_s16_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s16_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s16_sra.txt"

# Sample 17
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_17_sra\"|" "${PARAMETER_PATH}/sra_process/par_s17_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s17_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s17_sra.txt"

# Sample 18
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_18_sra\"|" "${PARAMETER_PATH}/sra_process/par_s18_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s18_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s18_sra.txt"

# Sample 19
sed -i "14s|.*|output_path=\"${SRA_RESULTS_PATH}/sample_19_sra\"|" "${PARAMETER_PATH}/sra_process/par_s19_sra.txt"
sed -i "104s|.*|conda_bin=\"${CONDA_INST_DIR}\"|" "${PARAMETER_PATH}/sra_process/par_s19_sra.txt"
sed -i "105s|.*|conda_sh=\"${CONDA_SH_PATH}\"|" "${PARAMETER_PATH}/sra_process/par_s19_sra.txt"
