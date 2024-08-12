#!/bin/bash
set -e

# Get the current directory
SRA_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARAMETER_PATH=$(dirname "${SRA_PATH}")
PS_PATH=$(dirname "${PARAMETER_PATH}")

if [ ! -d "${PS_PATH}/results/samples_sra" ]; then
  mkdir "${PS_PATH}/results/samples_sra"
fi

python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s1_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s1_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s2_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s2_sea.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s3_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s3_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s4_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s4_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s5_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s5_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s6_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s6_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s7_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s7_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s8_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s8_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s9_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s9_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s10_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s10_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s11_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s11_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s12_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s12_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s13_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s13_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s14_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s14_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s15_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s15_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s16_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s16_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s17_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s17_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s18_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s18_sra.txt"
python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s19_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s19_sra.txt"
