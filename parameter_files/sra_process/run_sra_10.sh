#!/bin/bash
set -e

# Get the current directory
SRA_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARAMETER_PATH=$(dirname "${SRA_PATH}")
PS_PATH=$(dirname "${PARAMETER_PATH}")

if [ ! -d "${PS_PATH}/results/samples_sra" ]; then
  mkdir "${PS_PATH}/results/samples_sra"
fi

python -u "${PS_PATH}/proteoseeker.py" -pfp "${PS_PATH}/parameter_files/sra_process/par_s10_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s10_sra.txt"
