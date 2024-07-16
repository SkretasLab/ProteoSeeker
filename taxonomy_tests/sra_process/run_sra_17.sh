#!/bin/bash
set -e

./taxonomy_tests/setup.sh

if [ ! -d "${PS_PATH}/results/samples_sra" ]; then
  mkdir "${PS_PATH}/results/samples_sra"
fi

python -u proteoseeker.py -pfp "${PS_PATH}/parameter_files/sra_process/par_s17_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s17_sra.txt"