#!/bin/bash
set -e

./taxonomy_tests/setup.sh

if [ ! -d "${PS_PATH}/results/samples_sra" ]; then
  mkdir "${PS_PATH}/results/samples_sra"
fi

python -u proteoseeker.py -pfp "${PS_PATH}/parameter_files/sra_process/par_s15_sra.txt" &> "${PS_PATH}/results/samples_sra/stdoe_s15_sra.txt"