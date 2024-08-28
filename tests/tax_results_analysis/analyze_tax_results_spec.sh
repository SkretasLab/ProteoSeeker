#!/bin/bash
set -e

# Get the current directory
TAX_ANALYSIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the tests directory of ProteoSeeker.
TESTS_DIR=$(dirname "${TAX_ANALYSIS_DIR}")
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${TAX_ANALYSIS_DIR}")
# Set the path of the directory with the results from the taxonomy evaluation.
TAX_RESULTS_DIR="${PS_DIR}/results_tax_eval"

# Kraken2 8
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_db8 -m k8,k8_0c1,k8_1c0,k8_10,k8_100,k8_ng -t False

# Kraken2 16
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_db16 -m k16,k16_0c1,k16_1c0,k16_10,k16_100,k16_ng -t False

# Kraken2 72
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_db72 -m k72,k72_0c1,k72_1c0,k72_10,k72_100,k72_ng -t False

# Kraken2 0.0%
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_f0c0 -m k8,k16,k72 -t False

# Kraken2 0.1%
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_f0c1 -m k8_0c1,k16_0c1,k72_0c1 -t False

# Kraken2 1.0%
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_f1c0 -m k8_1c0,k16_1c0,k72_1c0 -t False

# Kraken2 10
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_f10 -m k8_10,k16_10,k72_10 -t False

# Kraken2 100
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_f100 -m k8_100,k16_100,k72_100 -t False

# Kraken2 non-gut
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_fng -m k8_ng,k16_ng,k72_ng -t False
