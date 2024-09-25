#!/bin/bash
set -e

# Get the current directory
TAX_ANALYSIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the tests directory of ProteoSeeker.
TESTS_DIR=$(dirname "${TAX_ANALYSIS_DIR}")
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${TESTS_DIR}")
# Get the path of the installation directory.
INSTALLATION_DIR="${PS_DIR}/installation"
# Set the path of the directory with the results from the taxonomy evaluation.
TAX_RESULTS_DIR="${PS_DIR}/results_tax_eval"

# Source conda.
source "${INSTALLATION_DIR}/find_conda.sh"
source $CONDA_SH_PATH

# Activate environment.
conda activate ps_result_analysis

# Kraken2 all
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_all -m k8,k16,k77,k8_0c01,k16_0c01,k77_0c01,k8_0c1,k16_0c1,k77_0c1,k8_1c0,k16_1c0,k77_1c0,k8_5c0,k16_5c0,k77_5c0,k8_10,k16_10,k77_10,k8_100,k16_100,k77_100,k8_500,k16_500,k77_500,k8_ng,k16_ng,k77_ng,k8_g,k16_g,k77_g -t False

# Selected combinations:
# k8_ng: Accuracy, F1 Score, Jaccard Index
# k8_500: F1 Score, L1 Norm
# k16_5c0: False Positive, Precision
# k77_5c0: False Positive, Precision
# k77_10: False Negative, L1 Norm, True Positive
# k8_500,k8_ng,k16_5c0,k77_5c0,k77_10
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_general -m k8_500,k8_ng,k16_5c0,k77_5c0,k77_10,cnr,mnr

# Deactivate environment.
conda deactivate