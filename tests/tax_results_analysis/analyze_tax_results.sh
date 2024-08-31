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
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_kraken_all -m k8,k16,k72,k8_0c1,k16_0c1,k72_0c1,k8_1c0,k16_1c0,k72_1c0,k8_10,k16_10,k72_10,k8_100,k16_100,k72_100,k8_ng,k16_ng,k72_ng -t False

# Kraken2 selected results and COMEBIn - nr and MetaBinner - nr results analysis
# Top method:
# kraken_8_100: Accuracy, F1 Score, Jaccard Index
# kraken_8_1c0: False Positive, Precision
# kraken_16_1c0: Precision
# kraken_72: True Positive, False Negative, Sensitivity
# kraken_72_ng: L1 Norm
# Selected methods: k8_100,k8_1c0,k72,k72_ng
python ps_br_analysis.py -p "${TAX_RESULTS_DIR}" -o analysis_plots_general -m k8_100,k8_1c0,k16_1c0,k72,k72_ng,cnr,mnr

# Deactivate environment.
conda deactivate