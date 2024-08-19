#!/bin/bash
set -e

# Kraken2 all
python ps_br_analysis.py -p /mnt/4529bb0c-30cc-4e67-8f04-e94a1b226730/Works/Enzymes_Metagenomes/results_tax_eval -o analysis_plots_kraken_all -m k8,k16,k72,k8_0c1,k16_0c1,k72_0c1,k8_1c0,k16_1c0,k72_1c0,k8_10,k16_10,k72_10,k8_100,k16_100,k72_100,k8_ng,k16_ng,k72_ng -t False

# Kraken2 selected results and COMEBIn - nr and MetaBinner - nr results analysis
# Top method:
# kraken_8_100: Accuracy, F1 Score, Jaccard Index
# kraken_8_1c0: Precision
# kraken_72: Sensitivity
# kraken_72_ng: L1 Norm
# Selected methods: k8_100,k8_1c0,k72,k72_ng
python ps_br_analysis.py -p /mnt/4529bb0c-30cc-4e67-8f04-e94a1b226730/Works/Enzymes_Metagenomes/results_tax_eval -o analysis_plots_general -m k8_100,k8_1c0,k72,k72_ng,cnr,mnr
