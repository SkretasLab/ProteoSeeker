#!/bin/bash
set -e

python ps_caal_analysis.py -p /mnt/4529bb0c-30cc-4e67-8f04-e94a1b226730/Works/Enzymes_Metagenomes/results_seek_tax_eval -b proteins_blastp_nr.tsv -o analysis_results
