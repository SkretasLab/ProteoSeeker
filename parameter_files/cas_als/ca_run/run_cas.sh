#!/bin/bash
set -e

# DRR163688 - Kraken2 non-gut:
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_betaca_ca_sra_dbs.txt" &> "results_cas_als/stdoe_sra_dbs_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_cnr_ca.txt" &> "results_cas_als/stdoe_cnr_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_mnr_ca.txt" &> "results_cas_als/stdoe_mnr_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k8_ca.txt" &> "results_cas_als/stdoe_k8_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k16_ca.txt" &> "results_cas_als/stdoe_k16_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k72_ca.txt" &> "results_cas_als/stdoe_k72_DRR163688.txt"

# DRR163688 - Kraken2 0.01%:
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k8_0c01_ca.txt" &> "results_cas_als/stdoe_k8_0c01_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k16_0c01_ca.txt" &> "results_cas_als/stdoe_k16_0c01_DRR163688.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/DRR163688/parameters_DRR163688_k72_0c01_ca.txt" &> "results_cas_als/stdoe_k72_0c01_DRR163688.txt"

# SRR3961740 - Kraken2 non-gut:
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_sra.txt" &> "results_cas_als/stdoe_sra_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_cnr_ca.txt" &> "results_cas_als/stdoe_cnr_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_mnr_ca.txt" &> "results_cas_als/stdoe_mnr_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k8_ca.txt" &> "results_cas_als/stdoe_k8_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k16_ca.txt" &> "results_cas_als/stdoe_k16_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k72_ca.txt" &> "results_cas_als/stdoe_k72_SRR3961740.txt"

# SRR3961740 - Kraken2 0.01%:
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k8_0c01_ca.txt" &> "results_cas_als/stdoe_k8_0c01_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k16_0c01_ca.txt" &> "results_cas_als/stdoe_k16_0c01_SRR3961740.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/ca_run/SRR3961740/parameters_SRR3961740_k72_0c01_ca.txt" &> "results_cas_als/stdoe_k72_0c01_SRR3961740.txt"

