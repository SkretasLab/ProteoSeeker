#!/bin/bash
set -e

# Alpha-amylase profile & amylase protein databases
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/parameters_alphaal_al_dbs.txt" &> "results_cas_als/stdoe_dbs_al.txt"

# SRR17771278
# SRA
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_sra.txt" &> "results_cas_als/stdoe_sra_al.txt"
# Runs - Kraken2 - non-gut
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_cnr_al.txt" &> "results_cas_als/stdoe_cnr_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_mnr_al.txt" &> "results_cas_als/stdoe_mnr_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k8_al.txt" &> "results_cas_als/stdoe_k8_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k16_al.txt" &> "results_cas_als/stdoe_k16_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k72_al.txt" &> "results_cas_als/stdoe_k72_al.txt"
# The HotZyme is not run because the taxonomy route of COMEBin/MetaBinner can not be applied when the reads are not available to ProteoSeeker. The taxonomy route of Kraken2 ideally should need reads but could be applied to contigs. Both routes should be able to be applied to compare them, hence none of the two is eventually applied.

# Runs - Kraken2 - 0.01%
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k8_0c01_al.txt" &> "results_cas_als/stdoe_k8_0c01_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k16_0c01_al.txt" &> "results_cas_als/stdoe_k16_0c01_al.txt"
python -u proteoseeker.py -pfp "parameter_files/cas_als/al_run/SRR17771278/parameters_SRR17771278_k72_0c01_al.txt" &> "results_cas_als/stdoe_k72_0c01_al.txt"

