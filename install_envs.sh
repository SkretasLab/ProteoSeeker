#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "${INSTALLATION_DIR}/find_conda.sh"
# Run tool installers.
echo -e "\nInstalling SRA tools..."
"${INSTALLATION_DIR}/sra_tools_install.sh &> output_errors/sra_tools_stdoe.txt"
echo -e "\nInstalling DIAMOND BLAST..."
"${INSTALLATION_DIR}/diamond_install.sh &> output_errors/diamond_blast_stdoe.txt"
echo -e "\nInstalling FastQC..."
"${INSTALLATION_DIR}/fastqc_install.sh &> output_errors/fastqc_stdoe.txt"
echo -e "\nInstalling BBDuk..."
"${INSTALLATION_DIR}/bbduk_install.sh &> output_errors/bbduk_stdoe.txt"
echo -e "\nInstalling Megahit..."
"${INSTALLATION_DIR}/megahit_install.sh &> output_errors/megahit_stdoe.txt"
echo -e "\nInstalling kraken2..."
"${INSTALLATION_DIR}/kraken_install.sh &> output_errors/kraken_stdoe.txt"
echo -e "\nInstalling COMEBin..."
"${INSTALLATION_DIR}/comebin_install.sh &> output_errors/comebin_stdoe.txt"
echo -e "\nInstalling MetaBinner..."
"${INSTALLATION_DIR}/metabinner_install.sh &> output_errors/metabinner_stdoe.txt"
echo -e "\nInstalling CD-HIT..."
"${INSTALLATION_DIR}/cd_hit_install.sh &> output_errors/cd_hit_stdoe.txt"
echo -e "\nInstalling HMMER..."
"${INSTALLATION_DIR}/hmmer_install.sh &> output_errors/hmmer_stdoe.txt"
echo -e "\nInstalling FragGeneScanRs..."
"${INSTALLATION_DIR}/fraggenescanrs_install.sh &> output_errors/fraggenescanrs_stdoe.txt"
echo -e "\nInstalling TaxonKit..."
"${INSTALLATION_DIR}/taxonkit.sh &> output_errors/taxonkit_stdoe.txt"
echo -e "\nInstalling Bowtie2..."
"${INSTALLATION_DIR}/bowtie.sh &> output_errors/bowtie_stdoe.txt"
echo -e "\nCreating the Phobius enviroment..."
"${INSTALLATION_DIR}/phobius_install.sh &> output_errors/phobius_stdoe.txt"
echo -e "\nInstalling the ProteoSeeker base enviroment..."
"${INSTALLATION_DIR}/proteoseeker_env_install.sh &> output_errors/proteoseeker_env_stdoe.txt"
