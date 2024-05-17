#!/bin/bash
# Check if conda can be found.
if ! command -v conda &> /dev/null
then
    echo "conda was not found, installation was aborted"
    exit 1
fi
# Run tool installers.
echo -e "\nInstalling SRA tools..."
./sra_tools_install.sh &> output_errors/sra_tools_stdoe.txt
echo -e "\nInstalling DIAMOND BLAST..."
./diamond_install.sh &> output_errors/diamond_blast_stdoe.txt
echo -e "\nInstalling FastQC..."
./fastqc_install.sh &> output_errors/fastqc_stdoe.txt
echo -e "\nInstalling BBDuk..."
./bbduk_install.sh &> output_errors/bbduk_stdoe.txt
echo -e "\nInstalling Megahit..."
./megahit_install.sh &> output_errors/megahit_stdoe.txt
echo -e "\nInstalling kraken2..."
./kraken_install.sh &> output_errors/kraken_stdoe.txt
echo -e "\nInstalling COMEBin..."
./comebin_install.sh &> output_errors/comebin_stdoe.txt
echo -e "\nInstalling MetaBinner..."
./metabinner_install.sh &> output_errors/metabinner_stdoe.txt
echo -e "\nInstalling CD-HIT..."
./cd_hit_install.sh &> output_errors/cd_hit_stdoe.txt
echo -e "\nInstalling HMMER..."
./hmmer_install.sh &> output_errors/hmmer_stdoe.txt
echo -e "\nInstalling Diamond..."
./diamond_install.sh &> output_errors/diamond_stdoe.txt
echo -e "\nInstalling FragGeneScanRs..."
./fraggenescanrs_install.sh &> output_errors/fraggenescanrs_stdoe.txt
echo -e "\nInstalling TaxonKit..."
./taxonkit.sh &> output_errors/taxonkit_stdoe.txt
echo -e "\nInstalling Bowtie2..."
./bowtie.sh &> output_errors/bowtie_stdoe.txt
echo -e "\nInstalling Phobius enviroment..."
./phobius_install.sh &> output_errors/phobius.txt
echo -e "\nInstalling the ProteoSeeker base enviroment..."
./proteoseeker_env_install.sh &> output_errors/proteoseeker_env_stdoe.txt
