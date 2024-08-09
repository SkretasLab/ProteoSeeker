#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
NR_DB_DIR="${PS_DIR}/nr_database"
NR_DB_GZ_FILE="${NR_DB_DIR}/nr.gz"
NR_DB_FILE="${NR_DB_DIR}/nr"

# Collect user input for downloading and extracting the nr database.
read -p "Select whether the nr protein database from NCBI will be downloaded and extracted. Type y/yes for download and extraction or n/no otherwise. Caution! The extracted nr database is approximately 300 GBs in size. Type your selection: " nr_selection
# Determine whether the nr installer will run or not.
if [ $nr_selection = "y" ] || [ $nr_selection = "yes" ]; then
    echo "Selected y/yes."
    # Create the nr_database dir if needed.
    if [ ! -d "${NR_DB_DIR}" ]; then
        mkdir "${NR_DB_DIR}"
    fi
    # Download.
    wget -P "${NR_DB_DIR}" https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    # Decompress the database.
    gunzip -dc "${NR_DB_GZ_FILE}" > "${NR_DB_FILE}"
elif [ $nr_selection = "n" ] || [ $nr_selection = "no" ]; then
    echo "Selected n/no. The nr protein database will not be downloaded."
else
    echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
