#!/bin/bash
set -e

# Collect user input for downloading and extracting the nr database.
read -p "Select whether the uniref100 protein database from the Uniprot database will be downloaded. Type y/yes for download or n/no otherwise. Caution! The Uniref100 database is approximately 96 GBs in size (decompressed). Type your selection: " uni_selection
# Determine whether the nr installer will run or not.
if [ $uni_selection = "y" ] || [ $uni_selection = "yes" ]; then
  # Get the current directory
  INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  # Get the parent directory which should be the main directory of ProteoSeeker.
  PS_DIR=$(dirname "$INSTALLATION_DIR")
  # Get the other directories.
  UNIREF_DBS_DIR="${PS_DIR}/uniref_db"
  UNIREF_DB_100_DIR="${UNIREF_DBS_DIR}/uniref_100_db"
  UNIREF_DB_100_GZ_FILE="${UNIREF_DB_100_DIR}/uniref100.fasta.gz"
  UNIREF_DB_100_FILE="${UNIREF_DB_100_DIR}/uniref100.fasta"

  # Collect user input for downloading and extracting the nr database.
  echo "Downloading Uniref100. The Uniref100 database is approximately 96 GBs in size (decompressed)."
  # Create the uniref_db dir if needed.
  if [ ! -d "${UNIREF_DBS_DIR}" ]; then
    mkdir "${UNIREF_DBS_DIR}"
  fi
  # Create the uniref_100_db dir if needed.
  if [ ! -d "${UNIREF_DB_100_DIR}" ]; then
    mkdir "${UNIREF_DB_100_DIR}"
  fi
  # Download.
  wget -P "${UNIREF_DB_100_DIR}" https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref_100/uniref_100.fasta.gz
  # Decompress.
  gunzip -dc "${UNIREF_DB_100_GZ_FILE}" > "${UNIREF_DB_100_FILE}"
  # Check for the database.
  if [[ -f "${UNIREF_DB_100_FILE}" ]]; then
      echo "Uniref100 was downloaded and extracted successfully."
  else
      echo "Uniref100 was not downloaded and extracted successfully."
  fi
elif [ $uni_selection = "n" ] || [ $uni_selection = "no" ]; then
  echo "Selected n/no. The nr protein database will not be downloaded."
else
  echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
