#!/bin/bash
set -e

# Collect user input for downloading and extracting the nr database.
read -p "Select whether the uniref90 protein database from the Uniprot database will be downloaded. Type y/yes for download or n/no otherwise. Caution! The Uniref90 database is approximately 87.7 GBs in size (decompressed). Type your selection: " uni_selection
# Determine whether the nr installer will run or not.
if [ $uni_selection = "y" ] || [ $uni_selection = "yes" ]; then
  # Get the current directory
  INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  # Get the parent directory which should be the main directory of ProteoSeeker.
  PS_DIR=$(dirname "$INSTALLATION_DIR")
  # Get the other directories.
  UNIREF_DBS_DIR="${PS_DIR}/uniref_db"
  UNIREF_DB_90_DIR="${UNIREF_DBS_DIR}/uniref_90_db"
  UNIREF_DB_90_GZ_FILE="${UNIREF_DB_90_DIR}/uniref90.fasta.gz"
  UNIREF_DB_90_FILE="${UNIREF_DB_90_DIR}/uniref90.fasta"

  # Collect user input for downloading and extracting the nr database.
  echo "Downloading Uniref90. The Uniref90 database is approximately 87.7 GBs in size (decompressed)."
  # Create the uniref_db dir if needed.
  if [ ! -d "${UNIREF_DBS_DIR}" ]; then
    mkdir "${UNIREF_DBS_DIR}"
  fi
  # Create the uniref_90_db dir if needed.
  if [ ! -d "${UNIREF_DB_90_DIR}" ]; then
    mkdir "${UNIREF_DB_90_DIR}"
  fi
  # Download.
  wget -P "${UNIREF_DB_90_DIR}" https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
  # Decompress.
  gunzip -dc "${UNIREF_DB_90_GZ_FILE}" > "${UNIREF_DB_90_FILE}"
  # Check for the database.
  if [[ -f "${UNIREF_DB_90_FILE}" ]]; then
      echo "Uniref90 was downloaded and extracted successfully."
  else
      echo "Uniref90 was not downloaded and extracted successfully."
  fi
elif [ $uni_selection = "n" ] || [ $uni_selection = "no" ]; then
  echo "Selected n/no. The nr protein database will not be downloaded."
else
  echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
