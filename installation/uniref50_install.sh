#!/bin/bash
set -e

# Collect user input for downloading and extracting the nr database.
read -p "Select whether the uniref50 protein database from the Uniprot database will be downloaded. Type y/yes for download or n/no otherwise. Caution! The Uniref50 database is approximately 26.1 GBs in size (decompressed). Type your selection: " uni_selection
# Determine whether the nr installer will run or not.
if [ $uni_selection = "y" ] || [ $uni_selection = "yes" ]; then
  echo "Selected y/yes."
  # Get the current directory
  INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  # Get the parent directory which should be the main directory of ProteoSeeker.
  PS_DIR=$(dirname "$INSTALLATION_DIR")
  # Get the other directories.
  UNIREF_DBS_DIR="${PS_DIR}/uniref_db"
  UNIREF_DB_50_DIR="${UNIREF_DBS_DIR}/uniref_50_db"
  UNIREF_DB_50_GZ_FILE="${UNIREF_DB_50_DIR}/uniref50.fasta.gz"
  UNIREF_DB_50_FILE="${UNIREF_DB_50_DIR}/uniref50.fasta"

  # Collect user input for downloading and extracting the nr database.
  echo "Downloading Uniref50. The Uniref50 database is approximately 26.1 GBs in size (decompressed)."
  # Create the uniref_db dir if needed.
  if [ ! -d "${UNIREF_DBS_DIR}" ]; then
    mkdir "${UNIREF_DBS_DIR}"
  fi
  # Create the uniref_50_db dir if needed.
  if [ ! -d "${UNIREF_DB_50_DIR}" ]; then
    mkdir "${UNIREF_DB_50_DIR}"
  fi
  # Download.
  wget -P "${UNIREF_DB_50_DIR}" https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
  # Decompress.
  gunzip -dc "${UNIREF_DB_50_GZ_FILE}" > "${UNIREF_DB_50_FILE}"
  # Check for the database.
  if [[ -f "${UNIREF_DB_50_FILE}" ]]; then
      echo "Uniref50 was downloaded and extracted successfully."
  else
      echo "Uniref50 was not downloaded and extracted successfully."
  fi
elif [ $uni_selection = "n" ] || [ $uni_selection = "no" ]; then
  echo "Selected n/no. The nr protein database will not be downloaded."
else
  echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
