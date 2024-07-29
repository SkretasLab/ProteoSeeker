#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
KRAKEN2_DIR="${PS_TOOLS_DIR}/kraken2"
KRAKEN2_DBS_DIR="${KRAKEN2_DIR}/kraken2_databases"
KRAKEN2_8_ST_DB_DIR="${KRAKEN2_DBS_DIR}/kraken2_8st_db"
KRAKEN2_8_ST_DB_FILE="${KRAKEN2_8_ST_DB_DIR}/k2_standard_08gb_20240112.tar.gz"

# Create the ps_tools dir if needed.
if [ ! -d "${PS_TOOLS_DIR}" ]; then
  mkdir "${PS_TOOLS_DIR}"
fi
# Create the kraken2 dir if needed.
if [ ! -d "${KRAKEN2_DIR}" ]; then
  mkdir "${KRAKEN2_DIR}"
fi
# Create the kraken_databases dir if needed.
if [ ! -d "${KRAKEN2_DBS_DIR}" ]; then
  mkdir "${KRAKEN2_DBS_DIR}"
fi
# Create the kraken2_8st_db dir if needed.
if [ ! -d "${KRAKEN2_8_ST_DB_DIR}" ]; then
  mkdir "${KRAKEN2_8_ST_DB_DIR}"
fi
# Download.
wget -P "${KRAKEN2_DBS_DIR}" https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz
# Decompress.
tar -xvzf "${KRAKEN2_8_ST_DB_FILE}" -C "${KRAKEN2_8_ST_DB_DIR}"