#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
NR_DB_FILE="${PS_DIR}/nr_database/nr"
UNI50_DB_FILE="${PS_DIR}/uniref_db/uniref_50_db"
UNI90_DB_FILE="${PS_DIR}/uniref_db/uniref_90_db"
UNI100_DB_FILE="${PS_DIR}/uniref_db/uniref_100_db"

PD_PATH=""
if [ -d "${NR_DB_FILE}" ]; then
  PD_PATH="nr_database/nr"
elif [ -d "${UNI50_DB_FILE}" ]; then
  PD_PATH="uniref_db/uniref_50_db"
elif [ -d "${UNI90_DB_FILE}" ]; then
  PD_PATH="uniref_db/uniref_90_db"
elif [ -d "${UNI100_DB_FILE}" ]; then
  PD_PATH="uniref_db/uniref_100_db"
fi