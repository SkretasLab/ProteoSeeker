#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the parent directory which should be the main directory of ProteoSeeker.
PS_DIR=$(dirname "${INSTALLATION_DIR}")
# Get the other directories.
PS_TOOLS_DIR="${PS_DIR}/ps_tools"
FGSRS_DIR="${PS_TOOLS_DIR}/fgsrs"
FGSRS_TAR_GZ_FILE="${FGSRS_DIR}/FragGeneScanRs-v1.1.0-x86_64-unknown-linux-musl.tar.gz"
FGSRS_FILE="${FGSRS_DIR}/FragGeneScanRs"

# Create the ps_tools dir if needed.
if [ ! -d "${PS_TOOLS_DIR}" ]; then
    mkdir "${PS_TOOLS_DIR}"
fi
# Create the fgsrs dir if needed.
if [ ! -d "${FGSRS_DIR}" ]; then
    mkdir "${FGSRS_DIR}"
fi
# Download a specific version.
wget -P "${FGSRS_DIR}" https://github.com/unipept/FragGeneScanRs/releases/download/v1.1.0/FragGeneScanRs-v1.1.0-x86_64-unknown-linux-musl.tar.gz
# Decompress.
tar -xvzf "${FGSRS_TAR_GZ_FILE}" -C "${FGSRS_DIR}"
# Check for the executable.
if [[ -f "${FGSRS_FILE}" ]]; then
    echo "FragGeneScanRs was installed successfully."
else
    echo "FragGeneScanRs was not installed successfully."
    exit 1
fi