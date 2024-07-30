#!/bin/bash
set -e

# Get the current directory
INSTALLATION_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check if conda can be found.
if ! command -v conda &> /dev/null
then
    echo "conda was not found, installation was aborted"
    exit 1
fi
# Run Pfam installer.
"${INSTALLATION_DIR}/pfam_install.sh"
# Run Swiss-Prot/UniprotKB installer.
"${INSTALLATION_DIR}/swissprot_install.sh"
# Run Kraken2 database installer.
"${INSTALLATION_DIR}/kraken_db_install.sh"