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
# Create the parameter files
"${INSTALLATION_DIR}/cp_file_phylok_p.sh"
"${INSTALLATION_DIR}/cp_file_phylomc_p.sh"
"${INSTALLATION_DIR}/cp_seek_c.sh"
"${INSTALLATION_DIR}/cp_seek_p.sh"
"${INSTALLATION_DIR}/cp_file_seek_phylok_c.sh"
"${INSTALLATION_DIR}/cp_file_seek_phylomc_c.sh"
"${INSTALLATION_DIR}/cp_file_seek_phylok_p.sh"
"${INSTALLATION_DIR}/cp_file_seek_phylomc_p.sh"