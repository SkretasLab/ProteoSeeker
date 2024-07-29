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
"${INSTALLATION_DIR}/par_seek_p.sh"
"${INSTALLATION_DIR}/par_seek_c.sh"
"${INSTALLATION_DIR}/par_tax_k_p.sh"
"${INSTALLATION_DIR}/par_tax_mc_p.sh"
"${INSTALLATION_DIR}/par_seek_tax_k_p.sh"
"${INSTALLATION_DIR}/par_seek_tax_mc_p.sh"
"${INSTALLATION_DIR}/par_seek_tax_k_c.sh"
"${INSTALLATION_DIR}/par_seek_tax_mc_c.sh"