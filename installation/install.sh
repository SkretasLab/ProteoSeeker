#!/bin/bash
# Check if conda can be found.
if ! command -v conda &> /dev/null
then
    echo "conda was not found, installation was aborted"
    exit 1
fi
# Run all the installers.
./install_envs.sh
./install_dbs.sh
./parameter_files.sh
