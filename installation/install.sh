#!/bin/bash
source ./find_conda.sh
# Run the installers.
./install_envs.sh
./install_dbs.sh
./parameter_files.sh