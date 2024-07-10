#!/bin/bash
# Check if conda can be found.
if ! command -v conda &> /dev/null
then
    echo "conda was not found, installation was aborted"
    exit 1
fi
# Create the parameter files
./cp_file_phylok_p.sh
./cp_file_phylomc_p.sh
./cp_seek_c.sh
./cp_seek_p.sh
./cp_file_seek_phylok_c.sh
./cp_file_seek_phylomc_c.sh
./cp_file_seek_phylok_p.sh
./cp_file_seek_phylomc_p.sh