#!/bin/bash
# Check if conda can be found.
if ! command -v conda &> /dev/null
then
    echo "conda was not found, installation was aborted"
    exit 1
fi
# Run Pfam installer.
./pfam_install.sh
# Run Swiss-Prot/UniprotKB installer.
./swissprot_install.sh
