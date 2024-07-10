#!/bin/bash
# Collect user input for downloading and extracting the nr database.
echo "Downloading Uniref90. The Uniref90 database is approximately 87.7 GBs in size (decompressed)."
# Determine whether the nr installer will run or not.
cd ..
if [ ! -d uniref_db ]; then
  mkdir uniref_db
fi
cd uniref_db
if [ ! -d uniref_90_db ]; then
  mkdir uniref_90_db
fi
cd uniref_90_db
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip -d uniref90.fasta.gz
cd ../../installation