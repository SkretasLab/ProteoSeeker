#!/bin/bash
# Collect user input for downloading and extracting the nr database.
echo "Downloading Uniref50. The Uniref50 database is approximately 26.1 GBs in size (decompressed)."
# Determine whether the nr installer will run or not.
cd ..
if [ ! -d uniref_db ]; then
  mkdir uniref_db
fi
cd uniref_db
if [ ! -d uniref_50_db ]; then
  mkdir uniref_50_db
fi
cd uniref_50_db
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
gunzip -d uniref50.fasta.gz
cd ../../installation