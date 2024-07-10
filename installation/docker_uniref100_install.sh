#!/bin/bash
# Collect user input for downloading and extracting the nr database.
echo "Downloading Uniref100. The Uniref100 database is approximately 96 GBs in size (decompressed)."
# Determine whether the nr installer will run or not.
cd ..
if [ ! -d uniref_db ]; then
  mkdir uniref_db
fi
cd uniref_db
if [ ! -d uniref_50_db ]; then
  mkdir uniref_100_db
fi
cd uniref_100_db
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref_100/uniref_100.fasta.gz
gunzip -d uniref_100.fasta.gz
cd ../../installation