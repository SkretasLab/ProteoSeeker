#!/bin/bash
# Collect user input for downloading and extracting the nr database.
read -p "Select whether the uniref50 protein database from the Uniprot database will be downloaded. Type y/yes for download or n/no otherwise. Caution! The Uniref50 database is approximately 26.1 GBs in size (decompressed). Type your selection: " uni_selection
# Determine whether the nr installer will run or not.
if [ $uni_selection = "y" ] || [ $uni_selection = "yes" ]; then
  echo "Selected y/yes."
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
  cd ../..
  cd installation
elif [ $uni_selection = "n" ] || [ $uni_selection = "no" ]; then
  echo "Selected n/no. The nr protein database will not be downloaded."
else
  echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
