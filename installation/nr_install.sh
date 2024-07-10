#!/bin/bash
# Collect user input for downloading and extracting the nr database.
read -p "Select whether the nr protein database from NCBI will be downloaded and extracted. Type y/yes for download and extraction or n/no otherwise. Caution! The extracted nr database is approximately 300 GBs in size. Type your selection: " nr_selection
# Determine whether the nr installer will run or not.
if [ $nr_selection = "y" ] || [ $nr_selection = "yes" ]; then
  echo "Selected y/yes."
  # Moving to the main folder.
  cd ..
  # nr database
  mkdir nr_database
  cd nr_database
  wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
  gunzip -d nr.gz
  cd ..
elif [ $nr_selection = "n" ] || [ $nr_selection = "no" ]; then
  echo "Selected n/no. The nr protein database will not be downloaded."
else
  echo "Improper selection. Using default selection (n/no) which is not to download the nr database."
fi
