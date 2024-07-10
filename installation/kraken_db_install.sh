#!/bin/bash
# In any case create a directory for the kraken database and download the minkraken database.
cd ..
if [ ! -d ps_tools ]; then
  mkdir ps_tools
fi
cd ps_tools
if [ ! -d kraken2 ]; then
  mkdir kraken2
fi
cd kraken2
if [ ! -d kraken_databases ]; then
  mkdir kraken_databases
fi
cd kraken_databases
# Download the minikraken 8GB database.
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
tar -xvzf k2_standard_08gb_20240112.tar.gz
rm k2_standard_08gb_20240112.tar.gz
cd ../../../installation