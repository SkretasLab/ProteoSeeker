#!/bin/bash
# TaxonKit
source ./find_conda.sh
source $CONDA_SH_PATH
conda create -n ps_taxonkit -y
conda activate ps_taxonkit
conda install -c bioconda taxonkit -y
conda install bioconda::csvtk -y
conda deactivate
# Download the necessary information from the Taxonomy database
mkdir taxonomy_db_info
cd taxonomy_db_info
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit
cd ..