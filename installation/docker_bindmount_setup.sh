#!/bin/bash
set -e

# Set the path for the directory to be used as a mount point.
BIND_DIR="${HOME}/proteoseeker_bindmount"

# Create the directory.
mkdir "${BIND_DIR}"

# Create the directories needed to run ProteoSeeker in the Docker image.
mkdir "${BIND_DIR}/input_files"
mkdir "${BIND_DIR}/protein_databases"
mkdir "${BIND_DIR}/results"
mkdir "${BIND_DIR}/parameter_files"
mkdir "${BIND_DIR}/phobius"

# Copy the template parameter files used to run ProteoSeeker in the Docker image.
cp "parameter_files/docker/par_seek_p.txt" \
"parameter_files/docker/par_seek_tax_mc_c.txt" \
"parameter_files/docker/par_tax_k_p.txt" \
"parameter_files/docker/par_seek_tax_k_p.txt" \
"parameter_files/docker/par_seek_tax_mc_p.txt" \
"parameter_files/docker/par_tax_mc_p.txt" \
"${BIND_DIR}/parameter_files"

# Copy the protein database used to run ProteoSeeker in the Docker image.
cp "parameter_files/docker/nr_head_1000000_rnapol.fasta" "${BIND_DIR}/protein_databases"