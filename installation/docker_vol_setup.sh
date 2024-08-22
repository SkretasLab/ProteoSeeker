#!/bin/bash
set -e

# Set the name for the docker volume.
PS_VOL_NAME="proteoseeker_vol"

# Create the Docker volume
docker volume create "${PS_VOL_NAME}"

# Find the full path of the Docker volume.
VOL_PATH=$(docker volume inspect "${PS_VOL_NAME}" | grep -oP '"Mountpoint":\s*"\K[^"]+')

# Create the directories needed to run ProteoSeeker in the Docker image.
mkdir "${VOL_PATH}/input_files"
mkdir "${VOL_PATH}/protein_databases"
mkdir "${VOL_PATH}/results"
mkdir "${VOL_PATH}/parameter_files"
mkdir "${VOL_PATH}/phobius"

# Copy the template parameter files used to run ProteoSeeker in the Docker image.
cp "parameter_files/docker/par_tax_m_p.txt" \
"parameter_files/docker/par_tax_k_p.txt" \
"parameter_files/docker/par_tax_c_p.txt" \
"parameter_files/docker/par_seek_tax_m_p.txt" \
"parameter_files/docker/par_seek_tax_k_p.txt" \
"parameter_files/docker/par_seek_tax_c_p.txt" \
"parameter_files/docker/par_seek_p.txt" \
"${VOL_PATH}/parameter_files"

# Copy the protein database used to run ProteoSeeker in the Docker image.
cp "parameter_files/docker/nr_part.fasta" "${VOL_PATH}/protein_databases"