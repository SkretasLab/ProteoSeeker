#!/bin/bash
set -e

# Set the path for the directory to be used as a mount point.
BIND_DIR="${HOME}/proteoseeker_bindmount"

# Run ProteoSeeker in interactive mode in docker.
read -p "Select which example will be run.\nChoose to run ProteoSeeker based on template file by selecting a number from 1, 2, 3, 4, 5, 6 and 7: " nr_selection
# Determine which example will be run.
if [ $nr_selection = "1" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_seek_tax_k_p.txt"
elif [ $nr_selection = "2" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_seek_tax_m_p.txt"
elif [ $nr_selection = "3" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_seek_tax_c_p.txt"
elif [ $nr_selection = "4" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_seek_p.txt"
elif [ $nr_selection = "5" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_tax_k_p.txt"
elif [ $nr_selection = "6" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_tax_m_p.txt"
elif [ $nr_selection = "7" ]; then
	sudo docker run -i -t --mount type=bind,src="${BIND_DIR}",target=/home/ps_data proteoseeker:main -pfp "/home/ps_data/parameter_files/par_tax_c_p.txt"
else
    echo "Improper selection. Exiting."
    exit 1
fi