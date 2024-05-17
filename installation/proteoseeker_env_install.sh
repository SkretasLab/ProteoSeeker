#!/bin/bash
# Base enviroment for ProteoSeeker with: Python, Openpyxl
conda create -n ps_env python=3.9.7 -y
source activate ps_env
python -m pip install openpyxl
python -c "import openpyxl"
if [ ! $? ]; then
  echo "The python package openpyxl was not installed successfully. Trying installation through anaconda."
  conda install anaconda::openpyxl -y
fi
python -c "import openpyxl"
if [ $? ]; then
  echo "The python package openpyxl was installed successfully."
else
  echo "The python package openpyxl was not installed successfully."
fi
conda deactivate
