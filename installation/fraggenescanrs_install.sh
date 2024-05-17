#!/bin/bash
# Moving to the main folder.
cd ..
if [ ! -d ps_tools ]; then
  mkdir ps_tools
fi
# FragGeneScanRs
cd ps_tools
mkdir fgsrs
cd fgsrs
wget https://github.com/unipept/FragGeneScanRs/releases/download/v1.1.0/FragGeneScanRs-v1.1.0-x86_64-unknown-linux-musl.tar.gz
tar -xvzf FragGeneScanRs-v1.1.0-x86_64-unknown-linux-musl.tar.gz
cd ../..
