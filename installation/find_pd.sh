#!/bin/bash
PD_PATH=""
if [ -d nr_database/nr ]; then
  PD_PATH="nr_database/nr"
elif [ -d uniref_db/uniref_50_db ]; then
  PD_PATH="uniref_db/uniref_50_db"
elif [ -d uniref_db/uniref_90_db ]; then
  PD_PATH="uniref_db/uniref_90_db"
elif [ -d uniref_db/uniref_100_db ]; then
  PD_PATH="uniref_db/uniref_100_db"
fi