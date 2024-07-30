#!/bin/bash
set -e

MTBPATH=$(conda run -n ps_metabinner which run_metabinner.sh)
MBINPATH="$(dirname "${MTBPATH}")" ; FILE="$(basename "${MTBPATH}")"
CTBPATH=$(conda run -n ps_comebin which run_comebin.sh)
CBINPATH="$(dirname "${CTBPATH}")" ; FILE="$(basename "${CTBPATH}")"