#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python

cd /src/hpc-predict/hpc-predict-io/data_interfaces/

set -x
python3 validate_data_interface.py --app ../../IMPACT/data_interface/impact_data_interface.yml:data-assimilation
set +x

