#!/bin/bash

set -euo pipefail

source /src/hpc-predict/flownet/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/flownet:/src/hpc-predict/hpc-predict-io/python

cd /src/hpc-predict/hpc-predict-io/data_interfaces/

set -x
python3 validate_data_interface.py --app flownet_data_interface.yml:inference
set +x

