#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python

cd /src/hpc-predict/hpc-predict-io/data_interfaces/

set -x
python3 validate_pipeline.py --pipeline flownet_data_interface.yml:inference cnn_segmenter_data_interface.yml:inference impact_data_interface.yml:data-assimilation
set +x

