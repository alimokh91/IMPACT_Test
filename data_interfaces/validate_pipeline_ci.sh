#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python

cd /src/hpc-predict/hpc-predict-io/data_interfaces/

set -x
python3 validate_pipeline.py --pipeline ../../flownet/data_interface/flownet_data_interface.yml:inference ../../segmenter/data_interface/cnn_segmenter_data_interface.yml:inference ../../IMPACT/data_interface/impact_data_interface.yml:data-assimilation
set +x

