#!/bin/bash

set -euo pipefail

source /src/hpc-predict/segmenter/random_walker_segmenter_for_mri_4d_flow/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/segmenter/random_walker_segmenter_for_mri_4d_flow:/src/hpc-predict/hpc-predict-io/python

cd /src/hpc-predict/hpc-predict-io/data_interfaces/

set -x
python3 validate_data_interface.py --app cnn_segmenter_data_interface.yml:inference
set +x

