#!/bin/bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Wrong number of arguments "$#" - expected data-assimilation-input, data-assimilation-output."
    exit 1
fi

container_input_file="$1"
container_output_file="$2"

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate
export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python

# Generate config file
# FIXME: Move input configuration on commandline here to a tracked file
set -x
python3 -u \
    -m mr_io_impact_config \
    --input-mri "${container_input_file}" \
    --output-mri "${container_output_file}" \
    --sr 2 2 2 \
    --padding 0.5 0.5 0.5 \
    --tr 2 \
    --output "$(dirname ${container_output_file})/config.txt" \
    --np 4
set +x
