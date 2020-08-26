#!/bin/bash

set -euo pipefail

HPC_PREDICT_IO_IMAGE=${HPC_PREDICT_IO_IMAGE:-'lukasgd/hpc-predict:io-deploy'}

# visualize result - here shown on the example of flownet, but 
INFERENCE_OUTPUT_FILE=$(realpath $1)

INFERENCE_OUTPUT_DIR=$(dirname ${INFERENCE_OUTPUT_FILE})

host_input_output_directory="${INFERENCE_OUTPUT_DIR}"

if [ -f ${host_input_output_directory} ]; then
    echo "Directory \"${host_input_output_directory}\" with HPC-PREDICT-IO HDF5 input data for visualization does not exist. Exiting..."
    exit 1
fi

echo "Host output directory: ${host_input_output_directory}"

# FIXME: Why is there hpc-predict-output and inference-output
shell_command=$(printf "%s" \
    " source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
    " export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python && " \
    " set -x && " \
    " python3 " \
    " /src/hpc-predict/hpc-predict-io/python/hdf5xdmf.py " \
    " ${INFERENCE_OUTPUT_FILE} ") 

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${INFERENCE_OUTPUT_DIR}:${INFERENCE_OUTPUT_DIR} --entrypoint bash "${HPC_PREDICT_IO_IMAGE}" -c "${shell_command}"
set +x
