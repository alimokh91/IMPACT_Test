#!/bin/bash

set -euxo pipefail

HPC_PREDICT_DATA_DIR="$(realpath "$1")"
HPC_PREDICT_MRI_PATH="$(realpath --relative-to="$1" "$2")"

# Development of visualization notebook in Docker
docker run --rm -v "${HPC_PREDICT_DATA_DIR}":/hpc-predict-data -v $(pwd)/../python/visualization:/src/hpc-predict/hpc-predict-io/python/visualization --env HPC_PREDICT_MRI_PATH="${HPC_PREDICT_MRI_PATH}" --network host --entrypoint bash ${HPC_PREDICT_IO_IMAGE} -c "export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python; source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && jupyter lab --allow-root /src/hpc-predict/hpc-predict-io/python/visualization/flow_mri_viewer.ipynb"

# Running visualization notebook in Docker
#docker run --rm -v "${HPC_PREDICT_DATA_DIR}":/hpc-predict-data --env HPC_PREDICT_MRI_PATH="${HPC_PREDICT_MRI_PATH}" --network host --entrypoint bash ${HPC_PREDICT_IO_IMAGE} -c "export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python; source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && jupyter lab --allow-root /src/hpc-predict/hpc-predict-io/python/visualization/flow_mri_viewer.ipynb"
