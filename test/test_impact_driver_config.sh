#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test
export PATH=${FORTRAN_TEST_BINARY_PATH}:${PATH}

# This must be passed to the container in an environment variable
# export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE}

HPC_PREDICT_IO_TEST_CONFIG_COMMAND=$(python /src/hpc-predict/hpc-predict-io/test/test_impact_args.py  --component config --test "$@")

echo "HPC_PREDICT_IO_TEST_CONFIG_COMMAND = ${HPC_PREDICT_IO_TEST_CONFIG_COMMAND}"

IFS='.' read -r test_module test_class <<<"$1"
cd ${CI_CACHE_FOLDER}/${test_module}/${test_class}

set -x
echo "$(pwd)"
eval ${HPC_PREDICT_IO_TEST_CONFIG_COMMAND}
set +x

