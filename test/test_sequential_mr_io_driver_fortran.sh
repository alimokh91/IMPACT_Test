#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

export FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test

export HPC_PREDICT_IO_TEST_FORTRAN_COMMAND=${FORTRAN_TEST_BINARY_PATH}/$(python /src/hpc-predict/hpc-predict-io/test/test_sequential_args.py --test "$@")

echo "HPC_PREDICT_IO_TEST_FORTRAN_COMMAND = ${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}"

cd ${CI_CACHE_FOLDER}

eval ${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}