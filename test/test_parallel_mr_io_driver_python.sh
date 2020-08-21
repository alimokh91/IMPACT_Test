#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

export FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test

# This must be passed to the container in an environment variable
# export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE}

python -m unittest -v "$@"
