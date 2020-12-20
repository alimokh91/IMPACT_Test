#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

export FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test

# This must be passed to the container in an environment variable
# export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE}

IFS='.' read -r test_module test_class <<<"$1"
[[ -d "${CI_CACHE_FOLDER}/decrypt" ]] && decrypt_dir="decrypt/" || decrypt_dir=""
cd ${CI_CACHE_FOLDER}/${decrypt_dir}${test_module}/${test_class}

set -x
echo "$(pwd)"
python -m unittest -v "$@"
set +x

