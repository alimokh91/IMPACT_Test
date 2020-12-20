#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test
export PATH=${FORTRAN_TEST_BINARY_PATH}:${PATH}

HPC_PREDICT_IO_TEST_FORTRAN_COMMAND=$(python /src/hpc-predict/hpc-predict-io/test/test_sequential_args.py --test "$@")

echo "HPC_PREDICT_IO_TEST_FORTRAN_COMMAND = ${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}"

IFS='.' read -r test_module test_class <<<"$1"
[[ -d "${CI_CACHE_FOLDER}/decrypt" ]] && decrypt_dir="decrypt/"
cd ${CI_CACHE_FOLDER}/${test_module}/${decrypt_dir}${test_class}

set -x
echo "$(pwd)"
eval ${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}
set +x
