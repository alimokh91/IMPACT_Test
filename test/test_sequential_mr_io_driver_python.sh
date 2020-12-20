#!/bin/bash

set -euo pipefail

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test

export FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test

IFS='.' read -r test_module test_class <<<"$1"
[[ -d "${CI_CACHE_FOLDER}/decrypt" ]] && decrypt_dir="decrypt/"
cd ${CI_CACHE_FOLDER}/${test_module}/${decrypt_dir}${test_class}

set -x
echo "$(pwd)"
python -m unittest -v "$@"
set +x


