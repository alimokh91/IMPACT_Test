#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Container image must be supplied as first parameter."
else
    container_image="$1"
fi

# Make sure that commit inside and outside container are the same
local_commit=$(git rev-parse HEAD)
container_commit=$(docker run --rm "${container_image}" git -C /src/hpc-predict rev-parse HEAD)
git_diff_local_container=$(git diff ${local_commit}..${container_commit} --name-status -- . ':!data/*')
if [ -z "{git_diff_local_container}" ]; then
  echo "Validating code consistency: Git checkout on local filesystem and in container differs only in data directory."
else
  echo "Code consistency error: Git checkout on local filesystem (${local_commit}) and in container ${container_image} (${container_commit}) differ by more than data-folder (must use consistent to commit data with DVC):"
  echo "${git_diff_local_container}" 
  #exit 1 
fi
