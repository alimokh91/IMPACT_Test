#!/bin/bash

set -euxo pipefail

export CONTAINER_IMAGE="lukasgd/hpc-predict:impact-test"

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${CONTAINER_IMAGE}) ]]; then
        docker build --rm=false -f ../../docker/Dockerfile-IMPACT -t ${CONTAINER_IMAGE} ../..
    fi
fi

mkdir tmp
function clean_up_test_data {
  rm -r tmp
}
trap clean_up_test_data EXIT

./test_impact_driver.sh

