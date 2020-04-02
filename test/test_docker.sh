#!/bin/bash

set -euxo pipefail

export CONTAINER_IMAGE="lukasgd/hpc-predict:io-test"
#export CONTAINER_IMAGE="lukasgd/hpc-predict-io-test-mpi-hdf5-debug"
#export CONTAINER_IMAGE="lukasgd/hpc-predict-io-test-debug"

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${CONTAINER_IMAGE}) ]]; then
        docker build --rm=false -f ../docker/Dockerfile-test -t ${CONTAINER_IMAGE} ..
    #    docker build --rm=false -f ../docker/Dockerfile-debug -t ${CONTAINER_IMAGE} ..
    fi
fi

mkdir tmp
function clean_up_test_data {
  rm -r tmp
}
trap clean_up_test_data EXIT

./test_sequential_mr_io_driver.sh
./test_parallel_mr_io_driver.sh

