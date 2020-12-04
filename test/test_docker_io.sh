#!/bin/bash

set -euxo pipefail

CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}

export BUILD_ENV_IMAGE="${CONTAINER_REGISTRY}/hpc-predict/io/build-env"
export DEPLOY_IMAGE="${CONTAINER_REGISTRY}/hpc-predict/io/deploy"

../docker/build_deploy.sh

export CONTAINER_IMAGE=${DEPLOY_IMAGE}

# MPI_MASTER_HOST=$(hostname)
# if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
#     if [[ -z $(docker images -q ${BUILD_ENV_IMAGE}) ]]; then
#         docker build -f ../docker/build-env/Dockerfile -t ${BUILD_ENV_IMAGE} ..
#     fi
#     if [[ -z $(docker images -q ${CONTAINER_IMAGE}) ]]; then
#         docker build -f ../docker/deploy/Dockerfile -t ${CONTAINER_IMAGE} ..
#     fi
# fi

mkdir -p tmp
#function clean_up_test_data {
#  rm -r tmp
#}
#trap clean_up_test_data EXIT

./test_sequential_mr_io_driver.sh
./test_parallel_mr_io_driver.sh

