#!/bin/bash

set -euxo pipefail

if [ -z "${DEPLOY_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

export BUILD_ENV_IMAGE=${BUILD_ENV_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/build-env"}
export DEPLOY_IMAGE=${DEPLOY_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/deploy"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${BUILD_ENV_IMAGE}) ]]; then
        docker build -f ../docker/build-env/Dockerfile -t ${BUILD_ENV_IMAGE} ../docker/build-env
    fi
    if [[ -z $(docker images -q ${DEPLOY_IMAGE}) ]]; then
        docker build --build-arg BUILD_ENV=${BUILD_ENV_IMAGE} -f ../docker/deploy/Dockerfile -t ${DEPLOY_IMAGE} ../..
    fi
fi

echo "Run \"export HPC_PREDICT_IO_IMAGE=${DEPLOY_IMAGE}\" to automatically use the built image with docker/run scripts."
