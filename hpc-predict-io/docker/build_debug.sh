#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"
COMMIT_SHA=$(git -C ../.. rev-parse HEAD)

if [ -z "${DEBUG_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

# Setting BUILD_ENV_IMAGE as an environment variable for this script allows building on different base image
export BUILD_ENV_IMAGE=${BUILD_ENV_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/build-env:${COMMIT_SHA}"}
export DEBUG_IMAGE=${DEBUG_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/debug:${COMMIT_SHA}"}
export DEBUG_IMAGE_UNTAGGED=${DEBUG_IMAGE_UNTAGGED:-"${CONTAINER_REGISTRY}/hpc-predict/io/debug"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    set -x
    if [[ -z $(docker images -q ${BUILD_ENV_IMAGE}) ]]; then
        docker build -f ../docker/build-env/Dockerfile -t ${BUILD_ENV_IMAGE} ../docker/build-env
    fi
    if [[ -z $(docker images -q ${DEBUG_IMAGE}) ]]; then
        docker build --build-arg BUILD_ENV=${BUILD_ENV_IMAGE} -f ../docker/debug/Dockerfile -t ${DEBUG_IMAGE} ../..
        if [[ -n ${DEBUG_IMAGE_UNTAGGED} ]]; then
          docker tag ${DEBUG_IMAGE} ${DEBUG_IMAGE_UNTAGGED}
        fi
    fi
    set +x
fi

echo "Run \"export HPC_PREDICT_IO_IMAGE=${DEBUG_IMAGE_UNTAGGED}\" or \"export HPC_PREDICT_IO_IMAGE=${DEBUG_IMAGE}\" to automatically use the built image with docker/run scripts."
cd -
