#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"
COMMIT_SHA=$(git -C ../.. rev-parse HEAD)

if [ -z "${DEPLOY_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

# Setting BUILD_ENV_IMAGE as an environment variable for this script allows building on different base image.
# In that case you probably als want to provide a custom DEPLOY_IMAGE name
export BUILD_ENV_IMAGE=${BUILD_ENV_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/build-env"}
export DEPLOY_IMAGE=${DEPLOY_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/io/deploy:${COMMIT_SHA}"}
export DEPLOY_IMAGE_UNTAGGED=${DEPLOY_IMAGE_UNTAGGED:-"${CONTAINER_REGISTRY}/hpc-predict/io/deploy"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    set -x
    if [[ -z $(docker images -q ${BUILD_ENV_IMAGE}) ]]; then
        docker build -f ../docker/build-env/Dockerfile -t ${BUILD_ENV_IMAGE} ../docker/build-env
    fi
    if [[ -z $(docker images -q ${DEPLOY_IMAGE}) ]]; then
        docker build --build-arg BUILD_ENV=${BUILD_ENV_IMAGE} -f ../docker/deploy/Dockerfile -t ${DEPLOY_IMAGE} ../..
        if [[ -n ${DEPLOY_IMAGE_UNTAGGED} ]]; then
          docker tag ${DEPLOY_IMAGE} ${DEPLOY_IMAGE_UNTAGGED}
        fi
    fi
    set +x
fi

echo "Usage ${DEPLOY_IMAGE_UNTAGGED}\" for testing the pipeline with DVC and ${DEPLOY_IMAGE} for committing (reproducible) long-term results."
echo "Run \"export HPC_PREDICT_IO_IMAGE=${DEPLOY_IMAGE_UNTAGGED}\" or \"export HPC_PREDICT_IO_IMAGE=${DEPLOY_IMAGE}\" to automatically use the built image with docker/run scripts."
cd -