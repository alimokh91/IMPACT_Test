#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"
COMMIT_SHA=$(git -C ../.. rev-parse HEAD)

if [ -z "${DEPLOY_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

export HPC_PREDICT_IO_DEPLOY=${HPC_PREDICT_IO_DEPLOY:-"${CONTAINER_REGISTRY}/hpc-predict/io/deploy:${COMMIT_SHA}"}
export DEPLOY_IMAGE=${DEPLOY_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/impact/deploy:${COMMIT_SHA}"}
export DEPLOY_IMAGE_UNTAGGED=${DEPLOY_IMAGE_UNTAGGED:-"${CONTAINER_REGISTRY}/hpc-predict/impact/deploy"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    set -x
    if [[ -z $(docker images -q ${HPC_PREDICT_IO_DEPLOY}) ]]; then
        echo "Error: Docker base image ${HPC_PREDICT_IO_DEPLOY} for HPC-PREDICT-IO not available"
        exit 1
    fi
    if [[ -z $(docker images -q ${DEPLOY_IMAGE}) ]]; then
        docker build --build-arg HPC_PREDICT_IO_DEPLOY=${HPC_PREDICT_IO_DEPLOY} -f deploy/Dockerfile -t ${DEPLOY_IMAGE} ..
        if [[ -n ${DEPLOY_IMAGE_UNTAGGED} ]]; then
          docker tag ${DEPLOY_IMAGE} ${DEPLOY_IMAGE_UNTAGGED}
        fi
    fi
    set +x
fi

echo "Usage ${DEPLOY_IMAGE_UNTAGGED}\" for testing the pipeline with DVC and ${DEPLOY_IMAGE} for committing (reproducible) long-term results."
echo "Run \"export HPC_PREDICT_IMPACT_IMAGE=${DEPLOY_IMAGE_UNTAGGED}\" or \"export HPC_PREDICT_IMPACT_IMAGE=${DEPLOY_IMAGE}\" to automatically use the built image with docker/run scripts."
cd -
