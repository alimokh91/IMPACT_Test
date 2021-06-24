#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"
COMMIT_SHA=$(git -C ../.. rev-parse HEAD)

if [ -z "${DEBUG_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

export HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG:-"${CONTAINER_REGISTRY}/hpc-predict/io/debug:${COMMIT_SHA}"}
export DEBUG_IMAGE=${DEBUG_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/impact/debug:${COMMIT_SHA}"}
export DEBUG_IMAGE_UNTAGGED=${DEBUG_IMAGE_UNTAGGED:-"${CONTAINER_REGISTRY}/hpc-predict/impact/debug"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    set -x
    if [[ -z $(docker images -q ${HPC_PREDICT_IO_DEBUG}) ]]; then
        echo "Error: Docker base image ${HPC_PREDICT_IO_DEBUG} for HPC-PREDICT-IO not available"
        exit 1
    fi
    if [[ -z $(docker images -q ${DEBUG_IMAGE}) ]]; then
        docker build --build-arg HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG} -f debug/Dockerfile -t ${DEBUG_IMAGE} ..
        if [[ -n ${DEBUG_IMAGE_UNTAGGED} ]]; then
          docker tag ${DEBUG_IMAGE} ${DEBUG_IMAGE_UNTAGGED}
        fi
    fi
    set +x
fi

echo "Run \"export HPC_PREDICT_IMPACT_IMAGE=${DEBUG_IMAGE_UNTAGGED}\" or \"export HPC_PREDICT_IMPACT_IMAGE=${DEBUG_IMAGE}\" to automatically use the built image with docker/run scripts."
cd -
