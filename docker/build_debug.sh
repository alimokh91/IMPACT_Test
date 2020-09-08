#!/bin/bash

set -euxo pipefail

if [ -z "${DEBUG_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-local}
fi

export HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG:-"${CONTAINER_REGISTRY}/hpc-predict/io/debug"}
export DEBUG_IMAGE=${DEBUG_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict/impact/debug"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${HPC_PREDICT_IO_DEBUG}) ]]; then
        echo "Error: Docker base image ${HPC_PREDICT_IO_DEBUG} for HPC-PREDICT-IO not available"
        exit 1
    fi
    if [[ -z $(docker images -q ${DEBUG_IMAGE}) ]]; then
        docker build --build-arg HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG} -f debug/Dockerfile -t ${DEBUG_IMAGE} ..
    fi
fi

echo "Run \"export HPC_PREDICT_IO_IMAGE=${DEBUG_IMAGE}\" to automatically use the built image with docker/run scripts."
