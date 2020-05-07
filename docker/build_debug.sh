#!/bin/bash

set -euxo pipefail

if [ -z "${DEPLOY_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-lukasgd}
fi

export HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG:-"${CONTAINER_REGISTRY}/hpc-predict:io-debug"}
export DEPLOY_IMAGE=${DEPLOY_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict:impact-debug"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${HPC_PREDICT_IO_DEBUG}) ]]; then
        echo "Error: Docker base image ${HPC_PREDICT_IO_DEBUG} for HPC-PREDICT-IO not available"
        exit 1
    fi
    if [[ -z $(docker images -q ${DEPLOY_IMAGE}) ]]; then
        docker build --build-arg HPC_PREDICT_IO_DEBUG=${HPC_PREDICT_IO_DEBUG} -f debug/Dockerfile -t ${DEPLOY_IMAGE} ..
    fi
fi

