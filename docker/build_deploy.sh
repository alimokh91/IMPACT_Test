#!/bin/bash

set -euxo pipefail

if [ -z "${DEPLOY_IMAGE-}" ]; then
    CONTAINER_REGISTRY=${CONTAINER_REGISTRY:-lukasgd}
fi

export HPC_PREDICT_IO_DEPLOY=${HPC_PREDICT_IO_DEPLOY:-"${CONTAINER_REGISTRY}/hpc-predict:io-deploy"}
export DEPLOY_IMAGE=${DEPLOY_IMAGE:-"${CONTAINER_REGISTRY}/hpc-predict:impact-deploy"}

MPI_MASTER_HOST=$(hostname)
if [[ -z $(echo ${MPI_MASTER_HOST} | grep -i daint) ]]; then
    if [[ -z $(docker images -q ${HPC_PREDICT_IO_DEPLOY}) ]]; then
        echo "Error: Docker base image ${HPC_PREDICT_IO_DEPLOY} for HPC-PREDICT-IO not available"
        exit 1
    fi
    if [[ -z $(docker images -q ${DEPLOY_IMAGE}) ]]; then
        docker build --build-arg HPC_PREDICT_IO_DEPLOY=${HPC_PREDICT_IO_DEPLOY} -f deploy/Dockerfile -t ${DEPLOY_IMAGE} ..
    fi
fi

