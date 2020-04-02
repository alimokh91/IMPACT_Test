#!/bin/bash

set -euo pipefail

#CONTAINER_IMAGE="lukasgd/hpc-predict-io-test"

### FIXME: rewrite this mess with functions returning the correct commands (not variables)

MPI_MASTER_HOST=$(hostname)
if echo ${MPI_MASTER_HOST} | grep -i gpucandoit > /dev/null; then
    echo "Version for gpucandoit currently not implemented - exiting."
    exit -1
elif echo ${MPI_MASTER_HOST} | grep -i daint > /dev/null; then

    MODULE_SARUS_LOADED=false
    if [ $( which sarus 2> /dev/null ) ]; then
        echo "Found module sarus..."
    else
        echo "Sarus not found. Loading module sarus..."
        module load sarus
        MODULE_SARUS_LOADED=true
    fi
    CONTAINER_IMAGE="load/${CONTAINER_IMAGE}"

    MOUNT_OPT="--mount=type=bind,source=${MOUNT_HOST_DIR},destination=${MOUNT_CONTAINER_DIR}"

    MPIEXEC_CMD=()
    MPIEXEC_CMD_COMPILE=()
    for i in $(seq 0 1); do
        MPIEXEC_CMD+=("srun -N ${MPI_NUM_PROCS[$i]} -n ${MPI_NUM_PROCS[$i]}")
        MPIEXEC_CMD_COMPILE+=("srun -N 1 -n 1")
    done
    CONTAINER_RUN_CMD="sarus run --mpi"
    CONTAINER_ENTRYPOINT=()
    for i in $(seq 0 1); do
        CONTAINER_ENTRYPOINT+=("bash -c")
    done
else
    MOUNT_OPT="-v ${MOUNT_HOST_DIR}:${MOUNT_CONTAINER_DIR}"

    MPIEXEC_CMD=()
    MPIEXEC_CMD_COMPILE=()
    for i in $(seq 0 1); do
        MPIEXEC_CMD+=("")
        MPIEXEC_CMD_COMPILE+=("")
    done
    CONTAINER_RUN_CMD="docker run -u $(id -u ${USER}):$(id -g ${USER})"
    #with strace: CONTAINER_RUN_CMD="docker run" #" -u $(id -u ${USER}):$(id -g ${USER})"
    CONTAINER_ENTRYPOINT=()
    for i in $(seq 0 1); do
        CONTAINER_ENTRYPOINT+=("mpiexec -np ${MPI_NUM_PROCS[$i]} bash -c")
    done
    CONTAINER_ENTRYPOINT_LOCAL="bash -c"
fi
