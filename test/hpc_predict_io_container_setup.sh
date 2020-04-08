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
    MPI_NUM_NODES=(${MPI_NUM_PROCS[0]} $(python -c "print(min(${SLURM_JOB_NUM_NODES}-1, ${MPI_NUM_PROCS[1]}))"))
    for i in $(seq 0 1); do
        MPIEXEC_CMD+=("srun -N ${MPI_NUM_NODES[$i]} -n ${MPI_NUM_PROCS[$i]}")
    done
    CONTAINER_RUN_CMD="sarus run --mpi"
    CONTAINER_ENTRYPOINT=()
    for i in $(seq 0 1); do
        CONTAINER_ENTRYPOINT+=("bash -c")
    done
    MPIEXEC_CMD_SINGLE="srun -N 1 -n 1"
    CONTAINER_ENTRYPOINT_LOCAL="bash -c"
else
    MOUNT_OPT="-v ${MOUNT_HOST_DIR}:${MOUNT_CONTAINER_DIR}"

    MPIEXEC_CMD=()
    for i in $(seq 0 1); do
        MPIEXEC_CMD+=("")
    done
    CONTAINER_RUN_CMD="docker run -u $(id -u ${USER}):$(id -g ${USER})"
    #with strace: CONTAINER_RUN_CMD="docker run" #" -u $(id -u ${USER}):$(id -g ${USER})"
    CONTAINER_ENTRYPOINT=()
    for i in $(seq 0 1); do
        CONTAINER_ENTRYPOINT+=("mpiexec -np ${MPI_NUM_PROCS[$i]} bash -c")
    done
    MPIEXEC_CMD_SINGLE=""
    CONTAINER_ENTRYPOINT_LOCAL="bash -c"
fi

# Set MPI environment variable names specific to MPI version
if echo ${MPI_MASTER_HOST} | grep -i daint  > /dev/null; then
  echo "Using MPICH..." 
  MPI_RANK=SLURM_PROCID
  MPI_SIZE=SLURM_NPROCS
else
  MPI_EXEC_VERSION_STDOUT=$(bash -l -c "mpiexec --version")
  if echo "${MPI_EXEC_VERSION_STDOUT}" | grep -i openmpi > /dev/null; then
    echo "Using OpenMPI..."
    MPI_RANK=PMIX_RANK
    MPI_SIZE=PMIX_SIZE
  elif echo "${MPI_EXEC_VERSION_STDOUT}" | grep -i mpich > /dev/null; then
    echo "Using MPICH..." 
    MPI_RANK=PMI_RANK
    MPI_SIZE=PMI_SIZE
  else
    echo "Failed to identify MPI installation from 'mpiexec --version' - exiting."
    exit 1
  fi
fi
