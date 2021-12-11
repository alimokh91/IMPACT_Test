#!/bin/bash

echo "Checking environment variables..."
if [[ -z "${HPC_PREDICT_IO_DIR}" ]]; then echo "Error: Environment variable HPC_PREDICT_IO_DIR is not set"; exit 1; fi
if [[ -z "${HDF5_DIR}" ]]; then echo "Error: Environment variable HDF5_DIR is not set"; exit 1; fi

PROJECT_EXECUTABLE=fortran/test/mr_io_test_impact_mri
MPI_PWD=$(pwd)/..

if [[ -z "${PROJECT_EXECUTABLE}" ]]; then echo "Error: Environment variable PROJECT_EXECUTABLE is not set"; exit 1; fi
if [[ -z "${MPI_PWD}" ]]; then echo "Error: Environment variable MPI_PWD is not set"; exit 1; fi


# Use proper username as well for SSH login
MPI_MASTER_HOST=(localhost)
MPI_WORKER_HOSTS=(localhost)
#MPI_MASTER_PROJECT_SOURCE_DIR="$(pwd)/.."

# Set MPI environment variable names specific to MPI version
if ssh ${MPI_MASTER_HOST} mpiexec --version | grep -i openmpi > /dev/null; then
  echo "Using OpenMPI..."
  MPI_RANK=PMIX_RANK
  MPI_SIZE=PMIX_SIZE
elif ssh ${MPI_MASTER_HOST} mpiexec --version | grep -i mpich > /dev/null; then
  echo "Using MPICH..." 
  MPI_RANK=PMI_RANK
  MPI_SIZE=PMI_SIZE
else 
  echo "Failed to identify MPI installation from 'mpiexec --version' - exiting."
  exit 1
fi

# Set CSSH executable correctly
case $(uname -s) in
  Linux*)  CSSH_EXEC=cssh;;
  Darwin*) CSSH_EXEC=csshX;;
  *) echo "Unknown operating system..." && exit 1
esac

set -euxo pipefail


echo "Checking if all MPI hosts are reachable via SSH..."
ssh ${MPI_MASTER_HOST} exit
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh ${h} exit
done

PROJECT_DEBUG_INFO_DIR=/tmp/hpc_predict_io_debug_info
NUM_MPI_PROCESSES_PER_HOST=4
NUM_CSSH_XTERM_ROWS=2

# TODO: run the kill commands through ssh on MPI_WORKER_HOSTS
function cleanup {
  echo "Cleaning up - killing previously started MPI processes (incl. master process)"
  for h in "${MPI_WORKER_HOSTS[@]}"; do
    mpi_worker_pids=()
    while IFS=' ' read -r key value; do
       mpi_worker_pids+=(${value})
    done < <(ssh ${h} "cat ${PROJECT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes/*")
    ssh ${h} "kill -9 ${mpi_worker_pids[@]}"
  done
  kill -9 ${MPI_MASTER_PID} 
  ssh ${MPI_MASTER_HOST} "rm -r ${PROJECT_DEBUG_INFO_DIR}"
}

trap cleanup EXIT

ssh ${MPI_MASTER_HOST} "mkdir -p ${PROJECT_DEBUG_INFO_DIR}"
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh ${h} "rm -r ${PROJECT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes || true;\
            rm -r ${PROJECT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes || true;\
            mkdir -p ${PROJECT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes;\
            mkdir -p ${PROJECT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes"
done

echo "Starting ${PROJECT_EXECUTABLE} over MPI..."
# Run mpiexec async without waiting for its completion
# Write PIDs and MPI ranks to file that is read again by cssh
#mpiexec -np $(( ${#MPI_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c "echo \"Hello from MPI rank \${${MPI_RANK}} on \$(hostname) with PID \$(echo \$\$)!\"; flock -x -w 5 ${PROJECT_DEBUG_INFO_DIR}/mpi_processes echo \"\${${MPI_RANK}} \$(echo \$\$)\" >> ${PROJECT_DEBUG_INFO_DIR}/mpi_processes; cd ../prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ${PROJECT_EXECUTABLE}" & 

sh_command="echo \\\"Hello from MPI rank \\\${${MPI_RANK}} on \\\$(hostname) with PID \\\$(echo \\\$\\\$)!\\\"; echo \\\"\\\${${MPI_RANK}} \\\$(echo \\\$\\\$)\\\" >> ${PROJECT_DEBUG_INFO_DIR}/\\\$(hostname)/mpi_processes/\\\${${MPI_RANK}}; cd ${MPI_PWD}; LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ${PROJECT_EXECUTABLE}"
echo "sh_command=\"${sh_command}\""
ssh ${MPI_MASTER_HOST} mpiexec -np $(( ${#MPI_WORKER_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c \"${sh_command}\"  &
MPI_MASTER_PID=$!

# Not needed for mpi_processes an associative array
#cat ${PROJECT_DEBUG_INFO_DIR}/mpi_processes | sort | tee ${PROJECT_DEBUG_INFO_DIR}/mpi_processes

xterm_hosts=()

for h in "${MPI_WORKER_HOSTS[@]}"; do
  for i in $(seq 1 ${NUM_MPI_PROCESSES_PER_HOST}); do
    echo "Adding ( ${h} )"
    xterm_hosts+=( "${h}" )
  done
done

echo "MPI_MASTER_HOST=${MPI_MASTER_HOST}"
echo "MPI_WORKER_HOSTS=${MPI_WORKER_HOSTS[@]}"
echo "NUM_MPI_PROCESSES_PER_HOST=${NUM_MPI_PROCESSES_PER_HOST}"
echo "xterm_hosts is ${xterm_hosts[@]}"

sleep 0.5

echo "Launching CSSH..."

# Initialize all 
${CSSH_EXEC} --config-file clusterssh_config --rows ${NUM_CSSH_XTERM_ROWS} "${xterm_hosts[@]}" #-a "cd src/IMPACT/debug && source xterm_attach.sh && exec bash"
