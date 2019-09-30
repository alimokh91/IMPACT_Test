#!/bin/bash

echo "Checking environment variables..."
if [[ -z "${HPC_PREDICT_IO_DIR}" ]]; then echo "Error: Environment variable HPC_PREDICT_IO_DIR is not set"; exit 1; fi
if [[ -z "${HDF5_DIR}" ]]; then echo "Error: Environment variable HDF5_DIR is not set"; exit 1; fi
if [[ -z "${MPI_MASTER_IMPACT_DIR}" ]]; then echo "Error: Environment variable MPI_MASTER_IMPACT_DIR is not set"; exit 1; fi


# Set MPI environment variable names specific to MPI version
if mpiexec --version | grep -i openmpi > /dev/null; then
  echo "Using OpenMPI..."
  MPI_RANK=PMIX_RANK
  MPI_SIZE=PMIX_SIZE
elif mpiexec --version | grep -i mpich > /dev/null; then
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

# Use proper username as well for SSH login
MPI_MASTER_HOST=(localhost)
MPI_WORKER_HOSTS=(localhost)
#MPI_MASTER_IMPACT_DIR="$(pwd)/.."

echo "Checking if all MPI hosts are reachable via SSH..."
ssh ${MPI_MASTER_HOST} exit
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh ${h} exit
done

IMPACT_DEBUG_INFO_DIR=/tmp/impact_debug_info
NUM_MPI_PROCESSES_PER_HOST=4
NUM_CSSH_XTERM_ROWS=2

# TODO: run the kill commands through ssh on MPI_WORKER_HOSTS
function cleanup {
  echo "Cleaning up - killing previously started MPI processes (incl. master process)"
  for h in "${MPI_WORKER_HOSTS[@]}"; do
    mpi_worker_pids=()
    while IFS=' ' read -r key value; do
       mpi_worker_pids+=(${value})
    done < <(ssh ${h} "cat ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes/*")
    ssh ${h} "kill -9 ${mpi_worker_pids[@]}"
  done
  kill -9 ${IMPACT_MPI_MASTER_PID} 
  ssh ${MPI_MASTER_HOST} "rm -r ${IMPACT_DEBUG_INFO_DIR}"
}

trap cleanup EXIT

ssh ${MPI_MASTER_HOST} "mkdir -p ${IMPACT_DEBUG_INFO_DIR}"
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh ${h} "rm -r ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes || true;\
            rm -r ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes || true;\
            mkdir -p ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes;\
            mkdir -p ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes"
done

echo "Starting IMPACT debugger version over MPI..."
# Run mpiexec async without waiting for its completion
# Write PIDs and MPI ranks to file that is read again by cssh
#mpiexec -np $(( ${#MPI_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c "echo \"Hello from MPI rank \${${MPI_RANK}} on \$(hostname) with PID \$(echo \$\$)!\"; flock -x -w 5 ${IMPACT_DEBUG_INFO_DIR}/mpi_processes echo \"\${${MPI_RANK}} \$(echo \$\$)\" >> ${IMPACT_DEBUG_INFO_DIR}/mpi_processes; cd ../prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ./impact_debug.exe" & 

sh_command="echo \\\"Hello from MPI rank \\\${${MPI_RANK}} on \\\$(hostname) with PID \\\$(echo \\\$\\\$)!\\\"; echo \\\"\\\${${MPI_RANK}} \\\$(echo \\\$\\\$)\\\" >> ${IMPACT_DEBUG_INFO_DIR}/\\\$(hostname)/mpi_processes/\\\${${MPI_RANK}}; cd ${MPI_MASTER_IMPACT_DIR}/prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ./impact_debug.exe"
echo "sh_command=\"${sh_command}\""
ssh ${MPI_MASTER_HOST} mpiexec -np $(( ${#MPI_WORKER_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c \"${sh_command}\"  &
IMPACT_MPI_MASTER_PID=$!

# Not needed for mpi_processes an associative array
#cat ${IMPACT_DEBUG_INFO_DIR}/mpi_processes | sort | tee ${IMPACT_DEBUG_INFO_DIR}/mpi_processes

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

