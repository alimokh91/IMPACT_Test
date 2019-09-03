#!/bin/bash

set -euxo pipefail

IMPACT_DEBUG_INFO_DIR=/tmp/impact_debug_info
MPI_HOSTS=(localhost)
NUM_MPI_PROCESSES_PER_HOST=4
NUM_CSSH_XTERM_ROWS=2

function cleanup {
  rm -r ${IMPACT_DEBUG_INFO_DIR}
}

trap cleanup EXIT

mkdir -p ${IMPACT_DEBUG_INFO_DIR}
rm ${IMPACT_DEBUG_INFO_DIR}/mpi_processes || true
rm -r ${IMPACT_DEBUG_INFO_DIR}/xterm_processes || true
mkdir -p ${IMPACT_DEBUG_INFO_DIR}/xterm_processes

# Run mpiexec async without waiting for its completion
# Write PIDs and MPI ranks to file that is read again by cssh
mpiexec -np $(( ${#MPI_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c "echo \"Hello from MPI rank \${PMI_RANK} on \$(hostname) with PID \$(echo \$\$)!\"; flock -x -w 5 ${IMPACT_DEBUG_INFO_DIR}/mpi_processes echo \"\${PMI_RANK} \$(echo \$\$)\" >> ${IMPACT_DEBUG_INFO_DIR}/mpi_processes; cd ../prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ./impact.exe" & 

# Not needed for mpi_processes an associative array
#cat ${IMPACT_DEBUG_INFO_DIR}/mpi_processes | sort | tee ${IMPACT_DEBUG_INFO_DIR}/mpi_processes

xterm_hosts=()

for h in "${MPI_HOSTS[@]}"; do
  for i in $(seq 1 ${NUM_MPI_PROCESSES_PER_HOST}); do
    echo "Adding ( ${h} )"
    xterm_hosts+=( "${h}" )
  done
done

echo "MPI_HOSTS=${MPI_HOSTS[@]}"
echo "NUM_MPI_PROCESSES_PER_HOST=${NUM_MPI_PROCESSES_PER_HOST}"
echo "xterm_hosts is ${xterm_hosts[@]}"

sleep 0.5

# Initialize all 
cssh --rows ${NUM_CSSH_XTERM_ROWS} "${xterm_hosts[@]}" -a "cd src/IMPACT/debug && source xterm_attach.sh"


