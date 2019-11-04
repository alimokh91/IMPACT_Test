#!/bin/bash

echo "Checking environment variables..."
#if [[ -z "${HPC_PREDICT_IO_DIR}" ]]; then echo "Error: Environment variable HPC_PREDICT_IO_DIR is not set"; exit 1; fi
if [[ -z "${HDF5_DIR}" ]]; then echo "Error: Environment variable HDF5_DIR is not set"; exit 1; fi

# Use proper username as well for SSH login
MPI_MASTER_HOST=(localhost) # (gpucandoit.artorg.unibe.ch) or (daint)
MPI_MASTER_IMPACT_DIR="$(pwd)/.."

if [[ -z "${MPI_MASTER_IMPACT_DIR}" ]]; then echo "Error: Environment variable MPI_MASTER_IMPACT_DIR is not set"; exit 1; fi

NUM_MPI_PROCESSES_PER_HOST=4
NUM_CSSH_XTERM_ROWS=2

SSH_MASTER_PROXY=""
SSH_WORKER_PROXY=""

if echo ${MPI_MASTER_HOST} | grep -i gpucandoit > /dev/null; then
  # for Dario/Derick - load correct modules 
  MPI_WORKER_HOSTS=(gpucandoit.artorg.unibe.ch)
  MPI_MASTER_IMPACT_DIR="/media/Fast_and_Furious/derick/IMPACT/" 
  MPI_EXEC_COMMAND="module load mpich-3.2 hdf5-1.10.1; mpiexec"
  MPI_EXEC_NUM_PROCESSES="-np $(( ${#MPI_WORKER_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} ))"
  #"ssh derick@gpucandoit.artorg.unibe.ch bash -l -c "\"module load mpich-3.2 hdf5-1.10.1; mpiexec\"""
elif echo ${MPI_MASTER_HOST} | grep -i daint > /dev/null; then
  MPI_WORKER_HOSTS=(nid04276 nid04277 nid04278 nid04279)
  
  MPI_MASTER_IMPACT_DIR="/scratch/snx3000/lukasd/IMPACT"

  MPI_EXEC_COMMAND="source shell_env_alloc.sh; source /apps/daint/UES/anfink/cpu/environment; srun "
  MPI_EXEC_NUM_PROCESSES="-N$(( ${#MPI_WORKER_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} ))"
  SSH_WORKER_PROXY="" #ProxyCommand=ssh -q -Y daint.cscs.ch -W %h:%p"
  # TODO: Debug job allocation, writing shell_env_alloc.sh on MPI_MASTER_HOST and MPI_WORKER_HOSTS initialization
  # Currently on Piz Daint perform the following steps:
  # 1. Put the following line in your ~/.ssh/config: 
  #    Host nid*
  #    User <your-username>
  #    ControlMaster no
  #    StrictHostKeyChecking no
  #    ProxyCommand ssh -q -W "%h:%p" daint
  # 2. salloc  --nodes=<desired-num-nodes> --partition=debug --time=00:30:00 -C gpu (assuming that you already have a valid config.txt/MRIs)
  # 3.. cd ~/; declare -x > shell_env_alloc.sh; cd -;
  # 4.. Display allocated nodes: scontrol show job <jobid-allocated-above>
  # 5. Update MPI_WORKER_HOSTS above accordingly
else
  MPI_WORKER_HOSTS=(localhost)
  MPI_EXEC_COMMAND="mpiexec"
  MPI_EXEC_NUM_PROCESSES="-np $(( ${#MPI_WORKER_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} ))"
fi

IMPACT_DEBUG_INFO_DIR="/tmp/impact_debug_info"

# Set MPI environment variable names specific to MPI version
if echo ${MPI_MASTER_HOST} | grep -i daint  > /dev/null; then 
  echo "Using MPICH..." 
  MPI_RANK=SLURM_PROCID
  MPI_SIZE=SLURM_NPROCS
else  
  MPI_EXEC_VERSION_STDOUT=$(ssh -o "${SSH_MASTER_PROXY}" ${MPI_MASTER_HOST} bash -l -c \"${MPI_EXEC_COMMAND} --version \")
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

# Set CSSH executable correctly
case $(uname -s) in
  Linux*)  CSSH_EXEC=cssh;;
  Darwin*) CSSH_EXEC=csshX;;
  *) echo "Unknown operating system..." && exit 1
esac

set -euxo pipefail


echo "Checking if all MPI hosts are reachable via SSH..."
ssh -o "${SSH_MASTER_PROXY}" ${MPI_MASTER_HOST} exit
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh -o "${SSH_WORKER_PROXY}" ${h} exit
done

# TODO: run the kill commands through ssh on MPI_WORKER_HOSTS
function cleanup {
  echo "Cleaning up - killing previously started MPI processes (incl. master process)"
  for h in "${MPI_WORKER_HOSTS[@]}"; do
    mpi_worker_pids=()
    while IFS=' ' read -r key value; do
       mpi_worker_pids+=(${value})
    done < <(ssh -o "${SSH_WORKER_PROXY}" ${h} "cat ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes/*")
    ssh -o "${SSH_WORKER_PROXY}" ${h} "kill -9 ${mpi_worker_pids[@]}"
  done
  kill -9 ${IMPACT_MPI_MASTER_PID} 
  ssh -o "${SSH_MASTER_PROXY}" ${MPI_MASTER_HOST} "rm -r ${IMPACT_DEBUG_INFO_DIR}"
}

trap cleanup EXIT

ssh -o "${SSH_MASTER_PROXY}" ${MPI_MASTER_HOST} "mkdir -p ${IMPACT_DEBUG_INFO_DIR}"
for h in "${MPI_WORKER_HOSTS[@]}"; do
  ssh -o "${SSH_WORKER_PROXY}" ${h} \
           "rm -r ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes || true;\
            rm -r ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes || true;\
            mkdir -p ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/xterm_processes;\
            mkdir -p ${IMPACT_DEBUG_INFO_DIR}/\$(hostname)/mpi_processes"
done

echo "Starting IMPACT debugger version over MPI..."
# Run mpiexec async without waiting for its completion
# Write PIDs and MPI ranks to file that is read again by cssh
#mpiexec -np $(( ${#MPI_HOSTS[@]}*${NUM_MPI_PROCESSES_PER_HOST} )) bash -c "echo \"Hello from MPI rank \${${MPI_RANK}} on \$(hostname) with PID \$(echo \$\$)!\"; flock -x -w 5 ${IMPACT_DEBUG_INFO_DIR}/mpi_processes echo \"\${${MPI_RANK}} \$(echo \$\$)\" >> ${IMPACT_DEBUG_INFO_DIR}/mpi_processes; cd ../prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ./impact_debug.exe" & 

SH_ESC="\\\\\\"
sh_command="echo ${SH_ESC}\"Hello from MPI rank ${SH_ESC}\${${MPI_RANK}} on ${SH_ESC}\$(hostname) with PID ${SH_ESC}\$(echo ${SH_ESC}\$${SH_ESC}\$)!${SH_ESC}\"; echo ${SH_ESC}\"${SH_ESC}\${${MPI_RANK}} ${SH_ESC}\$(echo ${SH_ESC}\$${SH_ESC}\$)${SH_ESC}\" >> ${IMPACT_DEBUG_INFO_DIR}/${SH_ESC}\$(hostname)/mpi_processes/${SH_ESC}\${${MPI_RANK}}; cd ${MPI_MASTER_IMPACT_DIR}/prog; export LD_LIBRARY_PATH=${HDF5_DIR}/lib; exec ./impact_debug.exe"
echo "sh_command=\"${sh_command}\""
ssh -o "${SSH_MASTER_PROXY}" ${MPI_MASTER_HOST} bash -l -c \"${MPI_EXEC_COMMAND} ${MPI_EXEC_NUM_PROCESSES} bash -c \\\"${sh_command}\\\"\"  &
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

echo "Copying xterm_attach.sh and gdb_init to ${MPI_MASTER_HOST}:${IMPACT_DEBUG_INFO_DIR}"

for h in "${MPI_WORKER_HOSTS[@]}"; do
  scp -o "${SSH_WORKER_PROXY}" xterm_attach.sh gdb_init "${h}:${IMPACT_DEBUG_INFO_DIR}/" 
done

echo "Launching CSSH..."

# Initialize all 
${CSSH_EXEC} --config-file clusterssh_config --rows ${NUM_CSSH_XTERM_ROWS} -o "${SSH_WORKER_PROXY}" "${xterm_hosts[@]}" #-a "cd src/IMPACT/debug && source xterm_attach.sh && exec bash"

