#!/bin/bash

set -euo pipefail

HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=$((2**3)) # 7
MPI_NUM_PROCS=(1 ${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE})

source hpc_predict_io_container_setup.sh
function container_shutdown {
  source hpc_predict_io_container_shutdown.sh
}
trap container_shutdown EXIT

function run_test() {
  echo "*** $1 started ***"
  IFS='.' read -r test_module test_class <<<"$1"
  mkdir -p ${CI_CACHE_FOLDER}/${test_module}/${test_class}


#  # Filenames
#  FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test
#  fortran_exec_name=$2
#
#  filename_exec=${FORTRAN_TEST_BINARY_PATH}/${fortran_exec_name}
#  IFS='.' read -r python_module test_class <<<"$1"
#  filename_mri=${test_class}.h5
#  # filename_out=${fortran_exec_name}_\$\{${MPI_RANK}\}.out
#  # filename_err=${fortran_exec_name}_\$\{${MPI_RANK}\}.err

  BASH_CMD=()
#  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
#    HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} \
#    HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS=${HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS} \
#    IMPACT_CONFIG_TEMPLATE_PATH=/src/hpc-predict/hpc-predict-io/python/config.txt.j2 \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
#    python -m unittest -v $1") #/src/hpc-predict/hpc-predict-io/test/test_parallel_mr_io_container.py")
#  #BASH_CMD+=("source ~/src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#  #  HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} python test/test_parallel_mr_io_container.py")

  CONTAINER_ENV="--env HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE}"
  BASH_CMD+=("/src/hpc-predict/hpc-predict-io/test/test_impact_driver_python.sh $1")

  pid=()

  i=0
  set -x
  ${MPIEXEC_CMD[$i]} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_ENV} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT[$i]} "${BASH_CMD[$i]}" &
#    bash -c "${BASH_CMD[$i]}" &
  set +x
  unset i
  pid+=("$!")

#  function run_python_cmd() {
#    set -x
#    echo $(${MPIEXEC_CMD_SINGLE} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT_LOCAL} \
#  "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
#    HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} \
#    HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS=${HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS} \
#    IMPACT_CONFIG_TEMPLATE_PATH=/src/hpc-predict/hpc-predict-io/python/config.txt.j2 \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
#    MPI_RANK=${MPI_RANK} MPI_SIZE=${MPI_SIZE} \
#    python -c \"$1\"")
#    set +x
#  }
#
#  IFS='.' read -r python_module test_class <<<"$1"
#  python_cmd_get_config_args="from test_impact_common import get_config_args; from ${python_module} import ${test_class}; print(get_config_args(${test_class}))"
#  config_args=$(run_python_cmd "${python_cmd_get_config_args}")
#  python_cmd_get_impact_args="from test_impact_common import get_impact_args; from ${python_module} import ${test_class}; print(get_impact_args(${test_class}))"
#  impact_args=$(run_python_cmd "${python_cmd_get_impact_args}")

  # Config generator
    set -x
#  ${MPIEXEC_CMD_SINGLE} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_ENV} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT_LOCAL} \
#  "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python \
#    python ${config_args}"
  CONFIG_CMD="/src/hpc-predict/hpc-predict-io/test/test_impact_driver_config.sh $1"
  ${MPIEXEC_CMD[0]} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_ENV} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT[0]} "${CONFIG_CMD}" \

    set +x

  # IMPACT
#  BASH_CMD+=("cd ${MOUNT_CONTAINER_DIR} && \
#    ${filename_exec} ${impact_args}")
##    source scl_source enable devtoolset-7 && \
##    strace -o ${fortran_exec_name}_\${${MPI_RANK}}.strace ${filename_exec} ${fortran_args}") # 1> ${filename_out} 2> ${filename_err}")
#
#  #BASH_CMD+=("source ~/src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#  #  HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} python test/test_parallel_mr_io_container.py")
#  #BASH_CMD+=("${filename_exec} ${filename_mri} 1> ${filename_out} 2> ${filename_err}")


  BASH_CMD+=("/src/hpc-predict/hpc-predict-io/test/test_impact_driver_impact.sh $1")
  i=1
  set -x
  ${MPIEXEC_CMD[$i]} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_ENV} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT[$i]} "${BASH_CMD[$i]}" &
#    bash -c "${BASH_CMD[$i]}" &
  set +x
  unset i
  pid+=("$!")

  for i in $(seq 0 1); do
      wait ${pid[$i]}
  done

  rmdir ${CI_CACHE_FOLDER}/${test_module}/${test_class}
  rmdir ${CI_CACHE_FOLDER}/${test_module}
  echo "*** $1 finished ***"
}

run_test 'test_impact_integration.TestImpactInput'           'mr_io_test_impact_input'
run_test 'test_impact_integration.TestImpactInputPadding'    'mr_io_test_impact_input'
run_test 'test_impact_integration.TestImpactMRI'             'mr_io_test_impact_mri'
