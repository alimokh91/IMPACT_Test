#!/bin/bash

MPI_NUM_PROCS=(1 1)

MOUNT_HOST_DIR="$(pwd)/tmp"
MOUNT_CONTAINER_DIR="/mnt/test"

source hpc_predict_io_container_setup.sh
function container_shutdown {
  source hpc_predict_io_container_shutdown.sh
}
trap container_shutdown EXIT

function run_test() {
  echo "*** $1 started ***"

  # Filenames
  FORTRAN_TEST_BINARY_PATH=/src/hpc-predict/hpc-predict-io/install/bin/test
  fortran_exec_name=$2

  filename_exec=${FORTRAN_TEST_BINARY_PATH}/${fortran_exec_name}
  filename_mri=${fortran_exec_name}.h5
#  filename_out=${fortran_exec_name}.out
#  filename_err=${fortran_exec_name}.err

  IFS='.' read -r python_module test_class <<<"$1"
  python_cmd_get_fortran_args="from ${python_module} import get_fortran_args, ${test_class}; print(get_fortran_args(${test_class}))"
  fortran_args=$(${MPIEXEC_CMD_SINGLE} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT_LOCAL} \
  "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
    HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} \
    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
    python -c \"${python_cmd_get_fortran_args}\"")

  BASH_CMD=()
  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
    HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} \
    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
      python -m unittest -v $1") #/src/hpc-predict/hpc-predict-io/test/test_sequential_mr_io_container.py")
  BASH_CMD+=("cd ${MOUNT_CONTAINER_DIR} && \
     ${filename_exec} ${fortran_args}")
#    source scl_source enable devtoolset-7 && \
#    strace -o ${fortran_exec_name}.strace ${filename_exec} ${fortran_args}") # ${filename_mri} 1> ${filename_out} 2> ${filename_err}")

  #BASH_CMD+=("source ~/src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
  #  HPC_PREDICT_IO_TEST_MRI_PATH=${filename_mri} python test/test_sequential_mr_io_container.py")
  #BASH_CMD+=("${filename_exec} ${filename_mri} 1> ${filename_out} 2> ${filename_err}")

  pid=()

  for i in $(seq 0 1); do
      set -x
      ${MPIEXEC_CMD[$i]} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT[$i]} "${BASH_CMD[$i]}" &
  #    bash -c "${BASH_CMD[$i]}" &
      set +x
      pid+=("$!")
  done

  for i in $(seq 0 1); do
      wait ${pid[$i]}
  done

  echo "*** $1 finished ***"
}

run_test 'test_sequential_mr_io_container.TestSpatialMRI'       'mr_io_test_reader_container'
run_test 'test_sequential_mr_io_container.TestSpaceTimeMRI'     'mr_io_test_reader_space_time_container'
run_test 'test_sequential_mr_io_container.TestFlowMRI'          'mr_io_test_reader_flow_container'
run_test 'test_sequential_mr_io_container.TestSegmentedFlowMRI' 'mr_io_test_reader_segmented_flow_container'

run_test 'test_sequential_mr_io_bidirectional_container.TestSpatialMRIBidirectional'       'mr_io_test_reader_writer'
run_test 'test_sequential_mr_io_bidirectional_container.TestSpaceTimeMRIBidirectional'     'mr_io_test_reader_writer_space_time'
run_test 'test_sequential_mr_io_bidirectional_container.TestFlowMRIBidirectional'          'mr_io_test_reader_writer_flow'
run_test 'test_sequential_mr_io_bidirectional_container.TestSegmentedFlowMRIBidirectional' 'mr_io_test_reader_writer_segmented_flow'

