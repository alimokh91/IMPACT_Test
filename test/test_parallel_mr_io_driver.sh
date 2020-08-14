#!/bin/bash

set -euo pipefail

HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=$((2**1))
MPI_NUM_PROCS=(1 ${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE})

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

  #filename_exec=${FORTRAN_TEST_BINARY_PATH}/${fortran_exec_name}
  # filename_mri=${fortran_exec_name}.h5
  # filename_out=${fortran_exec_name}_\$\{${MPI_RANK}\}.out
  # filename_err=${fortran_exec_name}_\$\{${MPI_RANK}\}.err

#  IFS='.' read -r python_module test_class <<<"$1"
#  python_cmd_get_fortran_args="from ${python_module} import get_fortran_args, ${test_class}; print(get_fortran_args(${test_class}))"
#  fortran_cmd=$(${MPIEXEC_CMD_SINGLE} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT_LOCAL} \
#  "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
#    MPI_RANK_VAR=${MPI_RANK_VAR} MPI_SIZE_VAR=${MPI_SIZE_VAR} \
#    HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE} \
#    python -c \"${python_cmd_get_fortran_args}\"")

  BASH_CMD=()
  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
    export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test && \
    export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE} && \
    python -m unittest -v $1")
#    strace -o ${fortran_exec_name}.strace \
#    python -m unittest -v $1") # with strace

  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
    export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test && \
    export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE} && \
    export HPC_PREDICT_IO_TEST_FORTRAN_COMMAND=${FORTRAN_TEST_BINARY_PATH}/\$(python /src/hpc-predict/hpc-predict-io/test/test_parallel_args.py --test $1) && \
    echo \"HPC_PREDICT_IO_TEST_FORTRAN_COMMAND = \${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}\" && \
    cd ${MOUNT_CONTAINER_DIR} && \
    eval \${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}")
#  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#    export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE=${HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE} && \
#    cd /src/hpc-predict/hpc-predict-io/test/ && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_EXEC=\$(python  test_parallel_args.py --test $1 --type exec) && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_ARGS=(\$(python test_parallel_args.py --test $1 --type args)) && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_OUT=\$( python  test_parallel_args.py --test $1 --type out) && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_ERR=\$( python  test_parallel_args.py --test $1 --type err) && \
#    cd ${MOUNT_CONTAINER_DIR} && \
#    set -x && \
#    ${FORTRAN_TEST_BINARY_PATH}/\${HPC_PREDICT_IO_TEST_FORTRAN_EXEC} \${HPC_PREDICT_IO_TEST_FORTRAN_ARGS[@]} 1> \${HPC_PREDICT_IO_TEST_FORTRAN_OUT} 2> \${HPC_PREDICT_IO_TEST_FORTRAN_ERR}")

#  BASH_CMD+=("cd ${MOUNT_CONTAINER_DIR} && \
#    ${FORTRAN_TEST_BINARY_PATH}/${fortran_cmd}")
##  BASH_CMD+=("cd ${MOUNT_CONTAINER_DIR} && \
##    strace -o ${fortran_exec_name}_\${${MPI_RANK}}.strace \
##    ${filename_exec} ${fortran_args}") # with strace

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

run_test 'test_parallel_mr_io_container.TestSpatialMRI'             'mr_io_test_parallel_reader_container'
run_test 'test_parallel_mr_io_container.TestSpaceTimeMRI'           'mr_io_test_parallel_reader_space_time_container'
run_test 'test_parallel_mr_io_container.TestFlowMRI'                'mr_io_test_parallel_reader_flow_container'
run_test 'test_parallel_mr_io_container.TestSegmentedFlowMRI'       'mr_io_test_parallel_reader_segmented_flow_container'
run_test 'test_parallel_mr_io_container.TestFlowMRIPadded'          'mr_io_test_parallel_reader_flow_padded_container'
run_test 'test_parallel_mr_io_container.TestSegmentedFlowMRIPadded' 'mr_io_test_parallel_reader_segmented_flow_padded_container'

run_test 'test_parallel_mr_io_bidirectional_container.TestSpatialMRIBidirectional'             'mr_io_test_parallel_reader_writer'
run_test 'test_parallel_mr_io_bidirectional_container.TestSpaceTimeMRIBidirectional'           'mr_io_test_parallel_reader_writer_space_time'
run_test 'test_parallel_mr_io_bidirectional_container.TestFlowMRIBidirectional'                'mr_io_test_parallel_reader_writer_flow'
run_test 'test_parallel_mr_io_bidirectional_container.TestSegmentedFlowMRIBidirectional'       'mr_io_test_parallel_reader_writer_segmented_flow'
run_test 'test_parallel_mr_io_bidirectional_container.TestFlowMRIPaddedBidirectional'                      'mr_io_test_parallel_reader_writer_flow_padded'
run_test 'test_parallel_mr_io_bidirectional_container.TestSegmentedFlowMRIPaddedBidirectional'             'mr_io_test_parallel_reader_writer_segmented_flow_padded'
run_test 'test_parallel_mr_io_bidirectional_container.TestFlowMRIPaddedToSpaceTimeBidirectional'           'mr_io_test_parallel_reader_writer_flow_padded_to_space_time'
run_test 'test_parallel_mr_io_bidirectional_container.TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional'  'mr_io_test_parallel_reader_writer_segmented_flow_padded_to_space_time'
