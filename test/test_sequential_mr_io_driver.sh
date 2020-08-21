#!/bin/bash

MPI_NUM_PROCS=(1 1)

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

#  filename_exec=${FORTRAN_TEST_BINARY_PATH}/${fortran_exec_name}
#  filename_mri=${fortran_exec_name}.h5
#  filename_out=${fortran_exec_name}.out
#  filename_err=${fortran_exec_name}.err

#  IFS='.' read -r python_module test_class <<<"$1"
#  set -x
#  python_cmd_get_fortran_args="from ${python_module} import get_fortran_args, ${test_class}; print(get_fortran_args(${test_class}))"
#  fortran_args=$(${MPIEXEC_CMD_SINGLE} ${CONTAINER_RUN_CMD} ${MOUNT_OPT} ${CONTAINER_IMAGE} ${CONTAINER_ENTRYPOINT_LOCAL} \
#  "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd ${MOUNT_CONTAINER_DIR} && \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
#    python -c \"${python_cmd_get_fortran_args}\"")
#  set +x

  BASH_CMD=()
#  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#    PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test \
#      python -m unittest -v $1")
#  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
#    export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test && \
#    export HPC_PREDICT_IO_TEST_FORTRAN_COMMAND=${FORTRAN_TEST_BINARY_PATH}/\$(python /src/hpc-predict/hpc-predict-io/test/test_sequential_args.py --test $1) && \
#    echo \"HPC_PREDICT_IO_TEST_FORTRAN_COMMAND = \${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}\" && \
#    cd ${MOUNT_CONTAINER_DIR} && \
#    eval \${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND}")
##    source <(echo \${HPC_PREDICT_IO_TEST_FORTRAN_COMMAND})")
##  BASH_CMD+=("source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && \
##    export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python:/src/hpc-predict/hpc-predict-io/test && \
##    cd /src/hpc-predict/hpc-predict-io/test && \
##    export HPC_PREDICT_IO_TEST_FORTRAN_EXEC=\$(python  test_sequential_args.py --test $1 --type exec) && \
##    export HPC_PREDICT_IO_TEST_FORTRAN_ARGS=(\$(python test_sequential_args.py --test $1 --type args)) && \
##    export HPC_PREDICT_IO_TEST_FORTRAN_OUT=\$( python  test_sequential_args.py --test $1 --type out) && \
##    export HPC_PREDICT_IO_TEST_FORTRAN_ERR=\$( python  test_sequential_args.py --test $1 --type err) && \
##    cd ${MOUNT_CONTAINER_DIR} && \
##    set -x && \
##    ${FORTRAN_TEST_BINARY_PATH}/\${HPC_PREDICT_IO_TEST_FORTRAN_EXEC} \${HPC_PREDICT_IO_TEST_FORTRAN_ARGS[@]} 1> \${HPC_PREDICT_IO_TEST_FORTRAN_OUT} 2> \${HPC_PREDICT_IO_TEST_FORTRAN_ERR}")

  BASH_CMD+=("/src/hpc-predict/hpc-predict-io/test/test_sequential_mr_io_driver_python.sh $1")
  BASH_CMD+=("/src/hpc-predict/hpc-predict-io/test/test_sequential_mr_io_driver_fortran.sh $1")

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

