cmake_minimum_required(VERSION 3.8)
project(hpc-predict-io VERSION 0.1.0 LANGUAGES Fortran)

try_compile(CAN_COMPILE_MPI ${CMAKE_CURRENT_BINARY_DIR}
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/tests/test_mpi.f90)
if (CAN_COMPILE_MPI)
    message(STATUS "Fortran compiler can compile a simple MPI program")
else()
    message(FATAL_ERROR "Fortran compiler is unable to compile an MPI program.
    Please set the variable FC to an MPI enabled fortran compiler, e.g.
    mpifort")
endif()

add_subdirectory(fortran)

enable_testing()
add_test(NAME sequential COMMAND ${CMAKE_COMMAND} -E env
    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_sequential_mr_io.py)
add_test(NAME sequential_bidirectional COMMAND ${CMAKE_COMMAND} -E env
    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_sequential_mr_io_bidirectional.py)
add_test(NAME parallel COMMAND ${CMAKE_COMMAND} -E env
    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_parallel_mr_io.py)
add_test(NAME parallel_bidirectional COMMAND ${CMAKE_COMMAND} -E env
    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_parallel_mr_io_bidirectional.py)
add_test(NAME impact_integration WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND ${CMAKE_COMMAND} -E env
    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_impact_integration.py)
#add_test(NAME data_integration WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND ${CMAKE_COMMAND} -E env
#    FORTRAN_TEST_BINARY_PATH=${CMAKE_INSTALL_PREFIX}/bin/test/
#    PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
#    python ${CMAKE_CURRENT_SOURCE_DIR}/test/test_data_integration.py)

