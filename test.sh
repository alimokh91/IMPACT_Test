#!/bin/bash

set -euxo pipefail

FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_sequential_mr_io.py
FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_sequential_mr_io_bidirectional.py
FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_parallel_mr_io.py
FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_parallel_mr_io_bidirectional.py
FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_impact_integration.py
#FORTRAN_TEST_BINARY_PATH=fortran/test/ PYTHONPATH=python python test/test_data_integration.py

