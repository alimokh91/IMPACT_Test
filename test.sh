#!/bin/bash

set -euxo pipefail

PYTHONPATH=python python test/test_sequential_mr_io.py
PYTHONPATH=python python test/test_sequential_mr_io_bidirectional.py
PYTHONPATH=python python test/test_parallel_mr_io.py
PYTHONPATH=python python test/test_parallel_mr_io_bidirectional.py
PYTHONPATH=python python test/test_impact_input.py
PYTHONPATH=python python test/test_impact_input_padding.py
