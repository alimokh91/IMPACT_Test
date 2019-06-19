#!/bin/bash

set -euxo pipefail

PYTHONPATH=python python test/test_sequential_mr_io.py
PYTHONPATH=python python test/test_sequential_mr_io_bidirectional.py
