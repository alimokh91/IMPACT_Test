#!/bin/bash

set -uxo pipefail

rm  *.mod \
    *.o 

cd test
./clean.sh
cd ..
