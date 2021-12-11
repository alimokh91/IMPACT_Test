#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Wrong number of arguments "$#" - expected data-assimilation-output."
    exit 1
fi

container_output_directory="$1"

set -x
cd "${container_output_directory}" && mpirun -np 4 /src/hpc-predict/IMPACT/prog/impact_debug.exe
set +x
