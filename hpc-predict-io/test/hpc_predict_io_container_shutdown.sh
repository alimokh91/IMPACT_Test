#!/bin/bash

if echo ${MPI_MASTER_HOST} | grep -i daint > /dev/null; then
    if ${MODULE_SARUS_LOADED}; then
        module unload sarus
    fi
fi
