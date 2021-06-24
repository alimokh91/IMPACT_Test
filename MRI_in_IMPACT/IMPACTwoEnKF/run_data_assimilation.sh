#!/bin/bash

set -euo pipefail

if [ "$#" -eq 3 ]; then
    time_stamp_host="$3"
else
    time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
fi

DATA_DIR="/project/c12/pcorso/hpc-predict/data/v1"
IMPACT_DIR="/project/c12/pcorso/hpc-predict/IMPACT"
cd $DATA_DIR
mkdir -p IMPACT
IMPACT_DATA_DIR="${DATA_DIR}/IMPACT"
echo "$IMPACT_DATA_DIR"
BASENAME_H5_FILE=$(basename -s .h5 -- ${DATA_DIR}/*.h5)
echo "$BASENAME_H5_FILE"
IMPACT_DATA_ASSIMILATION_RESULT="${BASENAME_H5_FILE%".h5"}_data_assimilated.h5"
echo "$IMPACT_DATA_ASSIMILATION_RESULT"

# Generate config file 

# Activate beforehand python env in hpc-predict-io/python

echo "${DATA_DIR}/${BASENAME_H5_FILE}.h5" 
echo "${IMPACT_DATA_DIR}/${IMPACT_DATA_ASSIMILATION_RESULT}"

PYTHONPATH=/project/c12/pcorso/hpc-predict/hpc-predict-io/python && \
cd $PYTHONPATH && \
python mr_io_impact_config.py \
--input-mri "${DATA_DIR}/${BASENAME_H5_FILE}.h5" \
--output-mri "${IMPACT_DATA_DIR}/${IMPACT_DATA_ASSIMILATION_RESULT}" \
--sr 2 2 2 \
--padding 0.2 0 0 \
--tr 2 \
--pulses 1 \
--output "${IMPACT_DATA_DIR}/config.txt" \
--np 36
    # removed due to redundancy: "--config /src/hpc-predict/hpc-predict-io/python/config.txt.j2 " \


# Run data assimilation on Daint

#shell_command=$(printf "%s" \
#    "cd \"${container_output_directory}\" && " \
#    "mpirun " \
#    "-np 16 " \
#    "/src/hpc-predict/IMPACT/prog/impact.exe")


