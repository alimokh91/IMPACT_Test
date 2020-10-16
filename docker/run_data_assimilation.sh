#!/bin/bash

set -euo pipefail

HPC_PREDICT_IMPACT_IMAGE=${HPC_PREDICT_IMPACT_IMAGE:-'cscs-ci/hpc-predict/impact/deploy'}
HPC_PREDICT_DATA_DIR=$(realpath $1)
CNN_SEGMENTER_INFERENCE_OUTPUT=$2

if [ "$#" -eq 3 ]; then
    time_stamp_host="$3"
else
    time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
fi

CNN_SEGMENTER_INFERENCE_RESULT=$(basename ${CNN_SEGMENTER_INFERENCE_OUTPUT})
CNN_SEGMENTER_INFERENCE_DIR=$(dirname ${CNN_SEGMENTER_INFERENCE_OUTPUT})

IMPACT_DATA_ASSIMILATION_RESULT="${CNN_SEGMENTER_INFERENCE_RESULT%".h5"}_data_assimilated.h5"

relative_input_directory="segmenter/cnn_segmenter/hpc_predict/v1/inference/${CNN_SEGMENTER_INFERENCE_DIR}"
host_input_directory="${HPC_PREDICT_DATA_DIR}/${relative_input_directory}"
container_input_directory="/hpc-predict-data/${relative_input_directory}"

if [ -f "${host_input_directory}" ]; then
    echo "Segmenter inference directory \"${host_input_directory}\" with input data for inference does not exist. Exiting..."
    exit 1
fi

relative_output_directory="impact/hpc_predict/v1/${time_stamp_host}"
host_output_directory="${HPC_PREDICT_DATA_DIR}/${relative_output_directory}"
container_output_directory="/hpc-predict-data/${relative_output_directory}"

mkdir -p "${host_output_directory}"
echo "Host output directory: \"${host_output_directory}\""

# Generate config file

#FIXME: input mri path to output of segmenter,  

shell_command=$(printf "%s" \
    "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
    "set -x && " \
    "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python -u " \
    "-m mr_io_impact_config " \
    "--input-mri \"${container_input_directory}/${CNN_SEGMENTER_INFERENCE_RESULT}\" " \
    "--output-mri \"${container_output_directory}/${IMPACT_DATA_ASSIMILATION_RESULT}\" " \
    "--sr 2 2 2 " \
    "--padding 0.5 0.5 0.5 " \
    "--tr 2 " \
    "--output \"${container_output_directory}/config.txt\" " \
    "--np 4")
    # removed due to redundancy: "--config /src/hpc-predict/hpc-predict-io/python/config.txt.j2 " \

set -x
docker run  -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}"
set +x

# Run data assimilation

shell_command=$(printf "%s" \
    "cd \"${container_output_directory}\" && " \
    "mpirun " \
    "-np 4 " \
    "/src/hpc-predict/IMPACT/prog/impact_debug.exe")

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}" 
set +x

