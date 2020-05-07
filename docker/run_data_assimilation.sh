#!/bin/bash

set -euxo pipefail

HPC_PREDICT_DATA_DIR=$(realpath $1)
CNN_SEGMENTER_INFERENCE_OUTPUT=$2

CNN_SEGMENTER_INFERENCE_RESULT=$(basename ${CNN_SEGMENTER_INFERENCE_OUTPUT})
CNN_SEGMENTER_INFERENCE_DIR=$(dirname ${CNN_SEGMENTER_INFERENCE_OUTPUT})

IMPACT_DATA_ASSIMILATION_RESULT=${CNN_SEGMENTER_INFERENCE_RESULT%".h5"}_data_assimilated.h5

relative_input_directory=segmenter/cnn_segmenter/hpc_predict/inference/${CNN_SEGMENTER_INFERENCE_DIR}
host_input_directory="${HPC_PREDICT_DATA_DIR}/${relative_input_directory}"
container_input_directory="/hpc-predict-data/${relative_input_directory}"

if [ -f ${host_input_directory} ]; then
    echo "Segmenter inference directory \"${host_input_directory}\" with input data for inference does not exist. Exiting..."
    exit 1
fi

time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
relative_output_directory="impact/hpc_predict/${time_stamp_host}"
host_output_directory="${HPC_PREDICT_DATA_DIR}/${relative_output_directory}"
container_output_directory="/hpc-predict-data/${relative_output_directory}"

mkdir -p ${host_output_directory}
echo "Host output directory: ${host_output_directory}"

# Generate config file

#FIXME: input mri path to output of segmenter,  

shell_command=$(printf "%s" \
    "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
    "set -x && " \
    "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python " \
    "-m mr_io_impact_config " \
    "--input-mri ${container_input_directory}/${CNN_SEGMENTER_INFERENCE_RESULT} " \
    "--output-mri ${container_output_directory}/${IMPACT_DATA_ASSIMILATION_RESULT} " \
    "--sr 2 2 2 " \
    "--padding 0.5 0.5 0.5 " \
    "--tr 2 " \
    "--config /src/hpc-predict/hpc-predict-io/python/config.txt.j2 " \
    "--output ${container_output_directory}/config.txt " \
    "--np 4")

docker run  -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}"

# Run data assimilation

shell_command=$(printf "%s" \
    "cd ${container_output_directory} && " \
    "mpirun " \
    "-np 4 " \
    "/src/hpc-predict/IMPACT/prog/impact_debug.exe")

docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}" 

