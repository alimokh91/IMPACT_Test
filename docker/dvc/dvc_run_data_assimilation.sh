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
container_input_file="${container_input_directory}/${CNN_SEGMENTER_INFERENCE_RESULT}"

if [ -f "${host_input_directory}" ]; then
    echo "Segmenter inference directory \"${host_input_directory}\" with input data for inference does not exist. Exiting..."
    exit 1
fi

relative_output_directory="impact/hpc_predict/v1/${time_stamp_host}"
host_output_directory="${HPC_PREDICT_DATA_DIR}/${relative_output_directory}"
container_output_directory="/hpc-predict-data/${relative_output_directory}"
container_output_file="${container_output_directory}/output/${IMPACT_DATA_ASSIMILATION_RESULT}"

shell_command=$(printf "%s" \
  " \"$(realpath --relative-to="${host_output_directory}" $(realpath $(dirname $0)))/dvc_check_code_consistency.sh\" \"${HPC_PREDICT_IMPACT_IMAGE}\" && " \
  " docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data \"${HPC_PREDICT_IMPACT_IMAGE}\" " \
  " /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_config.sh "${container_input_file}" "${container_output_file}" && " \
  " docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data \"${HPC_PREDICT_IMPACT_IMAGE}\" " \
  " /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_impact.sh "${container_output_directory}/output" ")



set -x
mkdir -p "${host_output_directory}/output"
cd "${host_output_directory}"
# Generate config file & run data assimilation
dvc run --no-exec -n "impact_data_assimilation_${time_stamp_host}" -d "${host_input_directory}" -o "${host_output_directory}/output" ${shell_command}
#dvc run -n "impact_data_assimilation_${time_stamp_host}" -d "${host_input_directory}" -o "${host_output_directory}/output" $(realpath --relative-to="${host_output_directory}" $(realpath $(dirname $0)))/dvc_check_code_consistency.sh "${HPC_PREDICT_IMPACT_IMAGE}" && mkdir -p "${host_output_directory}/output" && docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data "${HPC_PREDICT_IMPACT_IMAGE}" /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_config.sh "${container_input_file}" "${container_output_file}" && docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data "${HPC_PREDICT_IMPACT_IMAGE}" /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_impact.sh "${container_output_directory}/output"
#dvc run -n "impact_data_assimilation_${time_stamp_host}" -d "${host_input_directory}" -o "${host_output_directory}/output" mkdir -p "${host_output_directory}/output" && docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data "${HPC_PREDICT_IMPACT_IMAGE}" /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_config.sh "${container_input_file}" "${container_output_file}" && docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data "${HPC_PREDICT_IMPACT_IMAGE}" /src/hpc-predict/IMPACT/docker/dvc/data_assimilation_impact.sh "${container_output_directory}/output"
cd -
set +x

