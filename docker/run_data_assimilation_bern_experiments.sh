#!/bin/bash

set -euo pipefail

HPC_PREDICT_IMPACT_IMAGE=${HPC_PREDICT_IMPACT_IMAGE:-'cscs-ci/hpc-predict/impact/deploy'}
HPC_PREDICT_DATA_DIR=$(realpath $1)

if [ "$#" -eq 2 ]; then
    time_stamp_host="$2"
else
    time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
fi

if [ ! -f "${HPC_PREDICT_DATA_DIR}/input_data/preprocessed/bern_experiments/v1/2020-01-01_00-00-00_dario/bern_experimental_dataset_segmented_flow_mri.h5" ]; then

    if [ ! -f "${HPC_PREDICT_DATA_DIR}/input_data/original/bern_experiments/v1/2020-01-01_00-00-00_dario/bern_exp_metadata.json" ]; then
        echo "Bern experiments original data not found - fetching now..."
        shell_command=$(printf "%s" \
            "set -x && " \
            "/hpc-predict-data/../fetch_scripts/fetch_bern_experiments_assimilation_in.sh " \
            "/hpc-predict-data/input_data/original/bern_experiments/v1/2020-01-01_00-00-00_dario ")

        set -x
        docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}"
        set +x
    fi

    echo "Bern experiments in HDF5 format of HPC-PREDICT-IO not found - converting original data now..."

    shell_command=$(printf "%s" \
        "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
        "set -x && " \
        "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python -u " \
        "/src/hpc-predict/hpc-predict-io/mri_datasource/convert_bern_expt_to_hpc_predict_with_downsampling.py " \
        "--input /hpc-predict-data/input_data/original/bern_experiments/v1/2020-01-01_00-00-00_dario/bern_exp_metadata.json " \
        "--output /hpc-predict-data/input_data/preprocessed/bern_experiments/v1/2020-01-01_00-00-00_dario/bern_experimental_dataset_segmented_flow_mri.h5 " \
        "--downsample 4" \
        "--single_rep -1")

    set -x
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}"
    set +x
fi

relative_output_directory="impact/bern_experiments/v1/${time_stamp_host}"
host_output_directory="${HPC_PREDICT_DATA_DIR}/${relative_output_directory}"
container_output_directory="/hpc-predict-data/${relative_output_directory}"

mkdir -p "${host_output_directory}"
echo "Host output directory: ${host_output_directory}"

# Generate config file

echo "Generate config.txt required by IMPACT based on MRI data"

shell_command=$(printf "%s" \
    "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
    "set -x && " \
    "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python -u " \
    "-m mr_io_impact_config " \
    "--input-mri /hpc-predict-data/input_data/preprocessed/bern_experiments/v1/2020-01-01_00-00-00_dario/bern_experimental_dataset_segmented_flow_mri.h5 " \
    "--output-mri \"${container_output_directory}/bern_experimental_dataset_assimilation_results.h5\" " \
    "--sr 4 4 4 " \
    "--padding 0.2 0.2 0.2 " \
    "--tr 2 " \
    "--pulses 10 " \
    "--output \"${container_output_directory}/config.txt\" " \
    "--np 8")
    # removed due to redundancy: "--config /src/hpc-predict/hpc-predict-io/python/config.txt.j2 " \

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}"
set +x

# Run data assimilation

echo "Launching data assimilation in IMPACT"

# FIXME: Generic docker/sarus run command 

shell_command=$(printf "%s" \
    "cd \"${container_output_directory}\" && " \
    "mpirun " \
    "-np 8 " \
    "/src/hpc-predict/IMPACT/prog/impact.exe")

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash "${HPC_PREDICT_IMPACT_IMAGE}" -c "${shell_command}" 
set +x
