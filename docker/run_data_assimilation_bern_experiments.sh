#!/bin/bash

set -euo pipefail

HPC_PREDICT_DATA_DIR=$(realpath $1)

if [ ! -f "${HPC_PREDICT_DATA_DIR}/input_data/preprocessed/bern_experiments/bern_experimental_dataset_segmented_flow_mri.h5" ]; then

    if [ ! -f "${HPC_PREDICT_DATA_DIR}/input_data/original/bern_experiments/bern_exp_metadata.json" ]; then
        echo "Bern experiments original data not found - fetching now..."
        shell_command=$(printf "%s" \
            "set -x && " \
            "/src/hpc-predict/hpc-predict-io/mri_datasource/fetch_bern_experimental_data.sh " \
            "/hpc-predict-data/input_data/original/bern_experiments ")

        set -x
        docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}"
        set +x
    fi

    echo "Bern experiments in HDF5 format of HPC-PREDICT-IO not found - converting original data now..."

    shell_command=$(printf "%s" \
        "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
        "set -x && " \
        "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python " \
        "/src/hpc-predict/hpc-predict-io/mri_datasource/convert_bern_expt_to_hpc_predict.py " \
        "--input /hpc-predict-data/input_data/original/bern_experiments/bern_exp_metadata.json " \
        "--output /hpc-predict-data/input_data/preprocessed/bern_experiments/bern_experimental_dataset_segmented_flow_mri.h5 ")

    set -x
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}"
    set +x
fi

time_stamp=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
relative_output_directory="impact/bern_experiments/${time_stamp}"
host_output_directory="${HPC_PREDICT_DATA_DIR}/${relative_output_directory}"
container_output_directory="/hpc-predict-data/${relative_output_directory}"

mkdir -p "${host_output_directory}"
echo "Host output directory: ${host_output_directory}"

# Generate config file

echo "Generate config.txt required by IMPACT based on MRI data"

shell_command=$(printf "%s" \
    "source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && " \
    "set -x && " \
    "PYTHONPATH=/src/hpc-predict/hpc-predict-io/python python " \
    "-m mr_io_impact_config " \
    "--input-mri /hpc-predict-data/input_data/preprocessed/bern_experiments/bern_experimental_dataset_segmented_flow_mri.h5 " \
    "--output-mri ${container_output_directory}/bern_experimental_dataset_assimilation_results.h5 " \
    "--sr 2 2 2 " \
    "--padding 0.5 0.5 0.5 " \
    "--tr 2 " \
    "--config /src/hpc-predict/hpc-predict-io/python/config.txt.j2 " \
    "--output ${container_output_directory}/config.txt " \
    "--np 4")

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}"
set +x

# Run data assimilation

echo "Launching data assimilation in IMPACT"

# FIXME: Generic docker/sarus run command 

shell_command=$(printf "%s" \
    "cd ${container_output_directory} && " \
    "mpirun " \
    "-np 4 " \
    "/src/hpc-predict/IMPACT/prog/impact_debug.exe")

set -x
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash lukasgd/hpc-predict:impact-deploy -c "${shell_command}" 
set +x
