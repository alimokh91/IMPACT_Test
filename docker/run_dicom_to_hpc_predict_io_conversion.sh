#!/bin/bash

set -euxo pipefail

HPC_PREDICT_DATA_DIR=$(realpath $1)
MRI_SAMPLE_DIR=$2

if [ "$#" -eq 3 ]; then
    time_stamp_host="$3"
else
    time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
fi

mri_data_root="input_data/original/mri/MRT Daten Bern tar"
output_root="input_data/preprocessed/mri/MRT Daten Bern tar/${time_stamp_host}"

# Process a single MRI sample directory (could be part of a loop over MRI sample directories, all within the same mri_data_root/output_root)
docker run --rm -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data ${HPC_PREDICT_IO_IMAGE} /src/hpc-predict/hpc-predict-io/docker/dicom_to_hpc_predict_io_conversion.sh "${mri_data_root}" "${output_root}" "${MRI_SAMPLE_DIR}"

# to be deleted
# docker run --rm -v ${HPC_PREDICT_DATA_DIR}:/hpc-predict-data --entrypoint bash ${HPC_PREDICT_IO_IMAGE} -c "export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python; export LC_ALL=C.UTF-8; export LANG=C.UTF-8; source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && mkdir -p \"/hpc-predict-data/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/${MRI_SAMPLE_DIR}\" && papermill -p hpc_predict_data_root /hpc-predict-data -y '{mri_samples:[${MRI_SAMPLE_DIR}]}' \"/src/hpc-predict/hpc-predict-io/mri_datasource/dicom/freiburg_dicom_reader.ipynb\" \"/hpc-predict-data/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/${MRI_SAMPLE_DIR}/freiburg_dicom_reader_papermilled.ipynb\""
