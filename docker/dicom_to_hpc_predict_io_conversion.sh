#!/bin/bash

set -euxo pipefail

if [ "$#" -ne 3 ]; then
    echo "Wrong number of arguments "$#" - expected mri_data_root, mri_sample_id, output_root."
    exit 1
fi

set -x
mri_data_root="$1" # "input_data/original/mri/MRT Daten Bern tar"
output_root="$2"   # "input_data/preprocessed/mri/MRT Daten Bern tar/${time_stamp_host}/output"
mri_sample_id="$3" # 10
set +x

export PYTHONPATH=/src/hpc-predict/hpc-predict-io/python
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate

set -x

mkdir -p "/hpc-predict-data/decrypt/${output_root}/${mri_sample_id}" || true

papermill -p hpc_predict_data_root /hpc-predict-data -p mri_data_root "${mri_data_root}" -p output_root "${output_root}" -y "{mri_samples:[${mri_sample_id}]}" "/src/hpc-predict/hpc-predict-io/mri_datasource/freiburg/freiburg_dicom_reader.ipynb" "/hpc-predict-data/decrypt/${output_root}/${mri_sample_id}/freiburg_dicom_reader_papermilled.ipynb"
