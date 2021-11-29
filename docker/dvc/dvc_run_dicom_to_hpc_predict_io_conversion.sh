#!/bin/bash

set -euo pipefail

HPC_PREDICT_DATA_DIR=$(realpath $1)
MRI_SAMPLE_DIR=$2 # e.g. 3

if [ "$#" -eq 3 ]; then
    time_stamp_host="$3"
else
    time_stamp_host=$(date +'%Y-%m-%d_%H-%M-%S')_$(hostname)
fi

mri_data_root="input_data/original/mri/MRT Daten Bern tar"
output_root="input_data/preprocessed/mri/MRT Daten Bern tar/${time_stamp_host}"

# find main repository directory (this works only if this is a submodule of the main repository)
# MAIN_REPO_DIR="$(dirname $(realpath $0))"
# while [[ ! -x ${MAIN_REPO_DIR}/ci/fetch_containers.sh && ${MAIN_REPO_DIR} != "/" ]] ; do
#     MAIN_REPO_DIR=$(dirname "${MAIN_REPO_DIR}")
# done
# if [[ ! -x "${MAIN_REPO_DIR}"/data/container_scripts/run_docker_or_sarus.sh ]] ; then
#     echo "Could not find the script container_scripts/run_docker_or_sarus.sh. Is this a separate clone, and not a full clone of the main repository?"
#     exit 1
# fi

# Add decrypt prefix in case of a data directory with encryption
if [ -d "${HPC_PREDICT_DATA_DIR}/encrypt" ]; then
  echo "Identified ${HPC_PREDICT_DATA_DIR} as directory with encryption."
#   pgrep_encfs=$(pgrep --list-full encfs)
#   if echo "${pgrep_encfs}" | grep -q "${HPC_PREDICT_DATA_DIR}/decrypt"; then
    echo "encfs is running on ${HPC_PREDICT_DATA_DIR}/decrypt."
    encrypt_dir="encrypt/"
    decrypt_dir="decrypt/"
    dvc_dir="config/"
#   else
#     echo "encfs seems not to be running - must be launched first."
#     exit 1
#   fi
else
  echo "Identified ${HPC_PREDICT_DATA_DIR} as directory without encryption."
  encrypt_dir=""
  decrypt_dir=""
  dvc_dir=""
fi

relative_input_file="${mri_data_root}/${MRI_SAMPLE_DIR}.tar"
host_encrypt_input_file="${HPC_PREDICT_DATA_DIR}/${encrypt_dir}${relative_input_file}"
host_decrypt_input_directory="$(dirname "${HPC_PREDICT_DATA_DIR}/${decrypt_dir}${relative_input_file}")"
container_decrypt_input_directory="$(dirname "/hpc-predict-data/${decrypt_dir}${relative_input_file}")"

if [ ! -f "${host_encrypt_input_file}" ]; then
    echo "MRI data source file \"${host_encrypt_input_file}\" with input data for conversion to HDF5 does not exist. Exiting..."
    exit 1
fi

relative_output_directory="${output_root}/${MRI_SAMPLE_DIR}"
host_encrypt_output_directory="${HPC_PREDICT_DATA_DIR}/${encrypt_dir}${relative_output_directory}"
host_decrypt_output_directory="${HPC_PREDICT_DATA_DIR}/${decrypt_dir}${relative_output_directory}"
container_decrypt_output_directory="/hpc-predict-data/${decrypt_dir}${relative_output_directory}"

host_dvc_directory="${HPC_PREDICT_DATA_DIR}/${dvc_dir}${relative_output_directory}"

# TODO: set -e & check if encfs is running (cf. above)
set -x
mkdir -p "${host_dvc_directory}"
mkdir -p "${host_decrypt_output_directory}"
set +x

#   " \"$(realpath --relative-to="${host_dvc_directory}" $(realpath $(dirname $0)))/dvc_check_code_consistency.sh\" \"${HPC_PREDICT_IO_IMAGE}\" && " \

# # The following two lines would restrict data access to certain subdirectories, but currently failing (see comment below, despite using allow_root or allow_other with encfs) 
#   " -v \"\$(pwd)/$(realpath --relative-to="${host_dvc_directory}" "${host_decrypt_input_directory}"):${container_decrypt_input_directory}:ro\" " \
#   " -v \"\$(pwd)/$(realpath --relative-to="${host_dvc_directory}" "${host_decrypt_output_directory}"):${container_decrypt_output_directory}\" \"${HPC_PREDICT_IO_IMAGE}\" " \

#   "  -v \"\$(pwd)/$(realpath --relative-to="${host_dvc_directory}" ${HPC_PREDICT_DATA_DIR}):/hpc-predict-data\" \"${HPC_PREDICT_IO_IMAGE}\" " \

# hpc-predict/data/v1/config/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11$ dvc repro
# +-------------------------------------------+
# |                                           |
# |     Update available 1.11.16 -> 2.1.0     |
# |      Run `pip install dvc --upgrade`      |
# |                                           |
# +-------------------------------------------+
# 
# WARNING: stage: '../../../../../original/mri/MRT Daten Bern tar/dvc.yaml:data_fetch_freiburg_mri_data_to_tar_2021-05-14_18-36-50_ThinkPad-X1' is frozen. Its dependencies are not going to be reproduced.
# Stage '../../../../../original/mri/MRT Daten Bern tar/dvc.yaml:data_fetch_freiburg_mri_data_to_tar_2021-05-14_18-36-50_ThinkPad-X1' didn't change, skipping
# Running stage 'dicom_to_hpc_predict_io_testing' with command:
# 	 docker run --rm  -v "$(pwd)/../../../../../../../decrypt/input_data/original/mri/MRT Daten Bern tar:/hpc-predict-data/decrypt/input_data/original/mri/MRT Daten Bern tar:ro"  -v "$(pwd)/../../../../../../../decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11:/hpc-predict-data/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11" "lukasgd/hpc-predict:io-deploy"  /src/hpc-predict/hpc-predict-io/docker/dicom_to_hpc_predict_io_conversion.sh "input_data/original/mri/MRT Daten Bern tar" "input_data/preprocessed/mri/MRT Daten Bern tar/testing" "11" 
# docker: Error response from daemon: error while creating mount source path '/home/lukasd/src/hpc-predict/hpc-predict/data/v1/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11': mkdir /home/lukasd/src/hpc-predict/hpc-predict/data/v1/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11: file exists.
# ERROR: failed to reproduce 'dvc.yaml': failed to run:  docker run --rm  -v "$(pwd)/../../../../../../../decrypt/input_data/original/mri/MRT Daten Bern tar:/hpc-predict-data/decrypt/input_data/original/mri/MRT Daten Bern tar:ro"  -v "$(pwd)/../../../../../../../decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11:/hpc-predict-data/decrypt/input_data/preprocessed/mri/MRT Daten Bern tar/testing/11" "lukasgd/hpc-predict:io-deploy"  /src/hpc-predict/hpc-predict-io/docker/dicom_to_hpc_predict_io_conversion.sh "input_data/original/mri/MRT Daten Bern tar" "input_data/preprocessed/mri/MRT Daten Bern tar/testing" "11" , exited with 125

shell_command=$(printf "%s" \
  " docker run --rm " \
  "  -v \"\$(pwd)/$(realpath --relative-to="${host_dvc_directory}" ${HPC_PREDICT_DATA_DIR}):/hpc-predict-data\" \"${HPC_PREDICT_IO_IMAGE}\" " \
  " /src/hpc-predict/hpc-predict-io/docker/dicom_to_hpc_predict_io_conversion.sh \"${mri_data_root}\" \"${output_root}\" \"${MRI_SAMPLE_DIR}\" ")

set -x
cd "${host_dvc_directory}"
dvc run --no-exec -n "dicom_to_hpc_predict_io_${time_stamp_host}" -d "$(realpath --relative-to="${host_dvc_directory}" "${host_encrypt_input_file}")" -o "$(realpath --relative-to="${host_dvc_directory}" "${host_encrypt_output_directory}")" "${shell_command}"
cd -
set +x
