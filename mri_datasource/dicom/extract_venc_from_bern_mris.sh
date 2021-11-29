#!/bin/bash

set -eou pipefail

if [ $# -ne 1 ]
  then
    echo "Error: Expecting HPC_PREDICT_DATA_DIR as single argument."
fi

HPC_PREDICT_DATA_DIR=$(realpath "$1")

# Add decrypt prefix in case of a data directory with encryption
if [ -d "${HPC_PREDICT_DATA_DIR}/decrypt" ]; then
  echo "Identified ${HPC_PREDICT_DATA_DIR} as directory with encryption."
  pgrep_encfs=$(pgrep --list-full encfs)
  if echo "${pgrep_encfs}" | grep -q "${HPC_PREDICT_DATA_DIR}/decrypt"; then
    echo "encfs is running on data directory - fetching to ${HPC_PREDICT_DATA_DIR}/decrypt."
    HPC_PREDICT_DATA_DIR="${HPC_PREDICT_DATA_DIR}/decrypt"
  else
    echo "encfs seems not to be running - must be launched first."
    exit 1
  fi
else
  echo "Identified ${HPC_PREDICT_DATA_DIR} as directory without encryption."
fi

# This could probably just be added to the Bern MRI dataset download script
venc_xy=$(gdcmdump -r --csa -i "${HPC_PREDICT_DATA_DIR}/input_data/original/mri/MRT Daten Bern/127/10005583" 2>/dev/null | grep 'sWiPMemBlock\.adFree\[8\]' | uniq | awk '{print $3;} ')
venc_z=$( gdcmdump -r --csa -i "${HPC_PREDICT_DATA_DIR}/input_data/original/mri/MRT Daten Bern/127/10005583" 2>/dev/null | grep 'sWiPMemBlock\.adFree\[9\]' | uniq | awk '{print $3;} ')

set -x
printf "{\n \"venc_xy\": ${venc_xy},\n \"venc_z\": ${venc_z},\n \"venc_units\": \"m/s\"\n}\n" #> "${HPC_PREDICT_DATA_DIR}/input_data/original/mri/MRT Daten Bern/venc.json"
