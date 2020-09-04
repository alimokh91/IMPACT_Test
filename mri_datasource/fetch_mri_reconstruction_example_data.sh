#!/bin/bash

set -eoux pipefail

DATASOURCE_DIR=./mri_reconstruction_example_source
DATASOURCE_FILENAME=mri_reconstruction_example.mat

if [ ! -d $DATASOURCE_DIR ]; then
  echo "Downloading files and extracting them to $(pwd)/${DATASOURCE_DIR} directory." 
  mkdir -p ${DATASOURCE_DIR}
  wget -O ${DATASOURCE_DIR}/${DATASOURCE_FILENAME} https://polybox.ethz.ch/index.php/s/xCPvKs2vCrCPaR9/download
else
  echo "${DATASOURCE_DIR} directory already exists - no files downloaded."
fi
