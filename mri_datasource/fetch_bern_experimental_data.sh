#!/bin/bash

set -eoux pipefail

DATASOURCE_DIR=./bern_data_experiments_source
DATASOURCE_FILENAME=bern_data_experiments.zip

if [ ! -d $DATASOURCE_DIR ]; then
  echo "Downloading files and extracting them to $(pwd)/${DATASOURCE_DIR} directory." 
  wget -O ${DATASOURCE_FILENAME} https://www.dropbox.com/s/ko4uv8tpee2o4pe/exp_data_all_timeslices.zip?dl=0
  unzip ${DATASOURCE_FILENAME} -d ${DATASOURCE_DIR}
  rm ${DATASOURCE_FILENAME}
else
  echo "${DATASOURCE_DIR} directory already exists - no files downloaded."
fi
