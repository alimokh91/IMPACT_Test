#!/bin/bash

set -eoux pipefail

DATASOURCE_DIR=./bern_data_experiments_source
DATASOURCE_FILENAME=bern_data_experiments.zip
DATASOURCE_METADATA_FILENAME=bern_exp_metadata.json

if [ ! -d $DATASOURCE_DIR ]; then
  mkdir -p ${DATASOURCE_DIR}

  echo "Downloading json file."
  wget -O ${DATASOURCE_METADATA_FILENAME} https://www.dropbox.com/s/otbu709nc6k8fsq/bern_exp_metadata.json?dl=0 
  mv ${DATASOURCE_METADATA_FILENAME} ${DATASOURCE_DIR}/ 

  echo "Downloading files and extracting them to $(pwd)/${DATASOURCE_DIR} directory." 
  wget -O ${DATASOURCE_FILENAME} https://www.dropbox.com/s/62gisr19w247thx/exp_data_all_timeslices.zip?dl=0 
  unzip ${DATASOURCE_FILENAME} -d ${DATASOURCE_DIR}
  rm ${DATASOURCE_FILENAME}
else
  echo "${DATASOURCE_DIR} directory already exists - no files downloaded."
fi
