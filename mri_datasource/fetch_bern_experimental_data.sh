#!/bin/bash

set -eoux pipefail

DATASOURCE_DIR=./bern_data_experiments_source
DATASOURCE_FILENAME=bern_data_experiments.zip

if [ ! -d $DATASOURCE_DIR ]; then
  echo "Downloading files and extracting them to $(pwd)/${DATASOURCE_DIR} directory." 
  wget -O ${DATASOURCE_FILENAME} https://www.dropbox.com/sh/dv94nnhlt9veo9m/AAABxGfQw9i_nLDNbTJPNqT1a?dl=1
  unzip ${DATASOURCE_FILENAME} -d ${DATASOURCE_DIR}
  rm ${DATASOURCE_FILENAME}
else
  echo "${DATASOURCE_DIR} directory already exists - no files downloaded."
fi
