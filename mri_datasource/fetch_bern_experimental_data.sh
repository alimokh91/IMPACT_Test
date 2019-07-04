#!/bin/bash

if [ ! -d "./bern_data_experiments_source" ]; then
  echo "Downloading files and extracting them to $(pwd)/data_experiments directory." 
  wget -O bern_data_experiments.zip https://www.dropbox.com/sh/dv94nnhlt9veo9m/AAABxGfQw9i_nLDNbTJPNqT1a?dl=1
  unzip bern_data_experiments.zip -d bern_data_experiments_source
  rm bern_data_experiments.zip
else
  echo "./data_experiments directory already exists - no files downloaded."
fi
