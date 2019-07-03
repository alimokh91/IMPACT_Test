#!/bin/bash

if [ ! -d "./data_experiments" ]; then
  echo "Downloading files and extracting them to $(pwd)/data_experiments directory." 
  wget -O data_experiments.zip https://www.dropbox.com/sh/dv94nnhlt9veo9m/AAABxGfQw9i_nLDNbTJPNqT1a?dl=1
  unzip data_experiments.zip -d data_experiments
  rm data_experiments.zip
else
  echo "./data_experiments directory already exists - no files downloaded."
fi
