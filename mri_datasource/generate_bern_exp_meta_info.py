import os
import glob
import argparse
import logging
import re
import numpy as np
import pdb
import json

# Parse data input and output directories
def parse_args():
    #input_path  = "./bern_data_experiments_source/"
    #output_file = "bern_data_experiments_metadata.json"
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from preprocessed HDF5-files (Bernese experimental dataset).')
    parser.add_argument('--input', type=str, required=True,
                    help='Directory containing experimental data from Bern (numpy files with coordinates/velocity)')
    parser.add_argument('--output', type=str,  required=True,
                    help='Output name of json file')
    parser.add_argument('--log', type=str, default="warn", help="Logging level")
    return parser.parse_args()



args = parse_args()
logging.basicConfig(level=args.log.upper())
if not os.path.exists(args.input):
    raise RuntimeError("The path {} does not exist. Exiting...".format(args.input))
if os.path.exists(args.output):
    raise RuntimeError("The file {} exists already. Exiting...".format(args.output))

# Read from input path the path to each phase folder and sort them
paths_to_phases = glob.glob(args.input + 'exp_data_all_timeslices/time*/')
paths_to_phases.sort()
time_slices = {}

# create a dictionary where the key is the time of each phase and the value is an array containing paths to all the files in this phase
for p in range(0,len(paths_to_phases)):
    time = re.search('time=(.*?)/',paths_to_phases[p]).group(1)
    paths_to_files = glob.glob(paths_to_phases[p] + 'data/*masked.npy')
    paths_to_files.sort()
    time_slices[time] = paths_to_files

heart_cycle_period = 60/70
exp_protocol={"time_slices":time_slices,"heart_cycle_period":heart_cycle_period}

with open(args.output,'w') as exp_protocol_file:
    #json.dump(time_slices,exp_protocol_file)
    json.dump(exp_protocol, exp_protocol_file, indent=4, sort_keys=True)
