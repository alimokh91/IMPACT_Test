import os
import glob
import argparse
import logging

# Parse data input and output directories
def parse_args():
    #input_path  = "./bern_data_experiments_source/"
    #output_path = "./bern_data_experiments_hpc_predict/"
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from preprocessed HDF5-files (Bernese experimental dataset).')
    parser.add_argument('--input', type=str, required=True,
                    help='Directory containing experimental data from Bern (numpy files with coordinates/velocity)')
    parser.add_argument('--output', type=str,  required=True,
                    help='Output directory for HDF5 files')
    parser.add_argument('--log', type=str, default="warn", help="Logging level")
    return parser.parse_args()


args = parse_args()
logging.basicConfig(level=args.log.upper())
output_filename = args.output + '/bern_experimental_dataset_flow_mri.h5'
if not os.path.exists(args.input):
    raise RuntimeError("The path {} does not exist. Exiting...".format(args.input))
#if os.path.exists(output_filename):
#    raise RuntimeError("The file {} exists already. Exiting...".format(output_filename))
if not os.path.exists(args.output):
    os.makedirs(args.output)
#import pdb
#pdb.set_trace()

# Define time slices (TODO: should be read from a metainformation file e.g. in JSON)
#time_slices = [{"time": 0., "files": glob.glob(args.input + '/*masked.npy')}]
time_slices = [{"time": 0.05, "files": glob.glob(args.input + '/time=0.05/*masked.npy')},
{"time": 0.06, "files": glob.glob(args.input + '/time=0.06/*masked.npy')}]
exp_protocol={"time_slices":time_slices, "period":0.7}
import json
with open("exp_protocol.json",'w') as exp_protocol_file:
    json.dump(exp_protocol,exp_protocol_file)
