import os
import numpy as np
import argparse
import logging
import scipy.io as sio
from mr_io import FlowMRI # Requires adding ../python to PYTHONPATH

# Parse data input and output directories
def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from Jonas\' 4D flow example provided.')
    parser.add_argument('--input', type=str,
                    help='Directory containing data4Dflow_mean_std.mat (with mean and coordinate-wise std deviation of velocity)')
    parser.add_argument('--voxel-width', type=int, nargs=3, default=[2.5, 2.45, 2.48],
                    help='Spatial voxel width')
    parser.add_argument('--time-width', type=int, default=0.04,
                    help='Time slice width')
    parser.add_argument('--output', type=str,
                    help='Output name for HDF5 file')
    return parser.parse_args()

args = parse_args()

mat_data = sio.loadmat(args.input)

logging.info("Read mat-file %s with fields {}".format(*mat_data.keys()), args.input)

# Guess geometry parameters if not supplied

spatial_voxel_width = args.voxel_width 
if args.voxel_width is None:
    spatial_voxel_width = [25./mat_data['results_v'].shape[1]]*3
    logging.warning("Spatial voxel width not supplied - guessing it as isotropic {} m".format(spatial_voxel_width))

time_slice_width = args.time_width
if args.time_width is None:
    time_slice_width = 60./70/mat_data['results_v'].shape[3]
    logging.warning("Time slice width not supplied - guessing it as {} s.".format(time_slice_width))

geometry = [ np.array([ i*spatial_voxel_width[j] for i in range(mat_data['results_v'].shape[j]) ]) \
            for j in range(3) ]
time = np.array([ i*time_slice_width for i in range(mat_data['results_v'].shape[3]) ])
velocity_cov = np.zeros( mat_data['results_v'].shape + mat_data['results_v'].shape[-1:] )
mat_data['results_std'][ mat_data['results_std'] < 0. ] = 0. # filter negative std deviations (cf. Jonas' email: negative standard deviations occur in regions where we have little signal (e.g. the lungs))
for i in range(mat_data['results_v'].shape[-1]):
    velocity_cov[:,:,:,:,i,i] = mat_data['results_std'][:,:,:,:,i]**2
    
output_dir = os.path.dirname(os.path.realpath(args.output))
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
flow_mri = FlowMRI(geometry=geometry,
                                time=time,
                                intensity=mat_data['I'].real,
                                velocity_mean=mat_data['results_v'],
                                velocity_cov=velocity_cov)
flow_mri.write_hdf5(args.output)
