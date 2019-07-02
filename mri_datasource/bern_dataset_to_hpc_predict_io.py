import h5py
import os
import numpy as np
import re
import argparse

from mr_io import SpaceTimeMRI

def main():
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from preprocessed HDF5-files (Bernese experimental dataset).')
    parser.add_argument('--input', type=str,
                    help='Directory containing HDF5 files with coordinates and mean/covariance of velocity field')
    parser.add_argument('--output', type=str,  default="mri_bern_experimental_dataset.h5",
                    help='Output configuration file for IMPACT')
    args = parser.parse_args()

    # Read existing data files
    with h5py.File(args.input + "/Coordinates.h5","r") as coordinates_f:
        x_coordinates = np.unique( coordinates_f["Coordinates"][()][0,:] )
        y_coordinates = np.unique( coordinates_f["Coordinates"][()][1,:] )
        z_coordinates = np.unique( coordinates_f["Coordinates"][()][2,:] )

        #import pdb;pdb.set_trace()
        #x = np.flip(coordinates_f["Coordinates"][()][0,:].reshape(1,len(z_coordinates),len(y_coordinates), len(x_coordinates) ).transpose(3,2,1,0), axis=1)
        #y = np.flip(coordinates_f["Coordinates"][()][1,:].reshape(1,len(z_coordinates),len(y_coordinates), len(x_coordinates) ).transpose(3,2,1,0), axis=1)
        #z = np.flip(coordinates_f["Coordinates"][()][2,:].reshape(1,len(z_coordinates),len(y_coordinates), len(x_coordinates) ).transpose(3,2,1,0), axis=1)
        
        assert(len(x_coordinates)*len(y_coordinates)*len(z_coordinates) == coordinates_f["Coordinates"][()].shape[1])
            
    mean_files =  [args.input + "/" + f for f in os.listdir(args.input) if re.match(r'Velocity_Mean_\d+.h5',f)]    
    means = np.ndarray(shape=(len(x_coordinates), len(y_coordinates), \
                             len(z_coordinates), len(mean_files), 3))
    for time_slice, fname in enumerate(mean_files):
        with h5py.File(fname,"r") as mean_f:
            means[:,:,:,time_slice,:] = np.flip( mean_f["Velocity_mean"][()].reshape(3, # correct unusual layout
                                        len(z_coordinates), 
                                        len(y_coordinates), 
                                        len(x_coordinates) ).transpose(3,2,1,0), axis=1)
        
    covariance_files =  [args.input + "/" + f for f in os.listdir(args.input) if re.match(r'Covariance_\d+.h5',f)]    
    assert(len(mean_files) == len(covariance_files))
    covariances = np.ndarray(shape=(len(x_coordinates), len(y_coordinates), \
                                    len(z_coordinates), len(covariance_files), 6))
    for time_slice, fname in enumerate(covariance_files):
        with h5py.File(fname,"r") as covariance_f:
            covariances[:,:,:,time_slice,:] = np.flip( covariance_f["Covariance"][()].reshape(6, # correct unusual layout
                                        len(z_coordinates), 
                                        len(y_coordinates), 
                                        len(x_coordinates) ).transpose(3,2,1,0), axis=1)
    
    # Write MRI using hpc-predict-io
    mri = SpaceTimeMRI(geometry=[x_coordinates,
                                 y_coordinates,
                                 z_coordinates],
                       time=np.array([i for i in range(len(mean_files))]),
                       voxel_feature=means)
    
    mri.write_hdf5(args.output)


if __name__ == '__main__':
    main()