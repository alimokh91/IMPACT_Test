import os
import h5py
import numpy as np
import glob
import argparse
import logging
import functools
from mr_io import SegmentedFlowMRI # Requires adding ../python to PYTHONPATH
import json
import pdb
import math

# Parse data input and output directories
def parse_args():
    #input_path  = "./bern_data_experiments_source/"
    #output_path = "./bern_data_experiments_hpc_predict/"
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from preprocessed HDF5-files (Bernese experimental dataset).')
    parser.add_argument('--input', type=str, required=True,
                    help='JSON file containing metadata about the  experimental data from Bern (numpy files with coordinates/velocity)')
    parser.add_argument('--output', type=str,  required=True,
                    help='Name of output file in hpc-predict-io HDF5 format')
    parser.add_argument('--downsample', default=1, type=int, required=False,
                    help='factor for downsampling the original data')		
    parser.add_argument('--single_rep', default=-1, type=int, required=False,
                    help='save velocity of a specific repetition instead of the phase-averaged velocity')		
    parser.add_argument('--log', type=str, default="warn", help="Logging level")
    return parser.parse_args()


args = parse_args()
logging.basicConfig(level=args.log.upper())
output_filename = os.path.realpath(args.output)
if args.downsample != 1 and args.downsample < 1:
    print('Downsampling factor cannot be smaller than one! Skip downsampling...')
    dw = 1
else:        
    print('downsampling by a factor of :',args.downsample**3)
    dw = args.downsample
if args.single_rep < 0 or args.single_rep > 23:
    print('saving mean velocity')
    single_rep = -1
    do_for_single_rep = False
else:        
    print('saving velocity of repetition:',args.single_rep)
    single_rep = args.single_rep
    do_for_single_rep = True
if not os.path.exists(args.input):
    raise RuntimeError("The file {} does not exist. Exiting...".format(args.input))
if os.path.exists(args.output):
    raise RuntimeError("The file {} exists already. Exiting...".format(args.output))
if args.output[-3:] != '.h5':
    raise RuntimeError("The file {} does not end with '.h5'. Exiting...".format(args.output))
if not os.path.exists(os.path.dirname(os.path.realpath(args.output))):
    os.makedirs(os.path.dirname(os.path.realpath(args.output)))

output_filename = os.path.realpath(args.output)
os.chdir(os.path.dirname(args.input))
set_boundary = True
# open JSON file with all the meta data of the experiment (time(=key), path to files(=values for each key))
with open(os.path.basename(args.input),'r') as exp_protocol_file:
    exp_protocol = json.load(exp_protocol_file)

time_slices = exp_protocol["time_slices"]
heart_cycle_period = exp_protocol["heart_cycle_period"]
time_str = list(time_slices.keys())
time_str.sort()

def read_velocity_time_slice(flist):
    # file list is the list of the n files representing the n repetitions of a specific phase.

    fname = flist[0]
    data = np.load(fname)
    data[np.isnan(data)] = 0.

    X = np.unique(data[:,0])
    Y = np.unique(data[:,1])
    Z = np.unique(data[:,2])
    I = np.size(X)
    J = np.size(Y)
    K = np.size(Z)

    coords_coarse = downsample_data(data[:,0:3], I,J,K)
    A = np.array(coords_coarse[:,0])
    x_size = A.size

    #For defining size of covariance matrix
    data_coarse = np.zeros((x_size,6))
    data_coarse[:,0:3] = coords_coarse

    # calculate mean
    #print('\ncalculate mean velocity field\n')
    k = 0;
    for fname in flist:

        # load data instantaneous field
        data = np.load(fname)
        data[np.isnan(data)] = 0.

        velocity_coarse = downsample_data(data[:,3:6], I,J,K)
        data_coarse[:,3:6] = velocity_coarse

        if k == 0:
            # allocate
            FU_mean = data_coarse[:,3].copy()
            FV_mean = data_coarse[:,4].copy()
            FW_mean = data_coarse[:,5].copy()
        else:
            # increment
            FU_mean = FU_mean + data_coarse[:,3]
            FV_mean = FV_mean + data_coarse[:,4]
            FW_mean = FW_mean + data_coarse[:,5]

        # iterations
        k = k+1
        #print(k)
    print('calculated mean')
    FU_mean = FU_mean / k
    FV_mean = FV_mean / k
    FW_mean = FW_mean / k

    # calculate rms and variance and covariance
    #--------------------------------------
    # old version w/o downsampling
    #--------------------------------------
    k = 0
    R = np.zeros((x_size,9))
    for fname in flist:

        #load data
        data = np.load(fname)
        data[np.isnan(data)] = 0.

        velocity_coarse = downsample_data(data[:,3:6], I,J,K)
        data_coarse[:,3:6] = velocity_coarse

        if k == 0:
            # allocate
    
            R[:,0] = (data_coarse[:,3] - FU_mean)**2
            R[:,1] = (data_coarse[:,3] - FU_mean)*(data_coarse[:,4] - FV_mean)
            R[:,2] = (data_coarse[:,3] - FU_mean)*(data_coarse[:,5] - FW_mean)
            R[:,3] = R[:,1]
            R[:,4] = (data_coarse[:,4] - FV_mean)**2
            R[:,5] = (data_coarse[:,4] - FV_mean)*(data_coarse[:,5] - FW_mean)
            R[:,6] = R[:,2]
            R[:,7] = R[:,5]
            R[:,8] = (data_coarse[:,5] - FW_mean)**2
        else:
            # increment
    
            R[:,0] = R[:,0] +  (data_coarse[:,3] - FU_mean)**2
            R[:,1] = R[:,1] + ((data_coarse[:,3] - FU_mean)*(data_coarse[:,4] - FV_mean))
            R[:,2] = R[:,2] + ((data_coarse[:,3] - FU_mean)*(data_coarse[:,5] - FW_mean))
            R[:,3] = R[:,1] #Symmetric
            R[:,4] = R[:,4] +  (data_coarse[:,4] - FV_mean)**2
            R[:,5] = R[:,5] + ((data_coarse[:,4] - FV_mean)*(data_coarse[:,5] - FW_mean))
            R[:,6] = R[:,2] #Symmetric
            R[:,7] = R[:,5] #Symmetric
            R[:,8] = R[:,8] + (data_coarse[:,5] - FW_mean)**2
            
        # iteration
        k = k+1
        #print(k)
    R = R/k 
    print('Max R value:',np.amax(R),' shape:', np.shape(R))

    #assemble in a spatial data structure
    boundary = np.ones([x_size])*(-1.)

    X = np.unique(coords_coarse[:,0])
    Y = np.unique(coords_coarse[:,1])
    Z = np.unique(coords_coarse[:,2])
    I = np.size(X)
    J = np.size(Y)
    K = np.size(Z)
    print(x_size,I*J*K,I,J,K)

    # use the first xz plane with data as an inflow boundary condition
    for j in range(J-5,J):
        for i in range(0,I):
            for k in range(0,K):
                boundary[i+I*j+J*I*k] = J-j

    # use data as wall boundary condition
    for i in range(0, x_size):
        if (FU_mean[i]+FV_mean[i]+FW_mean[i]).sum() == 0:
           boundary[i] = 10.

    if do_for_single_rep:

        #just take the vel of a single repetition
        fname = flist[single_rep]
        # load data instantaneous field
        data = np.load(fname)
        data[np.isnan(data)] = 0.

        X = np.unique(data[:,0])
        Y = np.unique(data[:,1])
        Z = np.unique(data[:,2])
        I = np.size(X)
        J = np.size(Y)
        K = np.size(Z)

        velocity_coarse = downsample_data(data[:,3:6], I,J,K)
        data_coarse[:,3:6] = velocity_coarse
        
        for i in range(0, x_size):
           if boundary[i] < 10.0:
              FU_mean[i] = data_coarse[i,3]
              FV_mean[i] = data_coarse[i,4]
              FW_mean[i] = data_coarse[i,5]
 
    # check xz plane with NaN data as inflow/outflow boundary condition

    #X = np.unique(coords_coarse[:,0])
    #Y = np.unique(coords_coarse[:,1])
    #Z = np.unique(coords_coarse[:,2])
    #I = np.size(X)
    #J = np.size(Y)
    #K = np.size(Z)
 
    #for j in range(J-5,J):
    #    mean_vel_y = 0.0
    #    for i in range(0,I):
    #        for k in range(0,K):
    #            mean_vel_y = mean_vel_y + FV_mean[i+I*j+J*I*k]
    #    for i in range(0,I):
    #        for k in range(0,K):
    #            if mean_vel_y == 0.0:
    #               boundary[i+I*j+J*I*k] = 20.
    #for j in range(0,5):
    #    mean_vel_y = 0.0
    #    for i in range(0,I):
    #        for k in range(0,K):
    #            mean_vel_y = mean_vel_y + FV_mean[i+I*j+J*I*k]
    #    for i in range(0,I):
    #        for k in range(0,K):
    #            if mean_vel_y == 0.0:
    #               boundary[i+I*j+J*I*k] = 20.

    print('found boundary..',boundary.sum())

    return boundary, np.stack([FU_mean, FV_mean, FW_mean], axis=1) , R


def read_coordinates(fname):
    """Reads coordinates from a single numpy file
    """
    coords = np.load(fname)[:,:3]
    coords = coords*1e-3
    if np.isnan(coords[:,:3]).any():
        raise ValueError("NaN-values in coordinates of file %s" % fname)
    return coords[:,0], coords[:,1], coords[:,2]

def read_velocity(fname):
    """Reads velocity from a single numpy file
    """
    velocity = np.load(fname)[:,3:]
    velocity[np.isnan(velocity)] = 0.
    return velocity

def check_if_coordinates_are_identical(flist):
    """Checks if coordinates in a list of numpy files are identical
    """
    fname_ref = flist[0]
    X_coord, Y_coord, Z_coord = read_coordinates(fname_ref)
    
    # Check coordinates for equality in all files
    for fname in flist:
        data = np.load(fname)
        data = data*1e-3
        if not np.equal(np.stack([X_coord, Y_coord, Z_coord], axis=1), data[:,:3]).all():
            raise ValueError("Inconsistent coordinates between files %s and %s" % (fname_ref, fname))
        else:
            print(fname)

def read_coordinates_checked(flist):
    """First checks that list of numpy files contain identical coordinates, then returns those coordinates
    :param list of numpy files containing experimental measurements
    :returns uni-dimensional non-deduplicated x-/y-/z-coordinate arrays (as read from file)
    """
    check_if_coordinates_are_identical(flist)
    fname = flist[0]
    return read_coordinates(fname)

def get_geometry_and_spatial_conversion(X_coord, Y_coord, Z_coord):
    """Get scalar-/vector-/matrix-field array transformations uni-dimensional non-deduplicated x-/y-/z-coordinate arrays
    :param X_coord: X coordinates
    :param Y_coord: Y coordinates
    :param Z_coord: Z coordinates
    :return: conversion functions for scalar, vector, matrix arrays and coordinates accessor functions
    """
    # Compute uniaxial coordinate values
    x, x_ids = np.unique(X_coord, return_inverse=True)
    y, y_ids = np.unique(Y_coord, return_inverse=True)
    z, z_ids = np.unique(Z_coord, return_inverse=True)
    
    def convert_scalar(coord):
        return np.flip(coord.reshape(z.size, y.size, x.size), axis=1).transpose(2,1,0)
    def convert_vector(vect):
        dim = vect.size//(z.size*y.size*x.size)
        return np.flip(vect.reshape(z.size, y.size, x.size, dim), axis=1).transpose(2,1,0,3)
    def convert_matrix(matr):
        dim = int(np.sqrt(matr.size/(z.size*y.size*x.size)))
        return np.flip(matr.reshape(z.size, y.size, x.size, dim, dim), axis=1).transpose(2,1,0,3,4)

    # Print coordinate ids to understand layout
    def log_coordinate_ids(x_id, y_id, z_id):
        logging.debug("  X-indices: ")
        logging.debug(x_id[:4,:,:])
        logging.debug("  ...")
        logging.debug(x_id[-2:,:,:])
        logging.debug("  Y-indices: ")
        logging.debug(y_id[:4,:,:])
        logging.debug("  ...")
        logging.debug(y_id[-2:,:,:])
        logging.debug("  Z-indices: ")
        logging.debug(z_id[:4,:,:])
        logging.debug("  ...")
        logging.debug(z_id[-2:,:,:])
    
    # transform coordinate indices
    x_ids_corr = convert_scalar(x_ids)
    y_ids_corr = convert_scalar(y_ids)
    z_ids_corr = convert_scalar(z_ids)
    
    logging.debug("\n\n********** Coordinate ids in data files **********\n\n")
    log_coordinate_ids(x_ids.reshape(z.size, y.size, x.size),
                       y_ids.reshape(z.size, y.size, x.size),
                       z_ids.reshape(z.size, y.size, x.size))

    logging.debug("\n\n********** Coordinate ids after correction **********\n\n")
    log_coordinate_ids(x_ids_corr, y_ids_corr, z_ids_corr)
    
    # End of debugging output

    if not (np.array_equal(x_ids_corr[:-1,:,:]+1, x_ids_corr[1:,:,:]) and np.array_equal(x_ids_corr[:,:-1,:-1], x_ids_corr[:,1:,1:])):
        raise AssertionError("X coordinates not sorted in ascending order along 0th axis")
    if not (np.array_equal(y_ids_corr[:,:-1,:]+1, y_ids_corr[:,1:,:]) and np.array_equal(y_ids_corr[:-1,:,:-1], y_ids_corr[1:,:,1:])):
        raise AssertionError("Y coordinates not sorted in ascending order along 1st axis")
    if not (np.array_equal(z_ids_corr[:,:,:-1]+1, z_ids_corr[:,:,1:]) and np.array_equal(z_ids_corr[:-1,:-1,:], z_ids_corr[1:,1:,:])):
        raise AssertionError("X coordinates not sorted in ascending order along 2nd axis")
    
    return [x, y, z], convert_scalar, convert_vector, convert_matrix

# convert to time slice
def convert_spatial_to_time_slice(voxel_feature):
    sh = voxel_feature.shape
    return voxel_feature.reshape(*sh[:3], 1, *sh[3:])

# downsample, I,J,K are number of gridpoints in each direction
def downsample_data(orig_data, I, J, K):
    if I*J*K != np.shape(orig_data)[0]:
        raise ValueError("wrong number of gridpoints: ", I*J*K, np.shape(orig_data[:,0]) )
    #elif dw == 1:
    #    return orig_data    

    L = I//dw
    M = (J-8)//dw
    N = K//dw
    coords = np.empty([L*M*N,3])
    for n in range(0, N):
        for m in range(0, M):
            for l in range(0, L):
                i = l*dw
                j = m*dw + 4
                k = n*dw
                coords[l+L*m+L*M*n,0] = orig_data[i+I*j+J*I*k,0]
                coords[l+L*m+L*M*n,1] = orig_data[i+I*j+J*I*k,1]
                coords[l+L*m+L*M*n,2] = orig_data[i+I*j+J*I*k,2]

    return coords

# Read coordinates from ALL files involved (checking that they are all identical)
#X_coord, Y_coord, Z_coord = read_coordinates_checked(functools.reduce(lambda x, y: x+y, [time_slice["files"] for time_slice in time_slices], []))

X_coord, Y_coord, Z_coord = read_coordinates_checked(functools.reduce(lambda x, y: x+y, [time_slices[times] for times in time_str], []))

# if downsampling is required
orig_coords = np.load(time_slices["0.05"][0])[:,:3]

X = np.unique(orig_coords[:,0])
Y = np.unique(orig_coords[:,1])
Z = np.unique(orig_coords[:,2])
I = np.size(X)
J = np.size(Y)
K = np.size(Z)

orig_coords[np.isnan(orig_coords)] = 0.
L = I//dw
M = (J-8)//dw  
N = K//dw 
X_coord = np.empty([L*M*N,1])
Y_coord = np.empty([L*M*N,1])
Z_coord = np.empty([L*M*N,1])
downsampled_data = downsample_data(orig_coords, I,J,K)
X_coord = downsampled_data[:,0]*1e-3     
Y_coord = downsampled_data[:,1]*1e-3
Z_coord = downsampled_data[:,2]*1e-3
print('shape of downsampled data: ',np.shape(downsampled_data))


# Compute conversion of file velocities/covariances to spatial layout used in hpc-predict-io without time dimension
geometry, convert_scalar_to_spatial, convert_vector_to_spatial, convert_matrix_to_spatial = get_geometry_and_spatial_conversion(X_coord, Y_coord, Z_coord)



# if False: # Examples of applying convert_scalar/convert_vector to scalar/vectorial quantities based  (can be left out without side-effect)
#     # This transformation  has to be applied to every scalar quantity based on the same grid (such as uniaxial coordinate)
#     X_coord_corr = convert_scalar_to_spatial(X_coord)
#
#     coords = np.stack([X_coord, Y_coord, Z_coord], axis=1)
#     velocity = read_velocity(flist[0])
#
#     # This transformation  has to be applied to every vectorial quantity based on the same grid (such as velocity field)
#     coords_corr = convert_vector_to_spatial(coords)
#     velocity_corr = convert_vector_to_spatial(velocity)

time = []
velocity_mean = []
velocity_cov = []
intensity = []
segmentation_prob = []
# Iteratively collect velocity measurements over each time slice
for times in time_str:
    print('reading files at time: ',times)
    time.append(float(times))
    segmentation_prob_from_file, vel_mean_from_file, vel_cov_from_file = read_velocity_time_slice(time_slices[times])
    segmentation_prob.append(convert_scalar_to_spatial(segmentation_prob_from_file))
    velocity_mean.append(convert_vector_to_spatial(vel_mean_from_file))
    velocity_cov.append(convert_matrix_to_spatial(vel_cov_from_file))

# Merge all time slices
time = np.array(time)
segmentation_prob = np.stack(segmentation_prob, axis=3)
velocity_mean = np.stack(velocity_mean, axis=3)
velocity_cov = np.stack(velocity_cov, axis=3)
print(np.shape(intensity))
# Write flow MRI
segmented_flow_mri = SegmentedFlowMRI(geometry=geometry,
                   time=time,
		   time_heart_cycle_period=heart_cycle_period,
                   intensity=np.zeros(velocity_mean.shape[:-1]),
                   segmentation_prob=segmentation_prob,
                   velocity_mean=velocity_mean,
                   velocity_cov=velocity_cov)
segmented_flow_mri.write_hdf5(output_filename)
