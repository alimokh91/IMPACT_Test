import os
import h5py
import numpy as np
import glob
import argparse
import logging
import functools
from mr_io import FlowMRI # Requires adding ../python to PYTHONPATH

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
if os.path.exists(output_filename):
    raise RuntimeError("The file {} exists already. Exiting...".format(output_filename))
if not os.path.exists(args.output):
    os.makedirs(args.output)

# Define time slices (TODO: should be read from a metainformation file e.g. in JSON)
time_slices = [{"time": 0., "files": glob.glob(args.input + '/*masked.npy')}]
# time_slices = [{"time": 0.05, "files": glob.glob(args.input + '/*masked.npy')},
#                {"time": 0.06, "files": glob.glob(args.input + '/*masked.npy')}]


def read_velocity_time_slice(flist):
    # file list is the list of the n files representing the n repetitions of a specific phase.

    # calculate mean
    #print('\ncalculate mean velocity field\n')
    k = 0;
    x_size = 0
    j=0
    for fname in flist:

        # load data instantaneous field
        data = np.load(fname)
        data[np.isnan(data)] = 0.

        #For defining size of covariance matrix
        A = np.array(data[:,0])
        x_size = A.size

        if k == 0:
            # allocate
            FU_mean = data[:,3].copy()
            FV_mean = data[:,4].copy()
            FW_mean = data[:,5].copy()

        else:
            # increment
            FU_mean = FU_mean + data[:,3]
            FV_mean = FV_mean + data[:,4]
            FW_mean = FW_mean + data[:,5]

        # iterations
        k = k+1
        #print(k)

    FU_mean = FU_mean / k
    FV_mean = FV_mean / k
    FW_mean = FW_mean / k

    # calculate rms and variance and covariance
    #print('\ncalculeting rms\n')
    k = 0
    R = np.zeros((x_size,9))
    X_coord = np.zeros((x_size,1))
    Y_coord = np.zeros((x_size,1))
    Z_coord = np.zeros((x_size,1))

    for fname in flist:

        # load data
        data = np.load(fname)
        data[np.isnan(data)] = 0.

        #To get X, Y and Z coordinates
        X_coord = np.array(data[:,0])
        Y_coord = np.array(data[:,1])
        Z_coord = np.array(data[:,2])

        if k == 0:
            # allocate

            R[:,0] = (data[:,3] - FU_mean)**2
            R[:,1] = (data[:,3] - FU_mean)*(data[:,4] - FV_mean)
            R[:,2] = (data[:,3] - FU_mean)*(data[:,5] - FW_mean)
            R[:,3] = R[:,1]
            R[:,4] = (data[:,4] - FV_mean)**2
            R[:,5] = (data[:,4] - FV_mean)*(data[:,5] - FW_mean)
            R[:,6] = R[:,2]
            R[:,7] = R[:,5]
            R[:,8] = (data[:,5] - FW_mean)**2

        else:
            # increment

            R[:,0] = R[:,0] + (data[:,3] - FU_mean)**2
            R[:,1] = R[:,1] + ((data[:,3] - FU_mean)*(data[:,4] - FV_mean))
            R[:,2] = R[:,2] + ((data[:,3] - FU_mean)*(data[:,5] - FW_mean))
            R[:,3] = R[:,1] #Symmetric
            R[:,4] = R[:,4] + (data[:,4] - FV_mean)**2
            R[:,5] = R[:,5] + ((data[:,4] - FV_mean)*(data[:,5] - FW_mean))
            R[:,6] = R[:,2] #Symmetric
            R[:,7] = R[:,5] #Symmetric
            R[:,8] = R[:,8] + (data[:,5] - FW_mean)**2

        # iteration
        k = k+1
        #print(k)

    R = R/k

    return np.stack([FU_mean, FV_mean, FW_mean], axis=1) , R

    # output_dir = os.path.realpath(args.output)
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    #
    # with h5py.File(args.output + '/Velocity_Mean.h5','w') as hf1:
    #     hf1.create_dataset('Velocity_mean', data=[FU_mean, FV_mean, FW_mean])
    #
    # with h5py.File(args.output + '/Covariance.h5','w') as hf2:
    #     hf2.create_dataset('Covariance', data=[R[:,0],R[:,4],R[:,8],R[:,1],R[:,2],R[:,5]])
    #
    # with h5py.File(args.output + '/Coordinates.h5','w') as hf3:
    #     hf3.create_dataset('Coordinates', data=[X_coord, Y_coord, Z_coord])
    #
    #
    #
    # print('X_min: ',np.amin(X_coord))
    # print('Y_min: ',np.amin(Y_coord))
    # print('Z_min: ',np.amin(Z_coord))
    #
    # print('X_max: ',np.amax(X_coord))
    # print('Y_max: ',np.amax(Y_coord))
    # print('Z_max: ',np.amax(Z_coord))


def read_coordinates(fname):
    """Reads coordinates from a single numpy file
    """
    coords = np.load(fname)[:,:3]
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
        if not np.equal(np.stack([X_coord, Y_coord, Z_coord], axis=1), data[:,:3]).all():
            raise ValueError("Inconsistent coordinates between files %s and %s" % (fname_ref, fname))

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


# Read coordinates from ALL files involved (checking that they are all identical)
X_coord, Y_coord, Z_coord = read_coordinates_checked(functools.reduce(lambda x, y: x+y, [time_slice["files"] for time_slice in time_slices], []))
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

# Iteratively collect velocity measurements over each time slice
for time_slice in time_slices:
    time.append(time_slice["time"])
    vel_mean_from_file, vel_cov_from_file = read_velocity_time_slice(time_slice["files"])
    velocity_mean.append(convert_vector_to_spatial(vel_mean_from_file))
    velocity_cov.append(convert_matrix_to_spatial(vel_cov_from_file))

# Merge all time slices
time = np.array(time)
velocity_mean = np.stack(velocity_mean, axis=3)
velocity_cov = np.stack(velocity_cov, axis=3)

# Write flow MRI
flow_mri = FlowMRI(geometry=geometry,
                   time=time,
                   intensity=np.zeros(velocity_mean.shape[:-1]),
                   velocity_mean=velocity_mean,
                   velocity_cov=velocity_cov)
flow_mri.write_hdf5(output_filename)
