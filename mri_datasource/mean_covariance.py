import os
import h5py
import numpy as np
import glob
import argparse
import logging
from mr_io import HPCPredictMRI # Requires adding ../python to PYTHONPATH

# Parse data input and output directories
def parse_args():
    #input_path  = "./bern_data_experiments_source/"
    #output_path = "./bern_data_experiments_hpc_predict/"
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate hpc-predict-io HDF5-message from preprocessed HDF5-files (Bernese experimental dataset).')
    parser.add_argument('--input', type=str, default="D:/Documents/Bern_Internship/bulk/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/data/",
                    help='Directory containing experimental data from Bern (numpy files with coordinates/velocity)')
    parser.add_argument('--output', type=str,  default="D:/HDF5_Results/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/",
                    help='Output directory for HDF5 files')
    parser.add_argument('--log', type=str, default="warn", help="Logging level")
    return parser.parse_args()



args = parse_args()
logging.basicConfig(level=args.log.upper())

# file list is the list of the n files representing the n repetitions of a specific phase. 
flist = glob.glob(args.input + '/*masked.npy')

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

output_dir = os.path.realpath(args.output)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

with h5py.File(args.output + '/Velocity_Mean.h5','w') as hf1:
    hf1.create_dataset('Velocity_mean', data=[FU_mean, FV_mean, FW_mean])

with h5py.File(args.output + '/Covariance.h5','w') as hf2:
    hf2.create_dataset('Covariance', data=[R[:,0],R[:,4],R[:,8],R[:,1],R[:,2],R[:,5]])

with h5py.File(args.output + '/Coordinates.h5','w') as hf3:
    hf3.create_dataset('Coordinates', data=[X_coord, Y_coord, Z_coord])


    
print('X_min: ',np.amin(X_coord))
print('Y_min: ',np.amin(Y_coord))
print('Z_min: ',np.amin(Z_coord))

print('X_max: ',np.amax(X_coord))
print('Y_max: ',np.amax(Y_coord))
print('Z_max: ',np.amax(Z_coord))


# Added parts

def read_coordinates(fname):
    coords = np.load(fname)[:,:3]
    if np.isnan(coords[:,:3]).any():
        raise ValueError("NaN-values in coordinates of file %s" % fname)
    return coords[:,0], coords[:,1], coords[:,2]

def read_velocity(fname):
    velocity = np.load(fname)[:,3:]
    velocity[np.isnan(velocity)] = 0.
    return velocity

def check_if_coordinates_are_identical(flist):
    fname_ref = flist[0]
    X_coord, Y_coord, Z_coord = read_coordinates(fname_ref)
    
    # Check coordinates for equality in all files
    for fname in flist:
        data = np.load(fname)
        if not np.equal(np.stack([X_coord, Y_coord, Z_coord], axis=1), data[:,:3]).all():
            raise ValueError("Inconsistent coordinates between files %s and %s" % (fname_ref, fname))

def get_spatial_conversion(X_coord, Y_coord, Z_coord):
    # Compute uniaxial coordinate values
    x, x_ids = np.unique(X_coord, return_inverse=True)
    y, y_ids = np.unique(Y_coord, return_inverse=True)
    z, z_ids = np.unique(Z_coord, return_inverse=True)
    
    def get_coordinates():
        return x, y, z
    def convert_scalar(coord):
        return np.flip(coord.reshape(z.size, y.size, x.size), axis=1).transpose(2,1,0)
    def convert_vector(vect):
        return np.flip(vect.reshape(z.size, y.size, x.size,-1), axis=1).transpose(2,1,0,3)
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

    return convert_scalar, convert_vector, convert_matrix, get_coordinates

check_if_coordinates_are_identical(flist)
X_coord, Y_coord, Z_coord = read_coordinates(flist[0])
convert_scalar_to_spatial, convert_vector_to_spatial, convert_matrix_to_spatial, get_coordinates = get_spatial_conversion(X_coord, Y_coord, Z_coord)

if False: # Examples of applying convert_scalar/convert_vector to scalar/vectorial quantities based  (can be left out without side-effect)
    # This transformation  has to be applied to every scalar quantity based on the same grid (such as uniaxial coordinate)
    X_coord_corr = convert_scalar_to_spatial(X_coord)
    
    coords = np.stack([X_coord, Y_coord, Z_coord], axis=1)
    velocity = read_velocity(flist[0])
    
    # This transformation  has to be applied to every vectorial quantity based on the same grid (such as velocity field)
    coords_corr = convert_vector_to_spatial(coords)
    velocity_corr = convert_vector_to_spatial(velocity)

geometry = get_coordinates()
time = np.array([0.])
velocity_mean = convert_vector_to_spatial(np.stack([FU_mean, FV_mean, FW_mean], axis=1))
velocity_cov = convert_matrix_to_spatial(R)
 
# convert to time slice
def convert_spatial_to_time_slice(voxel_feature):
    sh = voxel_feature.shape
    return voxel_feature.reshape(*sh[:3], 1, *sh[3:])
         
hpc_predict_mri = HPCPredictMRI(geometry=geometry,
                                time=time,
                                intensity=np.zeros((*velocity_mean.shape[:-1],1)),
                                velocity_mean=convert_spatial_to_time_slice(velocity_mean),
                                velocity_cov=convert_spatial_to_time_slice(velocity_cov))
hpc_predict_mri.write_hdf5(args.output + '/bern_experimental_dataset_hpc_predict_mri.h5')
