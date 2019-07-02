import h5py

import numpy as np
from numpy import arctan2, sqrt, pi

import glob

# file list is the list of the n files representing the n repetitions of a specific phase. 
flist = glob.glob('D:/Documents/Bern_Internship/bulk/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/data/*masked.npy')

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

with h5py.File('D:/HDF5_Results/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/Velocity_Mean.h5','w') as hf1:
    hf1.create_dataset('Velocity_mean', data=[FU_mean, FV_mean, FW_mean])

with h5py.File('D:/HDF5_Results/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/Covariance.h5','w') as hf2:
    hf2.create_dataset('Covariance', data=[R[:,0],R[:,4],R[:,8],R[:,1],R[:,2],R[:,5]])

with h5py.File('D:/HDF5_Results/Cam_Date=180306_Time=163841_ImgPreproc_FastMART_TomoPIV_DirectCorrelation_48x48x48_75ov=unknown/Coordinates.h5','w') as hf3:
    hf3.create_dataset('Coordinates', data=[X_coord, Y_coord, Z_coord])

    
print('X_min: ',np.amin(X_coord))
print('Y_min: ',np.amin(Y_coord))
print('Z_min: ',np.amin(Z_coord))

print('X_max: ',np.amax(X_coord))
print('Y_max: ',np.amax(Y_coord))
print('Z_max: ',np.amax(Z_coord))

'''