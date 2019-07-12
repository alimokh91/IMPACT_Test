import h5py
import os
import logging

class SpatialMRI:
    group_name = "spatial-mri"  

    def __init__(self, voxel_feature):
        # a numpy array
        if not(len(voxel_feature.shape) == 3):
            raise ValueError("voxel_feature must be a 3D-array (instead it has shape {}).".format(voxel_feature.shape))
        self.voxel_feature = voxel_feature
        

    def write_hdf5(self, path):
        # file handling
        if os.path.isfile(path):
            raise FileExistsError("Tried to open file %s, which exists already" % path)
    
        with h5py.File(path, "w") as f:
            voxel_feature_transposed = self.voxel_feature.transpose()
            grp = f.create_group(SpatialMRI.group_name)
            ds = grp.create_dataset("voxel_feature", voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)

    def read_hdf5(path):
        with h5py.File(path, "r") as f:
            # here comes the actual deserialization code
            return SpatialMRI(voxel_feature=f[SpatialMRI.group_name]["voxel_feature"][()].transpose())


def validate_spacetime_coordinates(geometry, time):
    for i in range(3):
        if not(len(geometry[i].shape) == 1):
            raise ValueError("geometry[{}] must be a 3D-array (instead it has shape {}).".format(i, geometry[i].shape))
    if not(len(time.shape) == 1):
        raise ValueError("time must be a 1D-array (instead it has shape {}).".format(time.shape))

def validate_spacetime_feature_coordinate_dims(geometry, time, voxel_feature):
    for i in range(3):
        if not(geometry[i].shape[0] == voxel_feature.shape[i]):
            raise ValueError("geometry[{}] and {}-th dimension in voxel_feature have inconsistent shape: {} vs. {}.".format(i, i, geometry[i].shape[0], voxel_feature.shape[i]))        
    if not(time.shape[0] == voxel_feature.shape[3]):
        raise ValueError("time and {}-th dimension in voxel_feature have inconsistent shape: {} vs. {}.".format(3, time.shape[0], voxel_feature.shape[3]))        

def validate_spacetime_scalar_feature(cls, geometry, time, voxel_feature):
    validate_spacetime_feature_coordinate_dims(geometry, time, voxel_feature)
    if not(len(voxel_feature.shape) == 4):
        raise ValueError("Spacetime scalar field must be a 4D-array (instead it has shape {}).".format(voxel_feature.shape))        
        
def validate_spacetime_vector_feature(cls, geometry, time, voxel_feature):
    validate_spacetime_feature_coordinate_dims(geometry, time, voxel_feature)
    if not(len(voxel_feature.shape) == 5):
        raise ValueError("Spacetime vector field must be a 5D-array (instead it has shape {}).".format(voxel_feature.shape))        
    if not(voxel_feature.shape[4] == 3):
        logging.warning("Constructing {} with non-3-dimensional vector field (instead {}-dimensional)".format(cls.__name__, voxel_feature.shape[4]))            
        
def validate_spacetime_matrix_feature(cls, geometry, time, voxel_feature):
    validate_spacetime_feature_coordinate_dims(geometry, time, voxel_feature)
    if not(len(voxel_feature.shape) == 6):
        raise ValueError("Spacetime matrix field must be a 6D-array (instead it has shape {}).".format(voxel_feature.shape))        
    if not(voxel_feature.shape[4] == 3) or not(voxel_feature.shape[5] == 3):
        logging.warning("Constructing {} with non-3x3-dimensional matrix field (instead {}x{}-dimensional)".format(cls.__name__, voxel_feature.shape[4], voxel_feature.shape[5]))                
        

def write_space_time_coordinates(grp, geometry, time):
    for i, coord_name in enumerate(["x_coordinates", "y_coordinates", "z_coordinates"]):
        grp.create_dataset(coord_name, geometry[i].shape, data=geometry[i], dtype=geometry[i].dtype)
    grp.create_dataset("t_coordinates", time.shape, data=time, dtype=time.dtype)
    
def write_space_time_voxel_scalar_feature(grp, name, voxel_feature):
    voxel_feature_transposed = voxel_feature.transpose((2,1,0,3))
    ds = grp.create_dataset(name, voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)
    
def write_space_time_voxel_feature(grp, name, voxel_feature):
    voxel_feature_transposed = voxel_feature.transpose((2,1,0,3,4))
    ds = grp.create_dataset(name, voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)
    
def write_space_time_voxel_matrix_feature(grp, name, voxel_feature):
    voxel_feature_transposed = voxel_feature.transpose((2,1,0,3,5,4))
    ds = grp.create_dataset(name, voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)

class SpaceTimeMRI:
    group_name = "space-time-mri"
    
    def __init__(self, geometry, time, voxel_feature):
        """Voxel-based parameters must be specified in (x,y,z,t,i)-order, Fortran will treat it in (i,t,x,y,z)-order.
           The index i is used as the component index (i.e. between 0..2 for mean and 0..5 for covariance of velocity field)
        """
        validate_spacetime_coordinates(geometry, time)
        validate_spacetime_vector_feature(SpaceTimeMRI, geometry, time, voxel_feature)
        self.geometry = geometry
        self.time = time
        self.voxel_feature = voxel_feature

    def write_hdf5(self, path):
        # file handling
        if os.path.isfile(path):
            raise FileExistsError("Tried to open file %s, which exists already" % path)
    
        with h5py.File(path, "w") as f:
            # here comes the actual serialization code (transposition to use Fortran memory layout)
            grp = f.create_group(SpaceTimeMRI.group_name)
            write_space_time_coordinates(grp, self.geometry, self.time)
            write_space_time_voxel_feature(grp, "voxel_feature", self.voxel_feature)

    def read_hdf5(path):
        with h5py.File(path, "r") as f:
            # here comes the actual deserialization code
            return SpaceTimeMRI(geometry=[f[SpaceTimeMRI.group_name][coord_name][()] \
                                          for coord_name in ["x_coordinates", "y_coordinates", "z_coordinates"]],
                                time=f[SpaceTimeMRI.group_name]["t_coordinates"][()],
                                voxel_feature=f[SpaceTimeMRI.group_name]["voxel_feature"][()].transpose((2,1,0,3,4)))


class HPCPredictMRI:
    group_name = "hpc-predict-mri"
    
    def __init__(self, geometry, time, intensity, velocity_mean, velocity_cov):
        """Voxel-based parameters must be specified in (x,y,z,t,i)-order, Fortran will treat it in (i,t,x,y,z)-order.
           The index i is used as the component index (i.e. between 0..2 for mean and 0..5 for covariance of velocity field)
        """
        validate_spacetime_coordinates(geometry, time)
        validate_spacetime_scalar_feature(HPCPredictMRI, geometry, time, intensity)
        validate_spacetime_vector_feature(HPCPredictMRI, geometry, time, velocity_mean)
        validate_spacetime_matrix_feature(HPCPredictMRI, geometry, time, velocity_cov)

        self.geometry = geometry
        self.time = time
        self.intensity= intensity
        self.velocity_mean = velocity_mean
        self.velocity_cov= velocity_cov

    def write_hdf5(self, path):
        # file handling
        if os.path.isfile(path):
            raise FileExistsError("Tried to open file %s, which exists already" % path)
    
        with h5py.File(path, "w") as f:
            # here comes the actual serialization code (transposition to use Fortran memory layout)
            grp = f.create_group(HPCPredictMRI.group_name)
            write_space_time_coordinates(grp, self.geometry, self.time)
            write_space_time_voxel_scalar_feature(grp, "intensity", self.intensity)
            write_space_time_voxel_feature(grp, "velocity_mean", self.velocity_mean)
            write_space_time_voxel_matrix_feature(grp, "velocity_cov", self.velocity_cov)

    def read_hdf5(path):
        with h5py.File(path, "r") as f:
            # here comes the actual deserialization code
            return HPCPredictMRI(geometry=[f[HPCPredictMRI.group_name][coord_name][()] \
                                          for coord_name in ["x_coordinates", "y_coordinates", "z_coordinates"]],
                                time=f[HPCPredictMRI.group_name]["t_coordinates"][()],
                                intensity=f[HPCPredictMRI.group_name]["intensity"][()].transpose((2,1,0,3)),
                                velocity_mean=f[HPCPredictMRI.group_name]["velocity_mean"][()].transpose((2,1,0,3,4)),
                                velocity_cov=f[HPCPredictMRI.group_name]["velocity_cov"][()].transpose((2,1,0,3,5,4)))


#TODO: Refactor this into class hierarchy with HPCPredictMRI
class SegmentedHPCPredictMRI:
    group_name = "segmented-hpc-predict-mri"
    
    def __init__(self, geometry, time, intensity, velocity_mean, velocity_cov, segmentation_prob):
        """Voxel-based parameters must be specified in (x,y,z,t,i)-order, Fortran will treat it in (i,t,x,y,z)-order.
           The index i is used as the component index (i.e. between 0..2 for mean and 0..5 for covariance of velocity field)
        """
        validate_spacetime_coordinates(geometry, time)
        validate_spacetime_scalar_feature(SegmentedHPCPredictMRI, geometry, time, intensity)
        validate_spacetime_vector_feature(SegmentedHPCPredictMRI, geometry, time, velocity_mean)
        validate_spacetime_matrix_feature(SegmentedHPCPredictMRI, geometry, time, velocity_cov)
        validate_spacetime_scalar_feature(SegmentedHPCPredictMRI, geometry, time, segmentation_prob)

        self.geometry = geometry
        self.time = time
        self.intensity= intensity
        self.velocity_mean = velocity_mean
        self.velocity_cov= velocity_cov
        self.segmentation_prob = segmentation_prob

    def write_hdf5(self, path):
        # file handling
        if os.path.isfile(path):
            raise FileExistsError("Tried to open file %s, which exists already" % path)
    
        with h5py.File(path, "w") as f:
            # here comes the actual serialization code (transposition to use Fortran memory layout)
            grp = f.create_group(SegmentedHPCPredictMRI.group_name)
            write_space_time_coordinates(grp, self.geometry, self.time)
            write_space_time_voxel_scalar_feature(grp, "intensity", self.intensity)
            write_space_time_voxel_feature(grp, "velocity_mean", self.velocity_mean)
            write_space_time_voxel_matrix_feature(grp, "velocity_cov", self.velocity_cov)
            write_space_time_voxel_scalar_feature(grp, "segmentation_prob", self.segmentation_prob)

    def read_hdf5(path):
        with h5py.File(path, "r") as f:
            # here comes the actual deserialization code
            return SegmentedHPCPredictMRI(geometry=[f[SegmentedHPCPredictMRI.group_name][coord_name][()] \
                                          for coord_name in ["x_coordinates", "y_coordinates", "z_coordinates"]],
                                time=f[SegmentedHPCPredictMRI.group_name]["t_coordinates"][()],
                                intensity=f[SegmentedHPCPredictMRI.group_name]["intensity"][()].transpose((2,1,0,3)),
                                velocity_mean=f[SegmentedHPCPredictMRI.group_name]["velocity_mean"][()].transpose((2,1,0,3,4)),
                                velocity_cov=f[SegmentedHPCPredictMRI.group_name]["velocity_cov"][()].transpose((2,1,0,3,5,4)),
                                segmentation_prob=f[SegmentedHPCPredictMRI.group_name]["segmentation_prob"][()].transpose((2,1,0,3)))

