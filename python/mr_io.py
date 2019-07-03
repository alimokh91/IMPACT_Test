import h5py
import os

class SpatialMRI:
    group_name = "spatial-mri"  

    def __init__(self, voxel_feature):
        # a numpy array
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

def write_space_time_coordinates(grp, geometry, time):
    for i, coord_name in enumerate(["x_coordinates", "y_coordinates", "z_coordinates"]):
        grp.create_dataset(coord_name, geometry[i].shape, data=geometry[i], dtype=geometry[i].dtype)
    grp.create_dataset("t_coordinates", time.shape, data=time, dtype=time.dtype)
    
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
        # a numpy array
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
    
    def __init__(self, geometry, time, velocity_mean, velocity_cov):
        """Voxel-based parameters must be specified in (x,y,z,t,i)-order, Fortran will treat it in (i,t,x,y,z)-order.
           The index i is used as the component index (i.e. between 0..2 for mean and 0..5 for covariance of velocity field)
        """
        # a numpy array
        self.geometry = geometry
        self.time = time
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
            write_space_time_voxel_feature(grp, "velocity_mean", self.velocity_mean)
            write_space_time_voxel_matrix_feature(grp, "velocity_cov", self.velocity_cov)

    def read_hdf5(path):
        with h5py.File(path, "r") as f:
            # here comes the actual deserialization code
            return HPCPredictMRI(geometry=[f[HPCPredictMRI.group_name][coord_name][()] \
                                          for coord_name in ["x_coordinates", "y_coordinates", "z_coordinates"]],
                                time=f[HPCPredictMRI.group_name]["t_coordinates"][()],
                                velocity_mean=f[HPCPredictMRI.group_name]["velocity_mean"][()].transpose((2,1,0,3,4)),
                                velocity_cov=f[HPCPredictMRI.group_name]["velocity_cov"][()].transpose((2,1,0,3,5,4)))
