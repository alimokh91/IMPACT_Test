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
      # here comes the actual serialization code (transposition to use Fortran memory layout)
      voxel_feature_transposed = self.voxel_feature.transpose()
      grp = f.create_group(SpatialMRI.group_name)
      ds = grp.create_dataset("voxel_feature", voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)

  def read_hdf5(path):
    with h5py.File(path, "r") as f:
      # here comes the actual deserialization code
      return SpatialMRI(voxel_feature=f[SpaceTimeMRI.group_name]["voxel_feature"][()].transpose())


class SpaceTimeMRI:
  group_name = "space-time-mri"

  def __init__(self, geometry, voxel_feature):
    """Voxel-based parameters must be specified in (x,y,z,t,i)-order, Fortran will treat it in (i,t,x,y,z)-order.
       The index i is used as the component index (i.e. between 0..2 for mean and 0..5 for covariance of velocity field)
    """
    # a numpy array
    self.geometry = geometry
    self.voxel_feature = voxel_feature

  def write_hdf5(self, path):
    # file handling
    if os.path.isfile(path):
      raise FileExistsError("Tried to open file %s, which exists already" % path)

    with h5py.File(path, "w") as f:
      # here comes the actual serialization code (transposition to use Fortran memory layout)
      voxel_feature_transposed = self.voxel_feature.transpose((2,1,0,3,4))
      grp = f.create_group(SpaceTimeMRI.group_name)
      for i, coord_name in enumerate(["x_coordinates", "y_coordinates", "z_coordinates"]):
          grp.create_dataset(coord_name, self.geometry[i].shape, data=self.geometry[i], dtype=self.geometry[i].dtype)
      ds = grp.create_dataset("voxel_feature", voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)

  def read_hdf5(path):
    with h5py.File(path, "r") as f:
      # here comes the actual deserialization code
      return SpaceTimeMRI(geometry=[f[SpaceTimeMRI.group_name][coord_name][()] \
                                    for coord_name in ["x_coordinates", "y_coordinates", "z_coordinates"]],
                          voxel_feature=f[SpaceTimeMRI.group_name]["voxel_feature"][()].transpose((2,1,0,3,4)))


