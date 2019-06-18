import h5py
import os

class MRI:
  group_name = "3d-mri"  

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
      grp = f.create_group(MRI.group_name)
      ds = grp.create_dataset("voxel_feature", voxel_feature_transposed.shape, data=voxel_feature_transposed, dtype=voxel_feature_transposed.dtype)

  def read_hdf5(path):
    with h5py.File(path, "r") as f:
      # here comes the actual deserialization code
      return MRI(voxel_feature=f[MRI.group_name]["voxel_feature"][()].transpose())



