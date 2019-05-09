import h5py
import os

class MRI:
  def __init__(self, voxel_feature):
    # an example numpy array
    self.voxel_feature = voxel_feature
  
  def write_hdf5(self, path):
    # file handling
    if os.path.isfile(path):
      raise FileExistsError("Tried to open file %s, which exists already" % path)
    
    with h5py.File(path, "w") as f:
      # here comes the actual serialization code
      f.create_dataset("voxel_feature", self.voxel_feature.shape, dtype=self.voxel_feature.dtype)

def read_hdf5(path):
  with h5py.File(path, "r") as f:
    # here comes the actual deserialization code
    return MRI(voxel_feature=f["voxel_feature"][()])

