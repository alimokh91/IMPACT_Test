import numpy as np
from mr_io import MRI

voxel_feature = np.random.rand(3,2)
mri = MRI(voxel_feature)
print(mri)
mri.write_hdf5("mr_io_test.h5")
