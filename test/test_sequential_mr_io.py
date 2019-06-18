import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import MRI

class Test2DMatrix(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(863,639)
        self.mri = MRI(voxel_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test.h5")
        
        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader","mr_io_test.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        
        fort_group_name, fort_dims_str, fort_array_str, _ = fort_stdout.split('\n')
        
        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
        
        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose()
        self.assertTrue(np.array_equal(fort_dims,self.mri.voxel_feature.shape))
        
        self.assertTrue(np.allclose(fort_array, self.mri.voxel_feature, rtol=1e-14))
        
    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test.h5") 

if __name__ == '__main__':
    unittest.main()

