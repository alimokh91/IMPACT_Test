import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI

class TestSpatialMRI(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(83,63,93)
        self.mri = SpatialMRI(voxel_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test.h5")
        
        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader","mr_io_test.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)

        fort_group_name, fort_dims_str, fort_array_str, _ = fort_stdout.split('\n')
        
        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
        
        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose()
        self.assertTrue(np.array_equal(fort_dims,self.mri.voxel_feature.shape))
        
        self.assertTrue(np.allclose(fort_array, self.mri.voxel_feature, rtol=1e-14))
        
    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test.h5") 



class TestSpaceTimeMRI(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        voxel_feature = np.random.rand(23,13,19,17,21)
        self.mri = SpaceTimeMRI(geometry, voxel_feature)

    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_space_time.h5")

        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_space_time","mr_io_test_space_time.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")   
        print(fort_stderr)

        fort_group_name, fort_x_dim_str, fort_x_coord_str,fort_y_dim_str, fort_y_coord_str, fort_z_dim_str, fort_z_coord_str, fort_dims_str, fort_array_str, _ = fort_stdout.split('\n')

        fort_x_dim = np.fromstring(fort_x_dim_str, dtype=int, sep=' ')
        fort_y_dim = np.fromstring(fort_y_dim_str, dtype=int, sep=' ')
        fort_z_dim = np.fromstring(fort_z_dim_str, dtype=int, sep=' ')

        fort_x_coord = np.fromstring(fort_x_coord_str, dtype=float, sep=' ')
        fort_y_coord = np.fromstring(fort_y_coord_str, dtype=float, sep=' ')
        fort_z_coord = np.fromstring(fort_z_coord_str, dtype=float, sep=' ')
     
        #print("Fortran x-coordinates: "); print(fort_x_coord)        
        #print("Python x-coordinates:  "); print(self.mri.geometry[0])
        #print("Fortran y-coordinates: "); print(fort_y_coord)   
        #print("Python y-coordinates:  "); print(self.mri.geometry[1])
        #print("Fortran z-coordinates: "); print(fort_z_coord)   
        #print("Python z-coordinates:  "); print(self.mri.geometry[2]) 

        self.assertEqual(fort_x_dim, fort_x_coord.shape[0])
        self.assertEqual(fort_y_dim, fort_y_coord.shape[0])
        self.assertEqual(fort_z_dim, fort_z_coord.shape[0])
 
        self.assertTrue(np.allclose(fort_x_coord, self.mri.geometry[0], rtol=1e-14))
        self.assertTrue(np.allclose(fort_y_coord, self.mri.geometry[1], rtol=1e-14))
        self.assertTrue(np.allclose(fort_z_coord, self.mri.geometry[2], rtol=1e-14))

        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')

        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose((2,1,0,3,4))
        # self.assertTrue(np.array_equal(fort_dims,self.mri.voxel_feature.shape))

        print("Fortran dims:"); print(fort_dims)
        print("Python dims: "); print(self.mri.voxel_feature.shape)
        self.assertTrue(np.allclose(fort_array, self.mri.voxel_feature, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_space_time.h5")

if __name__ == '__main__':
    unittest.main()

