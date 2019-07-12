import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, HPCPredictMRI, SegmentedHPCPredictMRI

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
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        voxel_feature = np.random.rand(23,13,19,17,21)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
 
    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_space_time.h5")
 
        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_space_time","mr_io_test_space_time.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
 
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")   
        print(fort_stderr)
 
        fort_group_name, fort_x_dim_str, fort_x_coord_str,fort_y_dim_str, fort_y_coord_str, fort_z_dim_str, fort_z_coord_str, fort_t_dim_str, fort_t_coord_str, fort_dims_str, fort_array_str, _ = fort_stdout.split('\n')
 
        fort_x_dim = np.fromstring(fort_x_dim_str, dtype=int, sep=' ')
        fort_y_dim = np.fromstring(fort_y_dim_str, dtype=int, sep=' ')
        fort_z_dim = np.fromstring(fort_z_dim_str, dtype=int, sep=' ')
        fort_t_dim = np.fromstring(fort_t_dim_str, dtype=int, sep=' ')
 
        fort_x_coord = np.fromstring(fort_x_coord_str, dtype=float, sep=' ')
        fort_y_coord = np.fromstring(fort_y_coord_str, dtype=float, sep=' ')
        fort_z_coord = np.fromstring(fort_z_coord_str, dtype=float, sep=' ')
        fort_t_coord = np.fromstring(fort_t_coord_str, dtype=float, sep=' ')
      
        #print("Fortran x-coordinates: "); print(fort_x_coord)        
        #print("Python x-coordinates:  "); print(self.mri.geometry[0])
        #print("Fortran y-coordinates: "); print(fort_y_coord)   
        #print("Python y-coordinates:  "); print(self.mri.geometry[1])
        #print("Fortran z-coordinates: "); print(fort_z_coord)   
        #print("Python z-coordinates:  "); print(self.mri.geometry[2]) 
 
        self.assertEqual(fort_x_dim, fort_x_coord.shape[0])
        self.assertEqual(fort_y_dim, fort_y_coord.shape[0])
        self.assertEqual(fort_z_dim, fort_z_coord.shape[0])
        self.assertEqual(fort_t_dim, fort_t_coord.shape[0])
  
        self.assertTrue(np.allclose(fort_x_coord, self.mri.geometry[0], rtol=1e-14))
        self.assertTrue(np.allclose(fort_y_coord, self.mri.geometry[1], rtol=1e-14))
        self.assertTrue(np.allclose(fort_z_coord, self.mri.geometry[2], rtol=1e-14))
        self.assertTrue(np.allclose(fort_t_coord, self.mri.time, rtol=1e-14))
 
        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
 
        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose((2,1,0,3,4))
        # self.assertTrue(np.array_equal(fort_dims,self.mri.voxel_feature.shape))
 
        self.assertTrue(np.allclose(fort_array, self.mri.voxel_feature, rtol=1e-14))
 
    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_space_time.h5")
 
 
class TestHPCPredictMRI(unittest.TestCase):
 
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity = np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        self.mri = HPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov)
 
    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_hpc_predict.h5")
 
        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_hpc_predict","mr_io_test_hpc_predict.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
 
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")   
        print(fort_stderr)
 
        fort_group_name, fort_x_dim_str, fort_x_coord_str,fort_y_dim_str, fort_y_coord_str, fort_z_dim_str, fort_z_coord_str, fort_t_dim_str, fort_t_coord_str,  fort_dims_str0, fort_array_str0, fort_dims_str, fort_array_str, fort_dims_str2, fort_array_str2, _ = fort_stdout.split('\n')
 
        fort_x_dim = np.fromstring(fort_x_dim_str, dtype=int, sep=' ')
        fort_y_dim = np.fromstring(fort_y_dim_str, dtype=int, sep=' ')
        fort_z_dim = np.fromstring(fort_z_dim_str, dtype=int, sep=' ')
        fort_t_dim = np.fromstring(fort_t_dim_str, dtype=int, sep=' ')
 
        fort_x_coord = np.fromstring(fort_x_coord_str, dtype=float, sep=' ')
        fort_y_coord = np.fromstring(fort_y_coord_str, dtype=float, sep=' ')
        fort_z_coord = np.fromstring(fort_z_coord_str, dtype=float, sep=' ')
        fort_t_coord = np.fromstring(fort_t_coord_str, dtype=float, sep=' ')
      
        #print("Fortran x-coordinates: "); print(fort_x_coord)        
        #print("Python x-coordinates:  "); print(self.mri.geometry[0])
        #print("Fortran y-coordinates: "); print(fort_y_coord)   
        #print("Python y-coordinates:  "); print(self.mri.geometry[1])
        #print("Fortran z-coordinates: "); print(fort_z_coord)   
        #print("Python z-coordinates:  "); print(self.mri.geometry[2]) 
 
        self.assertEqual(fort_x_dim, fort_x_coord.shape[0])
        self.assertEqual(fort_y_dim, fort_y_coord.shape[0])
        self.assertEqual(fort_z_dim, fort_z_coord.shape[0])
        self.assertEqual(fort_t_dim, fort_t_coord.shape[0])
  
        self.assertTrue(np.allclose(fort_x_coord, self.mri.geometry[0], rtol=1e-14))
        self.assertTrue(np.allclose(fort_y_coord, self.mri.geometry[1], rtol=1e-14))
        self.assertTrue(np.allclose(fort_z_coord, self.mri.geometry[2], rtol=1e-14))
        self.assertTrue(np.allclose(fort_t_coord, self.mri.time, rtol=1e-14))
 
        fort_dims0 = np.fromstring(fort_dims_str0, dtype=int, sep=' ')
        fort_array = np.fromstring(fort_array_str0, dtype=float, sep=' ').reshape(np.flip(fort_dims0)).transpose((2,1,0,3))
        self.assertTrue(np.allclose(fort_array, self.mri.intensity, rtol=1e-14))
 
        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose((2,1,0,3,4))
        self.assertTrue(np.allclose(fort_array, self.mri.velocity_mean, rtol=1e-14))
 
        fort_dims2 = np.fromstring(fort_dims_str2, dtype=int, sep=' ')
        fort_array2 = np.fromstring(fort_array_str2, dtype=float, sep=' ').reshape(np.flip(fort_dims2)).transpose((2,1,0,3,5,4))
        self.assertTrue(np.allclose(fort_array2, self.mri.velocity_cov, rtol=1e-14))
 
    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_hpc_predict.h5")


class TestSegmentedHPCPredictMRI(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity = np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        segmentation_prob = np.random.rand(23,13,19,17)
        self.mri = SegmentedHPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov, segmentation_prob)

    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_segmented_hpc_predict.h5")

        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_segmented_hpc_predict","mr_io_test_segmented_hpc_predict.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")   
        print(fort_stderr)

        fort_group_name, fort_x_dim_str, fort_x_coord_str,fort_y_dim_str, fort_y_coord_str, fort_z_dim_str, fort_z_coord_str, fort_t_dim_str, fort_t_coord_str,  fort_dims_str0, fort_array_str0, fort_dims_str, fort_array_str, fort_dims_str2, fort_array_str2, fort_dims_str3, fort_array_str3, _ = fort_stdout.split('\n')

        fort_x_dim = np.fromstring(fort_x_dim_str, dtype=int, sep=' ')
        fort_y_dim = np.fromstring(fort_y_dim_str, dtype=int, sep=' ')
        fort_z_dim = np.fromstring(fort_z_dim_str, dtype=int, sep=' ')
        fort_t_dim = np.fromstring(fort_t_dim_str, dtype=int, sep=' ')

        fort_x_coord = np.fromstring(fort_x_coord_str, dtype=float, sep=' ')
        fort_y_coord = np.fromstring(fort_y_coord_str, dtype=float, sep=' ')
        fort_z_coord = np.fromstring(fort_z_coord_str, dtype=float, sep=' ')
        fort_t_coord = np.fromstring(fort_t_coord_str, dtype=float, sep=' ')
     
        #print("Fortran x-coordinates: "); print(fort_x_coord)        
        #print("Python x-coordinates:  "); print(self.mri.geometry[0])
        #print("Fortran y-coordinates: "); print(fort_y_coord)   
        #print("Python y-coordinates:  "); print(self.mri.geometry[1])
        #print("Fortran z-coordinates: "); print(fort_z_coord)   
        #print("Python z-coordinates:  "); print(self.mri.geometry[2]) 

        self.assertEqual(fort_x_dim, fort_x_coord.shape[0])
        self.assertEqual(fort_y_dim, fort_y_coord.shape[0])
        self.assertEqual(fort_z_dim, fort_z_coord.shape[0])
        self.assertEqual(fort_t_dim, fort_t_coord.shape[0])
 
        self.assertTrue(np.allclose(fort_x_coord, self.mri.geometry[0], rtol=1e-14))
        self.assertTrue(np.allclose(fort_y_coord, self.mri.geometry[1], rtol=1e-14))
        self.assertTrue(np.allclose(fort_z_coord, self.mri.geometry[2], rtol=1e-14))
        self.assertTrue(np.allclose(fort_t_coord, self.mri.time, rtol=1e-14))

        fort_dims0 = np.fromstring(fort_dims_str0, dtype=int, sep=' ')
        fort_array0 = np.fromstring(fort_array_str0, dtype=float, sep=' ').reshape(np.flip(fort_dims0)).transpose((2,1,0,3))
        self.assertTrue(np.allclose(fort_array0, self.mri.intensity, rtol=1e-14))

        fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
        fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose((2,1,0,3,4))
        self.assertTrue(np.allclose(fort_array, self.mri.velocity_mean, rtol=1e-14))

        fort_dims2 = np.fromstring(fort_dims_str2, dtype=int, sep=' ')
        fort_array2 = np.fromstring(fort_array_str2, dtype=float, sep=' ').reshape(np.flip(fort_dims2)).transpose((2,1,0,3,5,4))
        self.assertTrue(np.allclose(fort_array2, self.mri.velocity_cov, rtol=1e-14))

        fort_dims3 = np.fromstring(fort_dims_str3, dtype=int, sep=' ')
        fort_array3 = np.fromstring(fort_array_str3, dtype=float, sep=' ').reshape(np.flip(fort_dims3)).transpose((2,1,0,3))
        self.assertTrue(np.allclose(fort_array3, self.mri.segmentation_prob, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_segmented_hpc_predict.h5")


if __name__ == '__main__':
    unittest.main()

