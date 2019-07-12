import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, HPCPredictMRI, SegmentedHPCPredictMRI

class TestSpatialMRIBidirectional(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(71,61,97)
        self.mri = SpatialMRI(voxel_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_in.h5")
        
        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_writer","mr_io_test_in.h5","mr_io_test_out.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        print(fort.stderr.decode("utf-8"))

        out_mri = SpatialMRI.read_hdf5("mr_io_test_out.h5")

        self.assertTrue(np.array_equal(out_mri.voxel_feature.shape, self.mri.voxel_feature.shape))        
        self.assertTrue(np.allclose(out_mri.voxel_feature, self.mri.voxel_feature, rtol=1e-14))
        
    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_in.h5")
        os.remove("mr_io_test_out.h5") 


class TestSpaceTimeMRIBidirectional(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(13)
        geometry = [np.random.rand(26), np.random.rand(11), np.random.rand(7)]
        voxel_feature = np.random.rand(26,11,7,13,19)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)

    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_space_time_in.h5")

        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_writer_space_time","mr_io_test_space_time_in.h5","mr_io_test_space_time_out.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        print(fort.stderr.decode("utf-8"))

        out_mri = SpaceTimeMRI.read_hdf5("mr_io_test_space_time_out.h5")
                
        for i in range(3):
            self.assertTrue(np.allclose(out_mri.geometry[i], self.mri.geometry[i], rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.time, self.mri.time, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.voxel_feature, self.mri.voxel_feature, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_space_time_in.h5")
        os.remove("mr_io_test_space_time_out.h5")


class TestHPCPredictMRIBidirectional(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity= np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        self.mri = HPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov)

    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_hpc_predict_in.h5")

        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_writer_hpc_predict","mr_io_test_hpc_predict_in.h5","mr_io_test_hpc_predict_out.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        print(fort.stderr.decode("utf-8"))

        out_mri = HPCPredictMRI.read_hdf5("mr_io_test_hpc_predict_out.h5")
        
        for i in range(3):
            self.assertTrue(np.allclose(out_mri.geometry[i], self.mri.geometry[i], rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.time, self.mri.time, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.intensity, self.mri.intensity, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.velocity_mean, self.mri.velocity_mean, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.velocity_cov, self.mri.velocity_cov, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_hpc_predict_in.h5")
        os.remove("mr_io_test_hpc_predict_out.h5")
        

class TestSegmentedHPCPredictMRIBidirectional(unittest.TestCase):

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity= np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        segmentation_prob = np.random.rand(23,13,19,17)
        self.mri = SegmentedHPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov, segmentation_prob)

    def test_communicator(self):
        # Write HDF5 from Python
        self.mri.write_hdf5("mr_io_test_segmented_hpc_predict_in.h5")

        # Read HDF5 from Fortran
        fort = sp.run(["fortran/mr_io_test_reader_writer_segmented_hpc_predict","mr_io_test_segmented_hpc_predict_in.h5","mr_io_test_segmented_hpc_predict_out.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        print(fort.stderr.decode("utf-8"))

        out_mri = SegmentedHPCPredictMRI.read_hdf5("mr_io_test_segmented_hpc_predict_out.h5")
        
        for i in range(3):
            self.assertTrue(np.allclose(out_mri.geometry[i], self.mri.geometry[i], rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.time, self.mri.time, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.intensity, self.mri.intensity, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.velocity_mean, self.mri.velocity_mean, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.velocity_cov, self.mri.velocity_cov, rtol=1e-14))
        self.assertTrue(np.allclose(out_mri.segmentation_prob, self.mri.segmentation_prob, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        os.remove("mr_io_test_segmented_hpc_predict_in.h5")
        os.remove("mr_io_test_segmented_hpc_predict_out.h5")


if __name__ == '__main__':
    unittest.main()

