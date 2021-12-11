import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, \
                  SpaceTimeMRI,\
                  FlowMRI, \
                  SegmentedFlowMRI
from test_common import filename_mri_in, filename_mri_out, remove_test_files, \
                        validate_group_name, validate_array


class TestSpatialMRIBidirectional(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_writer"
    spatial_feature_shape = (71,61,97)

    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape)
        self.mri = SpatialMRI(scalar_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri_in(self))

        out_mri = SpatialMRI.read_hdf5(filename_mri_out(self))
        
        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.scalar_feature, out_mri.scalar_feature)        
        
    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRIBidirectional(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_writer_space_time"
    spatial_feature_shape = (26,11,7)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(13)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        vector_feature = np.random.rand(*type(self).spatial_feature_shape, 13, 19)
        self.mri = SpaceTimeMRI(geometry, time, vector_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri_in(self))
            
        out_mri = SpaceTimeMRI.read_hdf5(filename_mri_out(self))

        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.vector_feature, out_mri.vector_feature)

    def tearDown(self):
        remove_test_files(self)
 
 
class TestFlowMRIBidirectional(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_writer_flow"
    spatial_feature_shape = (23,13,19)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 17)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 17, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 17, 3, 5)
        self.mri = FlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov)
 
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri_in(self))
            
        out_mri = FlowMRI.read_hdf5(filename_mri_out(self))

        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
 
    def tearDown(self):
        remove_test_files(self)
         
 
class TestSegmentedFlowMRIBidirectional(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_writer_segmented_flow"
    spatial_feature_shape = (23,13,19)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 17)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 17, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 17, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 17)
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean,
                                    velocity_cov, segmentation_prob)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri_in(self))
            
        out_mri = SegmentedFlowMRI.read_hdf5(filename_mri_out(self))
     
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
        validate_array(self, self.mri.segmentation_prob, out_mri.segmentation_prob)
 
    def tearDown(self):
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

