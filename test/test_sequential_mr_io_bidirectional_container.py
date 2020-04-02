import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, \
                  SpaceTimeMRI,\
                  FlowMRI, \
                  SegmentedFlowMRI
from test_common import validate_group_name, \
                        validate_array


def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri_in, test_cls.filename_mri_out, test_cls.filename_err, test_cls.filename_out]
    for f in files:
        os.remove(f)

def get_fortran_args(test_cls):
    fort_args = "%s %s  1> %s 2> %s" % \
                  (test_cls.filename_mri_in,
                   test_cls.filename_mri_out,
                   test_cls.filename_out,
                   test_cls.filename_err)
    return fort_args


class TestSpatialMRIBidirectional(unittest.TestCase):

    spatial_feature_shape = (71,61,97)

    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"

    mri_group_name = "spatial-mri"

    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape)
        self.mri = SpatialMRI(scalar_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(type(self).filename_mri_in)

        out_mri = SpatialMRI.read_hdf5(type(self).filename_mri_out)
        
        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.scalar_feature, out_mri.scalar_feature)        
        
    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRIBidirectional(unittest.TestCase):

    spatial_feature_shape = (26,11,7)

    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"

    mri_group_name = "space-time-mri"

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
        self.mri.write_hdf5(type(self).filename_mri_in)
            
        out_mri = SpaceTimeMRI.read_hdf5(type(self).filename_mri_out)

        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.vector_feature, out_mri.vector_feature)

    def tearDown(self):
        remove_test_files(self)
 
 
class TestFlowMRIBidirectional(unittest.TestCase):
 
    spatial_feature_shape = (23,13,19)

    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"

    mri_group_name = "flow-mri"

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
        self.mri.write_hdf5(type(self).filename_mri_in)
            
        out_mri = FlowMRI.read_hdf5(type(self).filename_mri_out)

        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
 
    def tearDown(self):
        remove_test_files(self)
         
 
class TestSegmentedFlowMRIBidirectional(unittest.TestCase):
 
    spatial_feature_shape = (23,13,19)

    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"

    mri_group_name = "segmented-flow-mri"

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
        self.mri.write_hdf5(type(self).filename_mri_in)
            
        out_mri = SegmentedFlowMRI.read_hdf5(type(self).filename_mri_out)
     
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

