import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI,\
                  SpaceTimeMRI,\
                  FlowMRI,\
                  SegmentedFlowMRI
from test_common import validate_group_name,\
                        validate_spatial_fort_array,\
                        validate_spacetime_scalar_fort_array,\
                        validate_spacetime_vector_fort_array,\
                        validate_spacetime_matrix_fort_array,\
                        read_out, print_err


def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri, test_cls.filename_err, test_cls.filename_out]
    for f in files:
        os.remove(f)

def get_fortran_args(test_cls):
    fort_args = "%s  1> %s 2> %s" % \
                  (test_cls.filename_mri,
                   test_cls.filename_out,
                   test_cls.filename_err)
    return fort_args


class TestSpatialMRI(unittest.TestCase):
  
    spatial_feature_shape=(83,63,93)

    # Filenames
    filename_mri = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH']  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"
 
    mri_group_name = "spatial-mri"
 
    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape)
        self.mri = SpatialMRI(scalar_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(type(self).filename_mri)

        out_lines = read_out(type(self).filename_out)
        print_err(type(self).filename_err)

        validate_group_name(self, out_lines[0])
        validate_spatial_fort_array(self, self.mri.scalar_feature, out_lines[1:3])

    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRI(unittest.TestCase):

    spatial_feature_shape=(23, 13, 19)

    # Filenames
    filename_mri = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH']  # + ".h5"
    filename_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".out"
    filename_err = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + ".err"

    mri_group_name = "space-time-mri"

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        vector_feature = np.random.rand(*type(self).spatial_feature_shape, 17, 21)
        self.mri = SpaceTimeMRI(geometry, time, vector_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(type(self).filename_mri)

        out_lines = read_out(type(self).filename_out)
        print_err(type(self).filename_err)

        validate_group_name(self, out_lines[0])

        validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
        validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
        validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
        validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])

        validate_spacetime_vector_fort_array(self, self.mri.vector_feature, out_lines[9:11])

    def tearDown(self):
        # Clean up file
        remove_test_files(self)


class TestFlowMRI(unittest.TestCase):

    spatial_feature_shape=(23, 13, 19)

    # Filenames
    filename_mri = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH']  # + ".h5"
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
        self.mri.write_hdf5(type(self).filename_mri)

        out_lines = read_out(type(self).filename_out)
        print_err(type(self).filename_err)

        validate_group_name(self, out_lines[0])

        validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
        validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
        validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
        validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])

        validate_spacetime_scalar_fort_array(self, self.mri.intensity, out_lines[9:11])
        validate_spacetime_vector_fort_array(self, self.mri.velocity_mean, out_lines[11:13])
        validate_spacetime_matrix_fort_array(self, self.mri.velocity_cov, out_lines[13:15])

    def tearDown(self):
        # Clean up file
        remove_test_files(self)


class TestSegmentedFlowMRI(unittest.TestCase):

    spatial_feature_shape=(23, 13, 19)

    # Filenames
    filename_mri = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH']  # + ".h5"
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
        self.mri.write_hdf5(type(self).filename_mri)

        out_lines = read_out(type(self).filename_out)
        print_err(type(self).filename_err)

        validate_group_name(self, out_lines[0])

        validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
        validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
        validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
        validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])

        validate_spacetime_scalar_fort_array(self, self.mri.intensity, out_lines[9:11])
        validate_spacetime_vector_fort_array(self, self.mri.velocity_mean, out_lines[11:13])
        validate_spacetime_matrix_fort_array(self, self.mri.velocity_cov, out_lines[13:15])
        validate_spacetime_scalar_fort_array(self, self.mri.segmentation_prob, out_lines[15:17])

    def tearDown(self):
        # Clean up file
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

