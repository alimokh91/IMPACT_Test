import unittest
import os
import numpy as np
from mr_io import SpatialMRI,\
                  SpaceTimeMRI,\
                  FlowMRI,\
                  SegmentedFlowMRI
from test_common import filename_mri, filename_out, filename_err, remove_test_files,\
                        validate_group_name,\
                        validate_spatial_fort_array,\
                        validate_spacetime_scalar_fort_array,\
                        validate_spacetime_vector_fort_array,\
                        validate_spacetime_matrix_fort_array,\
                        read_out, print_err


class TestSpatialMRI(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_container"
    spatial_feature_shape = (83,63,93)

    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape)
        self.mri = SpatialMRI(scalar_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        out_lines = read_out(filename_out(self))
        print_err(filename_err(self))

        validate_group_name(self, out_lines[0])
        validate_spatial_fort_array(self, self.mri.scalar_feature, out_lines[1:3])

    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRI(unittest.TestCase):

    fortran_exec = "mr_io_test_reader_space_time_container"
    spatial_feature_shape = (23, 13, 19)

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
        self.mri.write_hdf5(filename_mri(self))

        out_lines = read_out(filename_out(self))
        print_err(filename_err(self))

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

    fortran_exec = "mr_io_test_reader_flow_container"
    spatial_feature_shape = (23, 13, 19)

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
        self.mri.write_hdf5(filename_mri(self))

        out_lines = read_out(filename_out(self))
        print_err(filename_err(self))

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

    fortran_exec = "mr_io_test_reader_segmented_flow_container"
    spatial_feature_shape = (23, 13, 19)

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
        self.mri.write_hdf5(filename_mri(self))

        out_lines = read_out(filename_out(self))
        print_err(filename_err(self))

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

