import unittest
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, FlowMRI, SegmentedFlowMRI
from test_common import filename_mri, filename_out, filename_err, remove_test_files,\
                        validate_group_name,\
                        spatial_hyperslab_dims_test,\
                        spatial_hyperslab_loc_test,\
                        validate_spatial_fort_array,\
                        read_out, print_err


def validate_dist_spatial_fort_array(test_inst, mpi_rank, array, out_lines):
    test_cls = type(test_inst)

    fort_file_dims = np.fromstring(out_lines[0], dtype=int, sep=' ') 
    fort_offset = np.fromstring(out_lines[1], dtype=int, sep=' ') 
    fort_dims = np.fromstring(out_lines[2], dtype=int, sep=' ') 
    fort_array= np.fromstring(out_lines[3], dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose()

    hyperslab_offset, hyperslab_shape = spatial_hyperslab_loc_test(test_cls, mpi_rank, test_inst.mpi_cart_dims, array)

    test_inst.assertTrue(np.array_equal(fort_file_dims, array.shape))
    test_inst.assertTrue(np.array_equal(fort_offset, hyperslab_offset))
    test_inst.assertTrue(np.array_equal(fort_dims, hyperslab_shape))
    test_inst.assertTrue(np.allclose(fort_array, 
                    np.reshape(array[ \
                               hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                               hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                               hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2]], hyperslab_shape), rtol=1e-14))


def validate_dist_spacetime_fort_array(test_inst, mpi_rank, array, out_lines, transpose_dims):
    test_cls = type(test_inst)

    hyperslab_offset, hyperslab_shape = spatial_hyperslab_loc_test(test_cls, mpi_rank, test_inst.mpi_cart_dims, array)
        
    fort_file_dims = np.fromstring(out_lines[0], dtype=int, sep=' ') 
    fort_offset = np.fromstring(out_lines[1], dtype=int, sep=' ') 
    fort_dims = np.fromstring(out_lines[2], dtype=int, sep=' ') 
    
    fort_file_time_dim = np.fromstring(out_lines[3], dtype=int, sep=' ')[0]
    fort_time_offset = np.fromstring(out_lines[4], dtype=int, sep=' ')[0]
    fort_time_dim = np.fromstring(out_lines[5], dtype=int, sep=' ')[0]

    if len(transpose_dims) > 4:
        fort_vector_dims = np.fromstring(out_lines[6], dtype=int, sep=' ')
    else:
        fort_vector_dims = np.array([], dtype=int)
    data_line = 7 if len(transpose_dims) > 4 else 6
    fort_array = np.fromstring(out_lines[data_line], dtype=float, sep=' ')

    if hasattr(test_cls, "num_pad_vox_lhs"):
        fort_array_lbound = np.fromstring(out_lines[data_line+1], dtype=int, sep=' ')
        fort_array_ubound = np.fromstring(out_lines[data_line+2], dtype=int, sep=' ')
    
    test_inst.assertTrue(np.array_equal(fort_file_dims, array.shape[:3]))
    test_inst.assertTrue(fort_file_time_dim == array.shape[3])
    if len(transpose_dims) > 4:
        test_inst.assertTrue(np.array_equal(fort_vector_dims, array.shape[4:]))

    if hasattr(test_cls, "num_pad_vox_lhs"):
        if all(hyperslab_shape[:3] > (0, 0, 0)):
            hyperslab_lbound = [((test_cls.num_pad_vox_lhs[i] + hyperslab_offset[i]) % test_cls.num_vox_per_proc[i]) + 1 for i in range(3)]
            hyperslab_ubound = [hyperslab_lbound[i] + hyperslab_shape[i] - 1 for i in range(3)]
        else:
            hyperslab_lbound = [1,1,1]
            hyperslab_ubound = [0,0,0]
        hyperslab_lbound = np.array(hyperslab_lbound + [hyperslab_offset[i] + 1 for i in range(3,3+len(hyperslab_offset[3:]))])
        hyperslab_ubound = np.array(hyperslab_ubound + [hyperslab_offset[i] + hyperslab_shape[i] for i in range(3,3+len(hyperslab_offset[3:]))])
        
#         print("MPI rank:         {}".format(mpi_rank))
#         print("hyperslab_lbound: {} vs {}".format(hyperslab_lbound, fort_array_lbound))
#         print("hyperslab_ubound: {} vs {}".format(hyperslab_ubound, fort_array_ubound))

        test_inst.assertTrue(np.array_equal(hyperslab_lbound[:3], fort_array_lbound[-3:]))
        test_inst.assertTrue(np.array_equal(hyperslab_ubound[:3], fort_array_ubound[-3:]))

        test_inst.assertTrue( hyperslab_lbound[3], fort_array_lbound[-3] )
        test_inst.assertTrue( hyperslab_ubound[3], fort_array_ubound[-3] )

        if len(transpose_dims) > 4:
            test_inst.assertTrue(np.array_equal(hyperslab_lbound[4:], fort_array_lbound[:-4]))
            test_inst.assertTrue(np.array_equal(hyperslab_ubound[4:], fort_array_ubound[:-4]))
            
    if not hasattr(test_cls, "num_pad_vox_lhs") or (hyperslab_shape[:3] > np.zeros((3,),dtype=int)).all():    
        test_inst.assertTrue(np.array_equal(list(fort_offset) + [fort_time_offset] + [0]*len(fort_vector_dims), hyperslab_offset))
        test_inst.assertTrue(np.array_equal(np.concatenate((fort_dims, np.array([fort_time_dim-fort_time_offset]), fort_vector_dims)), hyperslab_shape))
        
        fort_array = fort_array.reshape(np.flip(np.concatenate((fort_vector_dims, np.array([fort_time_dim-fort_time_offset]), fort_dims)))).transpose(transpose_dims)
        if len(transpose_dims) == 4: # scalar
            test_inst.assertTrue(np.allclose(fort_array, 
                            np.reshape(array[ \
                                       hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                       hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                       hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                       hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3]], hyperslab_shape), rtol=1e-14))    
        elif len(transpose_dims) == 5: # vector
            test_inst.assertTrue(np.allclose(fort_array, 
                            np.reshape(array[ \
                                       hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                       hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                       hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                       hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3],
                                       hyperslab_offset[4]:hyperslab_offset[4]+hyperslab_shape[4]], hyperslab_shape), rtol=1e-14))    
        elif len(transpose_dims) == 6: # matrix
            test_inst.assertTrue(np.allclose(fort_array, 
                            np.reshape(array[ \
                                       hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                       hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                       hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                       hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3],
                                       hyperslab_offset[4]:hyperslab_offset[4]+hyperslab_shape[4],
                                       hyperslab_offset[5]:hyperslab_offset[5]+hyperslab_shape[5]], hyperslab_shape), rtol=1e-14))    
    else:
        test_inst.assertTrue(np.array_equal(fort_dims, hyperslab_shape[:3]))
        test_inst.assertTrue(len(fort_array) == 0)
        
        
def validate_dist_spacetime_scalar_fort_array(test_inst, mpi_rank, array, out_lines):
    validate_dist_spacetime_fort_array(test_inst, mpi_rank, array, out_lines, transpose_dims=(2,1,0,3))

def validate_dist_spacetime_vector_fort_array(test_inst, mpi_rank, array, out_lines):
    validate_dist_spacetime_fort_array(test_inst, mpi_rank, array, out_lines, transpose_dims=(2,1,0,3,4))

def validate_dist_spacetime_matrix_fort_array(test_inst, mpi_rank, array, out_lines):
    validate_dist_spacetime_fort_array(test_inst, mpi_rank, array, out_lines, transpose_dims=(2,1,0,3,5,4))


class TestSpatialMRI(unittest.TestCase):

    fortran_exec = "mr_io_test_parallel_reader_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**4
    spatial_feature_shape=(91,31,71)

    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        self.mri = SpatialMRI(scalar_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.scalar_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])
            validate_dist_spatial_fort_array(self, mpi_rank, self.mri.scalar_feature, out_lines[1:5])
            
    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRI(unittest.TestCase):

    fortran_exec = "mr_io_test_parallel_reader_space_time_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)

    def setUp(self):
        # Initialize the MRI datan
        time = np.random.rand(7)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        vector_feature = np.random.rand(*type(self).spatial_feature_shape, 7, 3)
        self.mri = SpaceTimeMRI(geometry, time, vector_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.vector_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])

            validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])

            validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.vector_feature, out_lines[9:17])

    def tearDown(self):
        remove_test_files(self)


class TestFlowMRI(unittest.TestCase):  # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_flow_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(7)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 7)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 7, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 7, 3, 5)
        self.mri = FlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])

            validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])

            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.intensity, out_lines[9:16])
            validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.velocity_mean, out_lines[16:24])
            validate_dist_spacetime_matrix_fort_array(self, mpi_rank, self.mri.velocity_cov, out_lines[24:32])

    def tearDown(self):
        remove_test_files(self)


class TestSegmentedFlowMRI(unittest.TestCase):  # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_segmented_flow_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(7)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 7)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 7, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 7, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 7)
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov,
                                    segmentation_prob)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])

            validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])

            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.intensity, out_lines[9:16])
            validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.velocity_mean, out_lines[16:24])
            validate_dist_spacetime_matrix_fort_array(self, mpi_rank, self.mri.velocity_cov, out_lines[24:32])
            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.segmentation_prob, out_lines[32:39])

    def tearDown(self):
        remove_test_files(self)


class TestFlowMRIPadded(unittest.TestCase):  # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_flow_padded_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)
    padding = (1.5, 1.4, 1.7)

    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(7)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 7)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 7, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 7, 3, 5)
        self.mri = FlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov)

        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)

        test_cls.num_vox = test_cls.spatial_feature_shape

        num_pad_vox = [int(np.ceil(2. * test_cls.padding[i] * test_cls.num_vox[i])) for i in range(3)]

        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] - 1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i] - test_cls.num_vox[i]) // 2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i] - test_cls.num_vox[i] + 1) // 2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i] // self.mpi_cart_dims[i] for i in range(3)]

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])

            validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])

            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.intensity, out_lines[9:18])  # 7+2
            validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.velocity_mean,
                                                      out_lines[18:28])  # 8+2
            validate_dist_spacetime_matrix_fort_array(self, mpi_rank, self.mri.velocity_cov,
                                                      out_lines[28:38])  # 8+2

    def tearDown(self):
        remove_test_files(self)


class TestSegmentedFlowMRIPadded(unittest.TestCase):  # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_segmented_flow_padded_container"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)
    padding = (1.5, 1.4, 1.7)

    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(7)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 7)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 7, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 7, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 7)
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov,
                                    segmentation_prob)

        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)

        test_cls.num_vox = intensity.shape[:3]

        num_pad_vox = [int(np.ceil(2. * test_cls.padding[i] * test_cls.num_vox[i])) for i in range(3)]

        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] - 1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i] - test_cls.num_vox[i]) // 2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i] - test_cls.num_vox[i] + 1) // 2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i] // self.mpi_cart_dims[i] for i in range(3)]

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran
        self.mri.write_hdf5(filename_mri(self))

        for mpi_rank in range(type(self).mpi_proc):
            out_lines = read_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

            validate_group_name(self, out_lines[0])

            validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])

            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.intensity, out_lines[9:18])  # 7+2
            validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.velocity_mean,
                                                      out_lines[18:28])  # 8+2
            validate_dist_spacetime_matrix_fort_array(self, mpi_rank, self.mri.velocity_cov,
                                                      out_lines[28:38])  # 8+2
            validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.segmentation_prob,
                                                      out_lines[38:47])  # 7+2

    def tearDown(self):
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

