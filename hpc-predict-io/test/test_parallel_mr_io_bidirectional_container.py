import unittest
import os
import numpy as np
from mr_io import SpatialMRI, \
                  SpaceTimeMRI, \
                  FlowMRI, \
                  SegmentedFlowMRI
from test_common import filename_mri_in, filename_mri_out, filename_out, filename_err, remove_test_files,\
                        spatial_hyperslab_dims_test, \
                        validate_group_name, \
                        validate_array, \
                        validate_replicated_mri_coordinates, \
                        validate_replicated_mri_vector_array


class TestSpatialMRIBidirectional(unittest.TestCase):

    fortran_exec = "mr_io_test_parallel_reader_writer"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape=(91,31,71)

    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        self.mri = SpatialMRI(scalar_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.scalar_feature)
            
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        self.mri.write_hdf5(filename_mri_in(self))

        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))

        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.scalar_feature, out_mri.scalar_feature)

    def tearDown(self):
        remove_test_files(self)
    
    
class TestSpaceTimeMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_space_time"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape = (17, 23, 19)

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(2)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        vector_feature = np.random.rand(*type(self).spatial_feature_shape, 2, 3)
        self.mri = SpaceTimeMRI(geometry, time, vector_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.vector_feature)

    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        self.mri.write_hdf5(filename_mri_in(self))
                
        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))
    
        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.vector_feature, out_mri.vector_feature)

    def tearDown(self):
        remove_test_files(self)
     
     
class TestFlowMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_flow"
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
        self.mri.write_hdf5(filename_mri_in(self))
                
        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))
    
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)

    def tearDown(self):
        remove_test_files(self)
    
   
class TestSegmentedFlowMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_segmented_flow"
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
        self.mri.write_hdf5(filename_mri_in(self))
               
        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))
   
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
        validate_array(self, self.mri.segmentation_prob, out_mri.segmentation_prob)
        
    def tearDown(self):
        remove_test_files(self)
 
 
class TestFlowMRIPaddedBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_flow_padded"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) #2**3
    spatial_feature_shape=(17,23,19)

    padding = (0.5, 0.4, 0.7)
   
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
        self.mri.write_hdf5(filename_mri_in(self))
              
        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))
  
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
       
    def tearDown(self):
        remove_test_files(self)    
 
 
class TestFlowMRIPaddedToSpaceTimeBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_flow_padded_to_space_time"
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE'])  # 2**3
    spatial_feature_shape = (17, 23, 19)

    padding = (0.5, 0.4, 0.7)

    tr = 3
    sr = (2, 3, 3)

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
        self.mri.write_hdf5(filename_mri_in(self))

        out_mri = SpaceTimeMRI.read_hdf5(filename_mri_out(self))
         
        validate_replicated_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.vector_feature)
       
    def tearDown(self):
        remove_test_files(self)


class TestSegmentedFlowMRIPaddedBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_segmented_flow_padded"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE'])  # 2**3
    spatial_feature_shape = (17, 23, 19)

    padding = (0.5, 0.4, 0.7)
  
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
        self.mri.write_hdf5(filename_mri_in(self))
             
        out_mri = type(self.mri).read_hdf5(filename_mri_out(self))
 
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
        validate_array(self, self.mri.segmentation_prob, out_mri.segmentation_prob)
      
    def tearDown(self):
        remove_test_files(self)  
        

class TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = "mr_io_test_parallel_reader_writer_segmented_flow_padded_to_space_time"
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE'])  # 2**3
    spatial_feature_shape = (17, 23, 19)

    padding = (0.5, 0.4, 0.7)

    tr = 3
    sr = (2, 3, 3)

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
        self.mri.write_hdf5(filename_mri_in(self))

        out_mri = SpaceTimeMRI.read_hdf5(filename_mri_out(self))
         
        validate_replicated_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.vector_feature)
       
    def tearDown(self):
        remove_test_files(self)

if __name__ == '__main__':
    unittest.main()

