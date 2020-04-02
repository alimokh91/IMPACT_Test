import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, \
                  SpaceTimeMRI, \
                  FlowMRI, \
                  SegmentedFlowMRI
from test_common import spatial_hyperslab_dims_test, \
                        validate_group_name, \
                        validate_array, \
                        validate_replicated_mri_coordinates, \
                        validate_replicated_mri_vector_array


def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri_in, test_cls.filename_mri_out] + \
          [test_cls.filename_out_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)] + \
          [test_cls.filename_err_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)]
    for f in files:
        os.remove(f)


def get_fortran_args(test_cls):
    from mr_io_domain_decomp import spatial_hyperslab_dims
    mpi_cart_dims = spatial_hyperslab_dims(test_cls.mpi_proc, test_cls.spatial_feature_shape)

    fort_args = "%d %d %d %s %s" % \
                   (*mpi_cart_dims,
                    test_cls.filename_mri_in,
                    test_cls.filename_mri_out)

    if test_cls in [TestFlowMRIPaddedBidirectional, TestSegmentedFlowMRIPaddedBidirectional,
                    TestFlowMRIPaddedToSpaceTimeBidirectional,
                    TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional]:
        test_cls.num_vox = test_cls.spatial_feature_shape

        num_pad_vox = [int(np.ceil(2. * test_cls.padding[i] * test_cls.num_vox[i])) for i in range(3)]

        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + mpi_cart_dims[i] - 1) // mpi_cart_dims[i]) \
                       * mpi_cart_dims[i] for i in range(3)]
        num_pad_vox_lhs = [(num_ext_vox[i] - test_cls.num_vox[i]) // 2 for i in range(3)]
        num_pad_vox_rhs = [(num_ext_vox[i] - test_cls.num_vox[i] + 1) // 2 for i in range(3)]
        num_vox_per_proc = [num_ext_vox[i] // mpi_cart_dims[i] for i in range(3)]

        fort_args += " %d %d %d %d %d %d " % (*num_pad_vox_lhs, *num_pad_vox_rhs)

    if test_cls in [TestFlowMRIPaddedToSpaceTimeBidirectional,
                    TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional]:
        fort_args += " %d %d %d %d " % (test_cls.tr, *test_cls.sr)

    fort_args += " 1> %s 2> %s" % \
                 (test_cls.filename_out_rank % ("${PMI_RANK}"),
                  test_cls.filename_err_rank % ("${PMI_RANK}"))

    return fort_args


class TestSpatialMRIBidirectional(unittest.TestCase):
         
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS']) #2**3
    spatial_feature_shape=(91,31,71)
      
    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"
      
    mri_group_name = "spatial-mri"
    
    def setUp(self):
        # Initialize the MRI data
        scalar_feature = np.random.rand(*type(self).spatial_feature_shape) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        self.mri = SpatialMRI(scalar_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.scalar_feature)
            
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        self.mri.write_hdf5(type(self).filename_mri_in)

        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)

        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.scalar_feature, out_mri.scalar_feature)

    def tearDown(self):
        remove_test_files(self)
    
    
class TestSpaceTimeMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
         
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS']) #2**3
    spatial_feature_shape=(67,43,29)
     
    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"

    mri_group_name = "space-time-mri"
    
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
        self.mri.write_hdf5(type(self).filename_mri_in)
                
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
    
        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.vector_feature, out_mri.vector_feature)

    def tearDown(self):
        remove_test_files(self)
     
     
class TestFlowMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
         
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS']) #2**3
    spatial_feature_shape=(67,43,29)

    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"

    mri_group_name = "flow-mri"
     
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
        self.mri = FlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)


    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        self.mri.write_hdf5(type(self).filename_mri_in)
                
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
    
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)

    def tearDown(self):
        remove_test_files(self)
    
   
class TestSegmentedFlowMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
        
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS']) #2**3
    spatial_feature_shape=(67,43,29)
     
    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"

    mri_group_name = "segmented-flow-mri"
   
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 11)
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov,
                                    segmentation_prob)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
   
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        self.mri.write_hdf5(type(self).filename_mri_in)
               
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
   
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
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS']) #2**3
    spatial_feature_shape=(67,43,29)
   
    # Filenames
    filename_mri_in  = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"
   
    mri_group_name = "flow-mri"
     
    padding = (0.5, 0.4, 0.7)
   
    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
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
        self.mri.write_hdf5(type(self).filename_mri_in)
              
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
  
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
       
    def tearDown(self):
        remove_test_files(self)    
 
 
class TestFlowMRIPaddedToSpaceTimeBidirectional(unittest.TestCase): # FIXME: coordinates test...
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS'])  # 2**3
    spatial_feature_shape = (67, 43, 29)

    # Filenames
    filename_mri_in = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"

    mri_group_name = "flow-mri"
     
    padding = (0.5, 0.4, 0.7)

    tr = 3
    sr = (2, 3, 3)

    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
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
        self.mri.write_hdf5(type(self).filename_mri_in)

        out_mri = SpaceTimeMRI.read_hdf5(type(self).filename_mri_out)
         
        validate_replicated_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.vector_feature)
       
    def tearDown(self):
        remove_test_files(self)


class TestSegmentedFlowMRIPaddedBidirectional(unittest.TestCase): # FIXME: coordinates test...
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS'])  # 2**3
    spatial_feature_shape = (67, 43, 29)

    # Filenames
    filename_mri_in = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"
  
    mri_group_name = "segmented-flow-mri"
    
    padding = (0.5, 0.4, 0.7)
  
    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 11)
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
        self.mri.write_hdf5(type(self).filename_mri_in)
             
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
 
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
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_MPI_FORTRAN_PROCS'])  # 2**3
    spatial_feature_shape = (67, 43, 29)

    # Filenames
    filename_mri_in = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_in.h5"  # + ".h5"
    filename_mri_out = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_out.h5"  # + ".h5"
    filename_out_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.out"
    filename_err_rank = os.environ['HPC_PREDICT_IO_TEST_MRI_PATH'][:-3] + "_%s.err"
   
    mri_group_name = "segmented-flow-mri"
     
    padding = (0.5, 0.4, 0.7)

    tr = 3
    sr = (2, 3, 3)

    def setUp(self):
        test_cls = type(self)

        # Initialize the MRI data
        time = np.random.rand(11)
        time_heart_cycle_period = np.max(time)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        geometry = [np.random.rand(type(self).spatial_feature_shape[0]), \
                    np.random.rand(type(self).spatial_feature_shape[1]), \
                    np.random.rand(type(self).spatial_feature_shape[2])]
        intensity = np.random.rand(*type(self).spatial_feature_shape, 11)
        velocity_mean = np.random.rand(*type(self).spatial_feature_shape, 11, 3)
        velocity_cov = np.random.rand(*type(self).spatial_feature_shape, 11, 3, 5)
        segmentation_prob = np.random.rand(*type(self).spatial_feature_shape, 11)
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
        self.mri.write_hdf5(type(self).filename_mri_in)

        out_mri = SpaceTimeMRI.read_hdf5(type(self).filename_mri_out)
         
        validate_replicated_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.vector_feature)
       
    def tearDown(self):
        remove_test_files(self)

if __name__ == '__main__':
    unittest.main()

