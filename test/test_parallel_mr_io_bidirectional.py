import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, FlowMRI, SegmentedFlowMRI
from test_common import spatial_hyperslab_dims_test, validate_group_name, validate_array

def write_hdf5_exec_fortran(test_inst):
    test_cls = type(test_inst)
 
    # Write HDF5 from Python
    test_inst.mri.write_hdf5(test_cls.filename_mri_in)
    

    ## Read HDF5 from Fortran in parallel
#   Run test
    fort_command = "%s %d %d %d %s %s " % \
                           (test_cls.filename_exec,
                            *test_inst.mpi_cart_dims,
                            test_cls.filename_mri_in,
                            test_cls.filename_mri_out) + \
                            (" %d %d %d %d %d %d " % (*test_cls.num_pad_vox_lhs, *test_cls.num_pad_vox_rhs) \
                             if hasattr(test_cls, "num_pad_vox_lhs") else "") + \
                            " 1> %s 2> %s " % (test_cls.filename_out_rank % ("${PMI_RANK}"),
                            test_cls.filename_err_rank % ("${PMI_RANK}"))
    fort = sp.run(["mpiexec","-np", "%d" % (test_cls.mpi_proc), \
                   "bash", "-c", fort_command],
                   stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    
#     # Debug test
#     fort_command = "fortran/test/%s %d %d %d %s %s" % \
#                            (test_cls.filename_exec,
#                             *test_inst.mpi_cart_dims,
#                             test_cls.filename_mri_in,
#                             test_cls.filename_mri_out) + \
#                             (" %d %d %d %d %d %d " % (*test_cls.num_pad_vox_lhs, *test_cls.num_pad_vox_rhs) \
#                              if hasattr(test_cls, "num_pad_vox_lhs") else "")
#     fort = sp.run(["mpiexec","-np", "%d" % (test_cls.mpi_proc), \
#                    "xterm", "-e", "gdb", "--args", 
#                    *fort_command.split(" ")],
#                    stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    
    print("Shell command returned out/err:")
    print(fort.stdout.decode("utf-8"))
    print(fort.stderr.decode("utf-8"))
         
    for mpi_rank in range(test_cls.mpi_proc):
        with open(test_cls.filename_err_rank % (mpi_rank), 'r') as err:
            print("Fortran command for rank %d returned err:" % (mpi_rank))
            print(err.read())


def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri_in, test_cls.filename_mri_out] + \
          [test_cls.filename_out_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)] + \
          [test_cls.filename_err_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)]
    for f in files:
        os.remove(f) 

class TestSpatialMRIBidirectional(unittest.TestCase):
        
    # number of Fortran MPI processes
    mpi_proc = 2**3
     
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_parallel_reader_writer"
   
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
     
    mri_group_name = "spatial-mri"
   
    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(91,31,71) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        #voxel_feature = np.reshape(np.array((9*[0]) + (9*[1])), (2,3,3))
        self.mri = SpatialMRI(voxel_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.voxel_feature)
           
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
               
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
   
        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.voxel_feature, out_mri.voxel_feature)        
   
    def tearDown(self):
        remove_test_files(self)
   
   
class TestSpaceTimeMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
        
    # number of Fortran MPI processes
    mpi_proc = 2**3
    
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_parallel_reader_writer_space_time"
   
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
     
    mri_group_name = "space-time-mri"
   
    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        voxel_feature = np.random.rand(67,43,29,11,3)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.voxel_feature)
   
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
               
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
   
        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.voxel_feature, out_mri.voxel_feature)
        
    def tearDown(self):
        remove_test_files(self)
    
    
class TestFlowMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
        
    # number of Fortran MPI processes
    mpi_proc = 2**3
    
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_parallel_reader_writer_flow"
   
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
     
    mri_group_name = "flow-mri"
    
    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        intensity = np.random.rand(67,43,29,11)        
        velocity_mean = np.random.rand(67,43,29,11,3)        
        velocity_cov = np.random.rand(67,43,29,11,3,5)        
        self.mri = FlowMRI(geometry, time, intensity, velocity_mean, velocity_cov)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
   
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
               
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
    mpi_proc = 2**3
    
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_parallel_reader_writer_segmented_flow"
  
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
    
    mri_group_name = "segmented-flow-mri"
  
    def setUp(self):
        # Initialize the MRI data
        # geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        intensity = np.random.rand(67,43,29,11)        
        velocity_mean = np.random.rand(67,43,29,11,3)        
        velocity_cov = np.random.rand(67,43,29,11,3,5)        
        segmentation_prob = np.random.rand(67,43,29,11)        
        self.mri = SegmentedFlowMRI(geometry, time, intensity, velocity_mean, velocity_cov, segmentation_prob)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
  
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
              
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
    mpi_proc = 2**3
  
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_parallel_reader_writer_flow_padded"

    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
  
    mri_group_name = "flow-mri"
    
    padding = (0.5, 0.4, 0.7)
  
    def setUp(self):
        test_cls = type(self)
        
        # Initialize the MRI data
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        intensity = np.random.rand(67,43,29,11)        
        velocity_mean = np.random.rand(67,43,29,11,3)        
        velocity_cov = np.random.rand(67,43,29,11,3,5)        
        self.mri = FlowMRI(geometry, time, intensity, velocity_mean, velocity_cov)
        
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
        
        test_cls.num_vox = intensity.shape[:3]
        
        num_pad_vox = [int(np.ceil(2.*test_cls.padding[i]*test_cls.num_vox[i])) for i in range(3)]
        
        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] -1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i]-test_cls.num_vox[i])//2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i]-test_cls.num_vox[i]+1)//2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i]//self.mpi_cart_dims[i] for i in range(3)]
      
  
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
             
        out_mri = type(self.mri).read_hdf5(type(self).filename_mri_out)
 
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
      
    def tearDown(self):
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

