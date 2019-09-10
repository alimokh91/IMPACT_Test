import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, FlowMRI, SegmentedFlowMRI
from test_common import spatial_hyperslab_dims_test, spatial_hyperslab_loc_test, validate_spatial_fort_array


def write_hdf5_read_in_fortran(test_inst):
    test_cls = type(test_inst)
 
    # Write HDF5 from Python
    test_inst.mri.write_hdf5(test_cls.filename_mri)

    ## Read HDF5 from Fortran in parallel
#   Run test
    fort_command = "fortran/%s %d %d %d %s " % \
                           (test_cls.filename_exec,
                            *test_inst.mpi_cart_dims,
                            test_cls.filename_mri) + \
                            (" %d %d %d %d %d %d " % (*test_cls.num_pad_vox_lhs, *test_cls.num_pad_vox_rhs) \
                             if hasattr(test_cls, "num_pad_vox_lhs") else "") + \
                            " 1> %s 2> %s" % \
                            (test_cls.filename_out_rank % ("${PMI_RANK}"),
                            test_cls.filename_err_rank % ("${PMI_RANK}"))
    fort = sp.run(["mpiexec","-np", "%d" % (test_cls.mpi_proc), \
                   "bash", "-c", fort_command],
                   stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    
#     # Debug test
#     fort_command = "fortran/%s %d %d %d %s" % \
#                             (test_cls.filename_exec,
#                              *test_inst.mpi_cart_dims,
#                              test_cls.filename_mri) + \
#                             (" %d %d %d %d %d %d " % (*test_cls.num_pad_vox_lhs, *test_cls.num_pad_vox_rhs) \
#                              if hasattr(test_cls, "num_pad_vox_lhs") else "")
#     print(fort_command)
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


def validate_group_name(test_inst, fort_group_name):
    test_inst.assertEqual(fort_group_name, type(test_inst).mri_group_name)        
    

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
    fort_array = np.fromstring(out_lines[7 if len(transpose_dims) > 4 else 6], dtype=float, sep=' ')
    
    test_inst.assertTrue(np.array_equal(fort_file_dims, array.shape[:3]))
    test_inst.assertTrue(fort_file_time_dim == array.shape[3])
    if len(transpose_dims) > 4:
        test_inst.assertTrue(np.array_equal(fort_vector_dims, array.shape[4:]))
    
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


def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri] + \
          [test_cls.filename_out_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)] + \
          [test_cls.filename_err_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)]
    for f in files:
        os.remove(f) 


class TestSpatialMRI(unittest.TestCase):
       
    # number of Fortran MPI processes
    mpi_proc = 2**4
   
    # Filenames
    filename_prefix = "mr_io_test_parallel_reader"
 
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
   
    mri_group_name = "spatial-mri"
 
    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(91,31,71) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        self.mri = SpatialMRI(voxel_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.voxel_feature)
 
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self)
   
        for mpi_rank in range(type(self).mpi_proc):
            filename_out = type(self).filename_out_rank % (mpi_rank)
   
            with open(filename_out,'r') as out:
                out_lines = [l.strip(' \n') for l in out.readlines()]
 
                validate_group_name(self, out_lines[0])
                 
                validate_dist_spatial_fort_array(self, mpi_rank, self.mri.voxel_feature, out_lines[1:5])
           
    def tearDown(self):
        remove_test_files(self)
  
  
class TestSpaceTimeMRI(unittest.TestCase): # FIXME: coordinates test...
       
    # number of Fortran MPI processes
    mpi_proc = 2**3
   
    # Filenames
    filename_prefix = "mr_io_test_parallel_reader_space_time"
 
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
   
    mri_group_name = "space-time-mri"
   
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        voxel_feature = np.random.rand(67,43,29,11,3)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.voxel_feature)
 
    def test_communicator(self):  
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self)
   
        for mpi_rank in range(type(self).mpi_proc):
            filename_out = type(self).filename_out_rank % (mpi_rank)
   
            with open(filename_out,'r') as out:
                out_lines = [l.strip(' \n') for l in out.readlines()]
   
                validate_group_name(self, out_lines[0])
 
                validate_spatial_fort_array(self, self.mri.time, out_lines[1:3])
                validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[3:5])
                validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[5:7])
                validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[7:9])
                 
                validate_dist_spacetime_vector_fort_array(self, mpi_rank, self.mri.voxel_feature, out_lines[9:17])
   
       
    def tearDown(self):
        remove_test_files(self)
               
   
class TestFlowMRI(unittest.TestCase): # FIXME: coordinates test...
    # number of Fortran MPI processes
    mpi_proc = 2**3
   
    # Filenames
    filename_prefix = "mr_io_test_parallel_reader_flow"
 
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
   
    mri_group_name = "flow-mri"
   
   
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        intensity = np.random.rand(67,43,29,11)        
        velocity_mean = np.random.rand(67,43,29,11,3)        
        velocity_cov = np.random.rand(67,43,29,11,3,5)        
        self.mri = FlowMRI(geometry, time, intensity, velocity_mean, velocity_cov)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
 
     
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self)
   
        for mpi_rank in range(type(self).mpi_proc):
            filename_out = type(self).filename_out_rank % (mpi_rank)
   
            with open(filename_out,'r') as out:
                out_lines = [l.strip(' \n') for l in out.readlines()]
   
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
              
 
class TestSegmentedFlowMRI(unittest.TestCase): # FIXME: coordinates test...
    # number of Fortran MPI processes
    mpi_proc = 2**3
   
    # Filenames
    filename_prefix = "mr_io_test_parallel_reader_segmented_flow"
 
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
   
    mri_group_name = "segmented-flow-mri"
   
   
    def setUp(self):
        # Initialize the MRI data
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
        write_hdf5_read_in_fortran(self)
   
        for mpi_rank in range(type(self).mpi_proc):
            filename_out = type(self).filename_out_rank % (mpi_rank)
   
            with open(filename_out,'r') as out:
                out_lines = [l.strip(' \n') for l in out.readlines()]
   
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

            
class TestFlowMRIPadded(unittest.TestCase): # FIXME: coordinates test...
    # number of Fortran MPI processes
    mpi_proc = 2**3
  
    # Filenames
    filename_prefix = "mr_io_test_parallel_reader_flow_padded"

    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
  
    mri_group_name = "flow-mri"
    
    padding = (1.5, 1.4, 1.7)
  
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
        write_hdf5_read_in_fortran(self)
  
        for mpi_rank in range(type(self).mpi_proc):
            filename_out = type(self).filename_out_rank % (mpi_rank)
  
            with open(filename_out,'r') as out:
                out_lines = [l.strip(' \n') for l in out.readlines()]

                validate_dist_spacetime_scalar_fort_array(self, mpi_rank, self.mri.intensity, out_lines[9:16])
  
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
            
            
if __name__ == '__main__':
    unittest.main()

