import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, HPCPredictMRI


def spatial_hyperslab_dims(cls, voxel_feature): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
    #import pdb; pdb.set_trace()
    mpi_size = cls.mpi_proc
    dims_mem = np.array(voxel_feature.shape[:3])
    block_dims = np.array((1,1,1))
    while (mpi_size > 1):
        refinement_axis = np.argmax(dims_mem)
        block_dims[refinement_axis] *= 2
        dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
        mpi_size //= 2
    dims_mem_boundary = np.array([ coord_shape % dims_mem[i] if coord_shape % dims_mem[i] != 0 else dims_mem[i] for i,coord_shape in enumerate(voxel_feature.shape[:3]) ]) 
    return block_dims, dims_mem, dims_mem_boundary

class TestSpatialMRI(unittest.TestCase):
     
    # number of Fortran MPI processes
    mpi_proc = 2**7
 
    # Filenames
    filename_prefix = "mr_io_test_parallel"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
 
    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(91,31,71) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        #voxel_feature = np.reshape(np.array((9*[0]) + (9*[1])), (2,3,3))
        #print(voxel_feature)
        self.mri = SpatialMRI(voxel_feature)
#         block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestSpatialMRI, self.mri.voxel_feature)
#         print(block_dims)
#         print(dims_mem)
#         print(dims_mem_boundary)
         
 
    def test_communicator(self):
 
        # Write HDF5 from Python
        self.mri.write_hdf5(TestSpatialMRI.filename_full)
 
        ## Read HDF5 from Fortran
        fort_command = "fortran/mr_io_test_parallel_reader %s  1> %s 2> %s" % \
                               (TestSpatialMRI.filename_full,
                                TestSpatialMRI.filename_out_rank % ("${PMI_RANK}"),
                                TestSpatialMRI.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestSpatialMRI.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
 
        # Read HDF5 from Fortran
        #fort = sp.run(["mpiexec","-np", "%d" % (TestSpatialMRI.mpi_proc), \
        #               "bash","-c","echo  \"spatial-mri\n 1 3 3\n %s \n \" 1> %s 2> %s" % \
        #                       (" ".join(9 * ["${PMI_RANK}"]), 
        #                        TestSpatialMRI.filename_out_rank % ("${PMI_RANK}"), 
        #                        TestSpatialMRI.filename_err_rank % ("${PMI_RANK}"))], 
        #                        stdout=sp.PIPE, stderr=sp.PIPE, check=True)
 
         
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)
        #import pdb; pdb.set_trace()
 
        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestSpatialMRI, self.mri.voxel_feature)
 
        for mpi_rank in range(TestSpatialMRI.mpi_proc):
            filename_out = TestSpatialMRI.filename_out_rank % (mpi_rank)
 
            with open(filename_out,'r') as out:
                out_lines = out.readlines()
                fort_full_dims_str = out_lines[1][:-1]
                fort_offset_str = out_lines[2][:-1] 
                fort_dims_str = out_lines[3][:-1]
                fort_array_str = out_lines[4][:-1]
                #print(" *** " + filename_out + " *** ")
                #print(fort_dims_str)
                #print(fort_array_str)
                fort_full_dims = np.fromstring(fort_full_dims_str, dtype=int, sep=' ')
                fort_offset = np.fromstring(fort_offset_str, dtype=int, sep=' ')
                fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
                fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose()
 
                #import pdb; pdb.set_trace()
                block_id = (mpi_rank % block_dims[0], (mpi_rank // block_dims[0]) % block_dims[1], mpi_rank // (block_dims[0]* block_dims[1]) )
                hyperslab_offset = np.array([dims_mem[i]*block_id[i] for i in range(3)])
                hyperslab_shape = np.array([dims_mem[i] if block_id[i] + 1 < block_dims[i] else dims_mem_boundary[i] for i in range(3)])
                #mpi_size = TestSpatialMRI.mpi_proc
                #hyperslab_shape = self.mri.voxel_feature.shape[:2] + ( (self.mri.voxel_feature.shape[2]+mpi_size-1)//mpi_size if mpi_rank + 1 != mpi_size else (self.mri.voxel_feature.shape[2] % mpi_size if self.mri.voxel_feature.shape[2] % mpi_size !=0 else self.mri.voxel_feature.shape[2] / mpi_size ),)
                #hyperslab_offset = (0,0) + (mpi_rank*((self.mri.voxel_feature.shape[2]+mpi_size-1)//mpi_size),)
                self.assertTrue(np.array_equal(fort_full_dims, self.mri.voxel_feature.shape))
                self.assertTrue(np.array_equal(fort_offset, hyperslab_offset))
                self.assertTrue(np.array_equal(fort_dims, hyperslab_shape))
                #import pdb; pdb.set_trace()
                self.assertTrue(np.allclose(fort_array, 
                                np.reshape(self.mri.voxel_feature[ \
                                           hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                           hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                           hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2]], hyperslab_shape), rtol=1e-14))
         
    def tearDown(self):
        # Clean up file
        files = [TestSpatialMRI.filename_full] + \
                [TestSpatialMRI.filename_out_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRI.mpi_proc)] + \
                [TestSpatialMRI.filename_err_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRI.mpi_proc)]
        for f in files:
            os.remove(f) 



class TestSpaceTimeMRI(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 2**7

    # Filenames
    filename_prefix = "mr_io_test_parallel_space_time"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"

    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        voxel_feature = np.random.rand(67,43,29,11,3)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestSpaceTimeMRI, self.mri.voxel_feature)
        print(block_dims)
        print(dims_mem)
        print(dims_mem_boundary)

    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestSpaceTimeMRI.filename_full)

        ## Read HDF5 from Fortran
        fort_command = "fortran/mr_io_test_parallel_reader_space_time %s  1> %s 2> %s" % \
                               (TestSpaceTimeMRI.filename_full,
                                TestSpaceTimeMRI.filename_out_rank % ("${PMI_RANK}"),
                                TestSpaceTimeMRI.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestSpaceTimeMRI.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        # Read HDF5 from Fortran
        #fort = sp.run(["mpiexec","-np", "%d" % (TestSpaceTimeMRI.mpi_proc), \
        #               "bash","-c","echo  \"spatial-mri\n 1 3 3\n %s \n \" 1> %s 2> %s" % \
        #                       (" ".join(9 * ["${PMI_RANK}"]), 
        #                        TestSpaceTimeMRI.filename_out_rank % ("${PMI_RANK}"), 
        #                        TestSpaceTimeMRI.filename_err_rank % ("${PMI_RANK}"))], 
        #                        stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)
        #import pdb; pdb.set_trace()

        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestSpaceTimeMRI, self.mri.voxel_feature)

        for mpi_rank in range(TestSpaceTimeMRI.mpi_proc):
            filename_out = TestSpaceTimeMRI.filename_out_rank % (mpi_rank)

            with open(filename_out,'r') as out:
                out_lines = out.readlines()
                
                for i in range(4):
                    fort_coord_dim = np.fromstring(out_lines[2*i+1][:-1], dtype=int, sep=' ')[0]
                    fort_coord_array = np.fromstring(out_lines[2*i+2][:-1], dtype=float, sep=' ')
                    self.assertTrue( (fort_coord_dim,) == fort_coord_array.shape )
                    self.assertTrue( np.allclose(fort_coord_array, self.mri.time if i == 0 else \
                                                                   self.mri.geometry[i-1], rtol=1e-14))
                                
                lo = 8                
                fort_full_dims_str = out_lines[lo+1][:-1]
                fort_offset_str = out_lines[lo+2][:-1] 
                fort_dims_str = out_lines[lo+3][:-1]
                fort_full_time_dim_str = out_lines[lo+4][:-1] 
                fort_time_offset_str = out_lines[lo+5][:-1] 
                fort_time_dim_str = out_lines[lo+6][:-1]
                fort_vector_dim_str = out_lines[lo+7][:-1]
                fort_array_str = out_lines[lo+8][:-1]
                #print(" *** " + filename_out + " *** ")
                #print(fort_dims_str)
                #print(fort_array_str)
                fort_full_dims = np.fromstring(fort_full_dims_str, dtype=int, sep=' ')
                fort_offset = np.fromstring(fort_offset_str, dtype=int, sep=' ')
                fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
                fort_full_time_dim = np.fromstring(fort_full_time_dim_str, dtype=int, sep=' ')[0]
                fort_time_offset = np.fromstring(fort_time_offset_str, dtype=int, sep=' ')[0]
                fort_time_dim = np.fromstring(fort_time_dim_str, dtype=int, sep=' ')[0]                
                fort_vector_dim = np.fromstring(fort_vector_dim_str, dtype=int, sep=' ')[0]                
                fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ')

                #import pdb; pdb.set_trace()
                fort_array = fort_array.reshape(np.flip([fort_vector_dim, fort_time_dim-fort_time_offset] + list(fort_dims))).transpose((2,1,0,3,4))
                
                block_id = (mpi_rank % block_dims[0], (mpi_rank // block_dims[0]) % block_dims[1], mpi_rank // (block_dims[0]* block_dims[1]) )
                hyperslab_offset = np.array([dims_mem[i]*block_id[i] for i in range(3)] + [0, 0])
                hyperslab_shape  = np.array([dims_mem[i] if block_id[i] + 1 < block_dims[i] else dims_mem_boundary[i] for i in range(3)] + [self.mri.voxel_feature.shape[i] for i in range(3,5)])
                #mpi_size = TestSpaceTimeMRI.mpi_proc
                #hyperslab_shape = self.mri.voxel_feature.shape[:2] + ( (self.mri.voxel_feature.shape[2]+mpi_size-1)//mpi_size if mpi_rank + 1 != mpi_size else (self.mri.voxel_feature.shape[2] % mpi_size if self.mri.voxel_feature.shape[2] % mpi_size !=0 else self.mri.voxel_feature.shape[2] / mpi_size ),)
                #hyperslab_offset = (0,0) + (mpi_rank*((self.mri.voxel_feature.shape[2]+mpi_size-1)//mpi_size),)
                self.assertTrue(np.array_equal(fort_full_dims, self.mri.voxel_feature.shape[:3]))
                self.assertTrue(fort_full_time_dim == self.mri.voxel_feature.shape[3])
                self.assertTrue(fort_vector_dim == self.mri.voxel_feature.shape[4])
                
                self.assertTrue(np.array_equal(list(fort_offset) + [fort_time_offset, 0], hyperslab_offset))
                self.assertTrue(np.array_equal(list(fort_dims) + [fort_time_dim-fort_time_offset, fort_vector_dim], hyperslab_shape))
                #import pdb; pdb.set_trace()
                self.assertTrue(np.allclose(fort_array, 
                                np.reshape(self.mri.voxel_feature[ \
                                           hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                           hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                           hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                           hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3],
                                           hyperslab_offset[4]:hyperslab_offset[4]+hyperslab_shape[4]], hyperslab_shape), rtol=1e-14))
    

    
    def tearDown(self):
        # Clean up file
        files = [TestSpaceTimeMRI.filename_full] + \
                [TestSpaceTimeMRI.filename_out_rank % (mpi_rank) for mpi_rank in range(TestSpaceTimeMRI.mpi_proc)] + \
                [TestSpaceTimeMRI.filename_err_rank % (mpi_rank) for mpi_rank in range(TestSpaceTimeMRI.mpi_proc)]
        for f in files:
            os.remove(f) 
            

class TestHPCPredictMRI(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 2**7

    # Filenames
    filename_prefix = "mr_io_test_parallel_hpc_predict"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"

    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        velocity_mean = np.random.rand(67,43,29,11,3)        
        velocity_cov = np.random.rand(67,43,29,11,3,5)        
        self.mri = HPCPredictMRI(geometry, time, velocity_mean, velocity_cov)
        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestHPCPredictMRI, self.mri.velocity_mean)
        print(block_dims)
        print(dims_mem)
        print(dims_mem_boundary)

    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestHPCPredictMRI.filename_full)

        ## Read HDF5 from Fortran
        fort_command = "fortran/mr_io_test_parallel_reader_hpc_predict %s  1> %s 2> %s" % \
                               (TestHPCPredictMRI.filename_full,
                                TestHPCPredictMRI.filename_out_rank % ("${PMI_RANK}"),
                                TestHPCPredictMRI.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestHPCPredictMRI.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        # Read HDF5 from Fortran
        #fort = sp.run(["mpiexec","-np", "%d" % (TestHPCPredictMRI.mpi_proc), \
        #               "bash","-c","echo  \"spatial-mri\n 1 3 3\n %s \n \" 1> %s 2> %s" % \
        #                       (" ".join(9 * ["${PMI_RANK}"]), 
        #                        TestHPCPredictMRI.filename_out_rank % ("${PMI_RANK}"), 
        #                        TestHPCPredictMRI.filename_err_rank % ("${PMI_RANK}"))], 
        #                        stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)
        #import pdb; pdb.set_trace()

        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestHPCPredictMRI, self.mri.velocity_mean)

        for mpi_rank in range(TestHPCPredictMRI.mpi_proc):
            filename_out = TestHPCPredictMRI.filename_out_rank % (mpi_rank)

            with open(filename_out,'r') as out:
                out_lines = out.readlines()
                
                for i in range(4):
                    fort_coord_dim = np.fromstring(out_lines[2*i+1][:-1], dtype=int, sep=' ')[0]
                    fort_coord_array = np.fromstring(out_lines[2*i+2][:-1], dtype=float, sep=' ')
                    self.assertTrue( (fort_coord_dim,) == fort_coord_array.shape )
                    self.assertTrue( np.allclose(fort_coord_array, self.mri.time if i == 0 else \
                                                                   self.mri.geometry[i-1], rtol=1e-14))
                                
                lo = 8                
                fort_full_dims_str = out_lines[lo+1][:-1]
                fort_offset_str = out_lines[lo+2][:-1] 
                fort_dims_str = out_lines[lo+3][:-1]
                fort_full_time_dim_str = out_lines[lo+4][:-1] 
                fort_time_offset_str = out_lines[lo+5][:-1] 
                fort_time_dim_str = out_lines[lo+6][:-1]
                fort_vector_dim_str = out_lines[lo+7][:-1]
                fort_array_str = out_lines[lo+8][:-1]
                #print(" *** " + filename_out + " *** ")
                #print(fort_dims_str)
                #print(fort_array_str)
                fort_full_dims = np.fromstring(fort_full_dims_str, dtype=int, sep=' ')
                fort_offset = np.fromstring(fort_offset_str, dtype=int, sep=' ')
                fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
                fort_full_time_dim = np.fromstring(fort_full_time_dim_str, dtype=int, sep=' ')[0]
                fort_time_offset = np.fromstring(fort_time_offset_str, dtype=int, sep=' ')[0]
                fort_time_dim = np.fromstring(fort_time_dim_str, dtype=int, sep=' ')[0]                
                fort_vector_dim = np.fromstring(fort_vector_dim_str, dtype=int, sep=' ')[0]                
                fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ')

                #import pdb; pdb.set_trace()
                fort_array = fort_array.reshape(np.flip([fort_vector_dim, fort_time_dim-fort_time_offset] + list(fort_dims))).transpose((2,1,0,3,4))
                
                block_id = (mpi_rank % block_dims[0], (mpi_rank // block_dims[0]) % block_dims[1], mpi_rank // (block_dims[0]* block_dims[1]) )
                hyperslab_offset = np.array([dims_mem[i]*block_id[i] for i in range(3)] + [0, 0])
                hyperslab_shape  = np.array([dims_mem[i] if block_id[i] + 1 < block_dims[i] else dims_mem_boundary[i] for i in range(3)] + [self.mri.velocity_mean.shape[i] for i in range(3,5)])
                #mpi_size = TestHPCPredictMRI.mpi_proc
                #hyperslab_shape = self.mri.velocity_mean.shape[:2] + ( (self.mri.velocity_mean.shape[2]+mpi_size-1)//mpi_size if mpi_rank + 1 != mpi_size else (self.mri.velocity_mean.shape[2] % mpi_size if self.mri.velocity_mean.shape[2] % mpi_size !=0 else self.mri.velocity_mean.shape[2] / mpi_size ),)
                #hyperslab_offset = (0,0) + (mpi_rank*((self.mri.velocity_mean.shape[2]+mpi_size-1)//mpi_size),)
                self.assertTrue(np.array_equal(fort_full_dims, self.mri.velocity_mean.shape[:3]))
                self.assertTrue(fort_full_time_dim == self.mri.velocity_mean.shape[3])
                self.assertTrue(fort_vector_dim == self.mri.velocity_mean.shape[4])
                
                self.assertTrue(np.array_equal(list(fort_offset) + [fort_time_offset, 0], hyperslab_offset))
                self.assertTrue(np.array_equal(list(fort_dims) + [fort_time_dim-fort_time_offset, fort_vector_dim], hyperslab_shape))
                #import pdb; pdb.set_trace()
                self.assertTrue(np.allclose(fort_array, 
                                np.reshape(self.mri.velocity_mean[ \
                                           hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                           hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                           hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                           hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3],
                                           hyperslab_offset[4]:hyperslab_offset[4]+hyperslab_shape[4]], hyperslab_shape), rtol=1e-14))
    
                # identical apart from lo -> 16 and velocity_mean -> velocity_cov
                lo = 16                
                fort_full_dims_str = out_lines[lo+1][:-1]
                fort_offset_str = out_lines[lo+2][:-1] 
                fort_dims_str = out_lines[lo+3][:-1]
                fort_full_time_dim_str = out_lines[lo+4][:-1] 
                fort_time_offset_str = out_lines[lo+5][:-1] 
                fort_time_dim_str = out_lines[lo+6][:-1]
                fort_vector_dim_str = out_lines[lo+7][:-1]
                fort_array_str = out_lines[lo+8][:-1]
                #print(" *** " + filename_out + " *** ")
                #print(fort_dims_str)
                #print(fort_array_str)
                fort_full_dims = np.fromstring(fort_full_dims_str, dtype=int, sep=' ')
                fort_offset = np.fromstring(fort_offset_str, dtype=int, sep=' ')
                fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')
                fort_full_time_dim = np.fromstring(fort_full_time_dim_str, dtype=int, sep=' ')[0]
                fort_time_offset = np.fromstring(fort_time_offset_str, dtype=int, sep=' ')[0]
                fort_time_dim = np.fromstring(fort_time_dim_str, dtype=int, sep=' ')[0]                
                fort_vector_dim = np.fromstring(fort_vector_dim_str, dtype=int, sep=' ')[:2]                
                fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ')

                #import pdb; pdb.set_trace()
                fort_array = fort_array.reshape(np.flip([*fort_vector_dim, fort_time_dim-fort_time_offset] + list(fort_dims))).transpose((2,1,0,3,5,4))
                
                block_id = (mpi_rank % block_dims[0], (mpi_rank // block_dims[0]) % block_dims[1], mpi_rank // (block_dims[0]* block_dims[1]) )
                hyperslab_offset = np.array([dims_mem[i]*block_id[i] for i in range(3)] + [0, 0, 0])
                hyperslab_shape  = np.array([dims_mem[i] if block_id[i] + 1 < block_dims[i] else dims_mem_boundary[i] for i in range(3)] + [self.mri.velocity_cov.shape[i] for i in range(3,6)])
                #mpi_size = TestHPCPredictMRI.mpi_proc
                #hyperslab_shape = self.mri.velocity_cov.shape[:2] + ( (self.mri.velocity_cov.shape[2]+mpi_size-1)//mpi_size if mpi_rank + 1 != mpi_size else (self.mri.velocity_cov.shape[2] % mpi_size if self.mri.velocity_cov.shape[2] % mpi_size !=0 else self.mri.velocity_cov.shape[2] / mpi_size ),)
                #hyperslab_offset = (0,0) + (mpi_rank*((self.mri.velocity_cov.shape[2]+mpi_size-1)//mpi_size),)
                self.assertTrue(np.array_equal(fort_full_dims, self.mri.velocity_cov.shape[:3]))
                self.assertTrue(fort_full_time_dim == self.mri.velocity_cov.shape[3])
                self.assertTrue( np.array_equal(fort_vector_dim, self.mri.velocity_cov.shape[4:6]) )
                
                self.assertTrue(np.array_equal(list(fort_offset) + [fort_time_offset, 0, 0], hyperslab_offset))
                self.assertTrue(np.array_equal(list(fort_dims) + [fort_time_dim-fort_time_offset, *fort_vector_dim], hyperslab_shape))

                self.assertTrue(np.allclose(fort_array, 
                                np.reshape(self.mri.velocity_cov[ \
                                           hyperslab_offset[0]:hyperslab_offset[0]+hyperslab_shape[0],
                                           hyperslab_offset[1]:hyperslab_offset[1]+hyperslab_shape[1],
                                           hyperslab_offset[2]:hyperslab_offset[2]+hyperslab_shape[2],
                                           hyperslab_offset[3]:hyperslab_offset[3]+hyperslab_shape[3],
                                           hyperslab_offset[4]:hyperslab_offset[4]+hyperslab_shape[4],
                                           hyperslab_offset[5]:hyperslab_offset[5]+hyperslab_shape[5]], hyperslab_shape), rtol=1e-14))
        

    
    def tearDown(self):
        # Clean up file
        files = [TestHPCPredictMRI.filename_full] + \
                [TestHPCPredictMRI.filename_out_rank % (mpi_rank) for mpi_rank in range(TestHPCPredictMRI.mpi_proc)] + \
                [TestHPCPredictMRI.filename_err_rank % (mpi_rank) for mpi_rank in range(TestHPCPredictMRI.mpi_proc)]
        for f in files:
            os.remove(f) 
            
if __name__ == '__main__':
    unittest.main()

