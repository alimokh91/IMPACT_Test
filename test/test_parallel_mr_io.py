import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI

class TestSpatialMRI(unittest.TestCase):
    
    # number of Fortran MPI processes
    mpi_proc = 2**4

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
        block_dims, dims_mem, dims_mem_boundary = self.hyperslab_dims()
        print(block_dims)
        print(dims_mem)
        print(dims_mem_boundary)
        
    def hyperslab_dims(self): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
        #import pdb; pdb.set_trace()
        mpi_size = TestSpatialMRI.mpi_proc
        dims_mem = np.array(self.mri.voxel_feature.shape)
        block_dims = np.array((1,1,1))
        while (mpi_size > 1):
            refinement_axis = np.argmax(dims_mem)
            block_dims[refinement_axis] *= 2
            dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
            mpi_size //= 2
        dims_mem_boundary = np.array([ coord_shape % dims_mem[i] if coord_shape % dims_mem[i] != 0 else dims_mem[i] for i,coord_shape in enumerate(self.mri.voxel_feature.shape) ]) 
        return block_dims, dims_mem, dims_mem_boundary

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

        block_dims, dims_mem, dims_mem_boundary = self.hyperslab_dims()

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



#class TestSpaceTimeMRI(unittest.TestCase):

    #def setUp(self):
        ## Initialize the MRI data
        #geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        #voxel_feature = np.random.rand(23,13,19,17,21)
        #self.mri = SpaceTimeMRI(geometry, voxel_feature)

    #def test_communicator(self):
        ## Write HDF5 from Python
        #self.mri.write_hdf5("mr_io_test_space_time.h5")

        ## Read HDF5 from Fortran
        #fort = sp.run(["fortran/mr_io_test_reader_space_time","mr_io_test_space_time.h5"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        ## Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        #fort_stdout = fort.stdout.decode("utf-8")
        #fort_stderr = fort.stderr.decode("utf-8")   
        #print(fort_stderr)

        #fort_group_name, fort_x_dim_str, fort_x_coord_str,fort_y_dim_str, fort_y_coord_str, fort_z_dim_str, fort_z_coord_str, fort_dims_str, fort_array_str, _ = fort_stdout.split('\n')

        #fort_x_dim = np.fromstring(fort_x_dim_str, dtype=int, sep=' ')
        #fort_y_dim = np.fromstring(fort_y_dim_str, dtype=int, sep=' ')
        #fort_z_dim = np.fromstring(fort_z_dim_str, dtype=int, sep=' ')

        #fort_x_coord = np.fromstring(fort_x_coord_str, dtype=float, sep=' ')
        #fort_y_coord = np.fromstring(fort_y_coord_str, dtype=float, sep=' ')
        #fort_z_coord = np.fromstring(fort_z_coord_str, dtype=float, sep=' ')
     
        ##print("Fortran x-coordinates: "); print(fort_x_coord)        
        ##print("Python x-coordinates:  "); print(self.mri.geometry[0])
        ##print("Fortran y-coordinates: "); print(fort_y_coord)   
        ##print("Python y-coordinates:  "); print(self.mri.geometry[1])
        ##print("Fortran z-coordinates: "); print(fort_z_coord)   
        ##print("Python z-coordinates:  "); print(self.mri.geometry[2]) 

        #self.assertEqual(fort_x_dim, fort_x_coord.shape[0])
        #self.assertEqual(fort_y_dim, fort_y_coord.shape[0])
        #self.assertEqual(fort_z_dim, fort_z_coord.shape[0])
 
        #self.assertTrue(np.allclose(fort_x_coord, self.mri.geometry[0], rtol=1e-14))
        #self.assertTrue(np.allclose(fort_y_coord, self.mri.geometry[1], rtol=1e-14))
        #self.assertTrue(np.allclose(fort_z_coord, self.mri.geometry[2], rtol=1e-14))

        #fort_dims = np.fromstring(fort_dims_str, dtype=int, sep=' ')

        #fort_array = np.fromstring(fort_array_str, dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose((2,1,0,3,4))
        ## self.assertTrue(np.array_equal(fort_dims,self.mri.voxel_feature.shape))

        #self.assertTrue(np.allclose(fort_array, self.mri.voxel_feature, rtol=1e-14))

    #def tearDown(self):
        ## Clean up file
        #os.remove("mr_io_test_space_time.h5")

if __name__ == '__main__':
    unittest.main()

