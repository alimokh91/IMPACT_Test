import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI

class TestSpatialMRIBidirectional(unittest.TestCase):
    
    # number of Fortran MPI processes
    mpi_proc = 2**7

    # Filenames
    filename_prefix = "mr_io_test_parallel"
    filename_in = filename_prefix + "_in.h5"
    filename_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"

    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(91,31,71) # splitting on the last dim in Fortran (must be the largest in size due to grid splitting!)
        #voxel_feature = np.reshape(np.array((9*[0]) + (9*[1])), (2,3,3))
        #print(voxel_feature)
        self.mri = SpatialMRI(voxel_feature)
        
    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestSpatialMRIBidirectional.filename_in)

        ## Read HDF5 from Fortran
        fort_command = "fortran/mr_io_test_parallel_reader_writer %s %s  1> %s 2> %s" % \
                               (TestSpatialMRIBidirectional.filename_in,
                                TestSpatialMRIBidirectional.filename_out,
                                TestSpatialMRIBidirectional.filename_out_rank % ("${PMI_RANK}"),
                                TestSpatialMRIBidirectional.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestSpatialMRIBidirectional.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
   
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)
        #import pdb; pdb.set_trace()

        out_mri = SpatialMRI.read_hdf5(TestSpatialMRIBidirectional.filename_out)

        self.assertTrue(np.array_equal(out_mri.voxel_feature.shape, self.mri.voxel_feature.shape))        
        self.assertTrue(np.allclose(out_mri.voxel_feature, self.mri.voxel_feature, rtol=1e-14))

    def tearDown(self):
        # Clean up file
        files = [TestSpatialMRIBidirectional.filename_in, TestSpatialMRIBidirectional.filename_out] + \
                [TestSpatialMRIBidirectional.filename_out_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRIBidirectional.mpi_proc)] + \
                [TestSpatialMRIBidirectional.filename_err_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRIBidirectional.mpi_proc)]
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

