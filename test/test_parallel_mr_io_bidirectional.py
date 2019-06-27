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
        files = [TestSpatialMRIBidirectional.filename_in, \
                 TestSpatialMRIBidirectional.filename_out] + \
                [TestSpatialMRIBidirectional.filename_out_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRIBidirectional.mpi_proc)] + \
                [TestSpatialMRIBidirectional.filename_err_rank % (mpi_rank) for mpi_rank in range(TestSpatialMRIBidirectional.mpi_proc)]
        for f in files:
            os.remove(f) 


class TestSpaceTimeMRIBidirectional(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 2**7

    # Filenames
    filename_prefix = "mr_io_test_parallel_space_time"
    filename_in = filename_prefix + "_in.h5"
    filename_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"

    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        geometry = [np.random.rand(67), np.random.rand(43), np.random.rand(29)]
        voxel_feature = np.random.rand(67,43,29,11,3)
        self.mri = SpaceTimeMRI(geometry, voxel_feature)

    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestSpaceTimeMRIBidirectional.filename_in)

        ## Read HDF5 from Fortran
        fort_command = "fortran/mr_io_test_parallel_reader_writer_space_time %s %s 1> %s 2> %s" % \
                               (TestSpaceTimeMRIBidirectional.filename_in,
                                TestSpaceTimeMRIBidirectional.filename_out,
                                TestSpaceTimeMRIBidirectional.filename_out_rank % ("${PMI_RANK}"),
                                TestSpaceTimeMRIBidirectional.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestSpaceTimeMRIBidirectional.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        # Read HDF5 from Fortran
        #fort = sp.run(["mpiexec","-np", "%d" % (TestSpaceTimeMRIBidirectionalBidirectional.mpi_proc), \
        #               "bash","-c","echo  \"spatial-mri\n 1 3 3\n %s \n \" 1> %s 2> %s" % \
        #                       (" ".join(9 * ["${PMI_RANK}"]), 
        #                        TestSpaceTimeMRIBidirectionalBidirectional.filename_out_rank % ("${PMI_RANK}"), 
        #                        TestSpaceTimeMRIBidirectionalBidirectional.filename_err_rank % ("${PMI_RANK}"))], 
        #                        stdout=sp.PIPE, stderr=sp.PIPE, check=True)

        
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)
        #import pdb; pdb.set_trace()

        out_mri = SpaceTimeMRI.read_hdf5(TestSpaceTimeMRIBidirectional.filename_out)

        for i in range(3):
            self.assertTrue(np.allclose(out_mri.geometry[i], self.mri.geometry[i], rtol=1e-14))

#        self.assertTrue(np.array_equal(out_mri.voxel_feature.shape, self.mri.voxel_feature.shape))        
        self.assertTrue(np.allclose(out_mri.voxel_feature, self.mri.voxel_feature, rtol=1e-14))
    
    def tearDown(self):
        # Clean up file
        files = [TestSpaceTimeMRIBidirectional.filename_in, \
                 TestSpaceTimeMRIBidirectional.filename_out] + \
                [TestSpaceTimeMRIBidirectional.filename_out_rank % (mpi_rank) for mpi_rank in range(TestSpaceTimeMRIBidirectional.mpi_proc)] + \
                [TestSpaceTimeMRIBidirectional.filename_err_rank % (mpi_rank) for mpi_rank in range(TestSpaceTimeMRIBidirectional.mpi_proc)]
        for f in files:
            os.remove(f) 




if __name__ == '__main__':
    unittest.main()

