import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpaceTimeMRI


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



class TestImpactInput(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 2**4

    # Filenames
    filename_prefix = "mr_io_test_impact_input"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"

    # python python/mr_io_impact_config.py --mri mr_io_test_space_time.h5 --sr 2 2 2 --tr 10 --config python/config.txt.j2 --output python/config.txt --np 8
    sr = [2, 2, 2] 
    tr = 10
    config_template = 'python/config.txt.j2'
    config_output = 'python/config.txt'

    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.random.rand(11)
        geometry = [np.random.rand(2**6), np.random.rand(2**5), np.random.rand(2**5)]
        voxel_feature = np.random.rand(2**6,2**5,2**4,11,3)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)

        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInput, self.mri.voxel_feature)
        print(block_dims)
        print(dims_mem)
        print(dims_mem_boundary)

    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestImpactInput.filename_full)

        ## Write configuration file for impact
        # python python/mr_io_impact_config.py --mri mr_io_test_space_time.h5 --sr 2 2 2 --tr 10 --config python/config.txt.j2 --output python/config.txt --np 8
        config_command = "python python/mr_io_impact_config.py --mri %s --sr %d %d %d --tr %d --config %s --output %s --np %d  1> %s 2> %s" % \
                               (TestImpactInput.filename_full,
                                TestImpactInput.sr[0],
                                TestImpactInput.sr[1],
                                TestImpactInput.sr[2],
                                TestImpactInput.tr,
                                TestImpactInput.config_template,
                                TestImpactInput.config_output,
                                TestImpactInput.mpi_proc,
                                TestImpactInput.filename_out_rank % ("config_writer"),
                                TestImpactInput.filename_err_rank % ("config_writer"))
        print(config_command)
        config_run = sp.run(["bash","-c",config_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
                               
        config_run_stdout = config_run.stdout.decode("utf-8")
        config_run_stderr = config_run.stderr.decode("utf-8")        
        print(config_run_stderr)
        print(config_run_stdout)

        import pdb; pdb.set_trace()

        fort_command = "fortran/mr_io_test_impact_input %s  1> %s 2> %s" % \
                               (TestImpactInput.config_output,
                                TestImpactInput.filename_out_rank % ("${PMI_RANK}"),
                                TestImpactInput.filename_err_rank % ("${PMI_RANK}"))
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestImpactInput.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)

        #import pdb; pdb.set_trace()

        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInput, self.mri.voxel_feature)

        for mpi_rank in range(TestImpactInput.mpi_proc):
            filename_out = TestImpactInput.filename_out_rank % (mpi_rank)

            with open(filename_out,'r') as out:
                out_lines = out.readlines()
                
                import pdb; pdb.set_trace()
                
                # validate coordinates
                for i in range(3):
                    fort_coord_dim = np.fromstring(out_lines[2*i+1][:-1], dtype=int, sep=' ')[0]
                    fort_coord_array = np.fromstring(out_lines[2*i+2][:-1], dtype=float, sep=' ')
                    self.assertTrue( (fort_coord_dim,) == fort_coord_array.shape )
                    self.assertTrue( np.allclose(fort_coord_array, self.mri.geometry[i], rtol=1e-14))
                                
                # further checks can be reused from parallel_reader_space_time test...
    
    def tearDown(self):
        # Clean up file
        files = [TestImpactInput.filename_full, TestImpactInput.config_output] + \
                [TestImpactInput.filename_out_rank % ("config_writer"), \
                 TestImpactInput.filename_err_rank % ("config_writer") ] + \
                [TestImpactInput.filename_out_rank % (mpi_rank) for mpi_rank in range(TestImpactInput.mpi_proc)] + \
                [TestImpactInput.filename_err_rank % (mpi_rank) for mpi_rank in range(TestImpactInput.mpi_proc)]
        for f in files:
            os.remove(f) 
            

if __name__ == '__main__':
    unittest.main()

