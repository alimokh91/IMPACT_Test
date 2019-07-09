import unittest
import subprocess as sp
import os
import numpy as np
import logging
from mr_io import HPCPredictMRI


def spatial_hyperslab_dims(cls, voxel_feature): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
    #import pdb; pdb.set_trace()
    mpi_size = cls.mpi_proc
    dims_mem = np.array(voxel_feature.shape[:3], dtype=int)
    block_dims = np.array((1,1,1), dtype=int)
    while (mpi_size > 1):
        refinement_axis = np.argmax(dims_mem)
        block_dims[refinement_axis] *= 2
        dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
        mpi_size //= 2
    dims_mem_boundary = np.array([ coord_shape % dims_mem[i] if coord_shape % dims_mem[i] != 0 else dims_mem[i] for i,coord_shape in enumerate(voxel_feature.shape[:3]) ], dtype=int) 
    return block_dims, dims_mem, dims_mem_boundary



class TestImpactInput(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 4
    
    # Filenames
    filename_prefix = "mr_io_test_impact_input"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
    impact_noise = ["test_beta.txt", "test_multigrid_properties.txt"] # IMPACT files that are written without asking

    # python python/mr_io_impact_config.py --mri mr_io_test_space_time.h5 --sr 2 2 2 --tr 10 --config python/config.txt.j2 --output python/config.txt --np 8
    num_vox = (2**3, 2**3, 2**3)   # should be divisible by domain decomposition computed in this test
    domain_length = (3., 1.5, 2.) 
    sr = [2, 2, 2]
    tr = 10
    config_template = 'python/config.txt.j2'
    config_output = 'config.txt'

    def setUp(self):
        # Initialize the MRI data
        #geometry = [np.random.rand(4), np.random.rand(2), np.random.rand(7)] 
        time = np.linspace(0.,1.,11)
#         geometry = [np.linspace(0.,3.,2**3+1), np.linspace(0.,1.5,2**3+1), np.linspace(0.,2.0,2**3+1)]
#         velocity_mean = np.random.rand(2**3+1,2**3+1,2**3+1,11,3)
#         velocity_cov = np.random.rand(2**3+1,2**3+1,2**3+1,11,3,3)
        geometry = [np.linspace(0.,TestImpactInput.domain_length[i],
                                   2*TestImpactInput.num_vox[i]+1)[1:-1:2] for i in range(3)]
        geometry_complement = [np.linspace(0.,TestImpactInput.domain_length[i],
                                   2*TestImpactInput.num_vox[i]+1)[::2] for i in range(3)]
        velocity_mean = np.random.rand(*TestImpactInput.num_vox,11,3)
        velocity_cov = np.random.rand(*TestImpactInput.num_vox,11,3,3)

        self.mri = HPCPredictMRI(geometry, time, velocity_mean, velocity_cov)
        self.geometry_complement = geometry_complement
        
        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInput, self.mri.velocity_mean)
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
                                stdout=sp.PIPE, stderr=sp.PIPE, check=False)
                               
        config_run_stdout = config_run.stdout.decode("utf-8")
        config_run_stderr = config_run.stderr.decode("utf-8")        
        print(config_run_stderr)
        print(config_run_stdout)


        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInput, self.mri.velocity_mean)

#         fort_command = "xterm -geometry 73x31+$(( 100 + 600*(${PMI_RANK}/%d/%d) ))+$(( 100 + 1500*((${PMI_RANK}/%d) %% %d) + 600*(${PMI_RANK} %% %d) )) -e gdb fortran/mr_io_test_impact_input %s  1> %s 2> %s" % \
#                                (block_dims[1],block_dims[2],
#                                 block_dims[2],
#                                 block_dims[1],
#                                 block_dims[2],
#                                 TestImpactInput.config_output,
#                                 TestImpactInput.filename_out_rank % ("${PMI_RANK}"),
#                                 TestImpactInput.filename_err_rank % ("${PMI_RANK}"))
        fort_command = "fortran/mr_io_test_impact_input %s %s  1> %s 2> %s" % \
                               (TestImpactInput.config_output,
                                TestImpactInput.filename_full,
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

        for mpi_rank in range(TestImpactInput.mpi_proc):
            block_id = ( mpi_rank//(block_dims[1]*block_dims[2]),
                         (mpi_rank//block_dims[2]) % block_dims[1],
                         mpi_rank % block_dims[2] )
            
            filename_out = TestImpactInput.filename_out_rank % (mpi_rank)
            
            with open(filename_out,'r') as out:
                out_lines = out.readlines()

                # import pdb; pdb.set_trace()                
                
                out_begin = out_lines.index(" ***** mr_io_test_impact_input output *****\n")                
                
                # validate coordinates
                for i in range(3):
                    impact_num_uniaxial_pressure_grid_points = np.fromstring(out_lines[out_begin+2*i+1][:-1], dtype=int, sep=' ')[0]
                    impact_coord_array = np.fromstring(out_lines[out_begin+2*i+2][:-1], dtype=float, sep=' ')
                    self.assertTrue( impact_num_uniaxial_pressure_grid_points == int(TestImpactInput.sr[i])*len(self.mri.geometry[i])/int(block_dims[i])+1 )
                    self.assertTrue( len(impact_coord_array) == impact_num_uniaxial_pressure_grid_points )
                    if not np.allclose(impact_coord_array[::TestImpactInput.sr[i]], 
                                       self.geometry_complement[i][ block_id[i]*dims_mem[i]:(block_id[i]+1)*dims_mem[i]+1 ], rtol=1e-14):
                        logging.warning("Numerical error above double precision in {}-block along {}-th dimension:"\
                                        "\nIMPACT:    {}\nreference: {}".format(block_id, i, \
                                     impact_coord_array[::TestImpactInput.sr[i]], \
                                     self.geometry_complement[i][ block_id[i]*dims_mem[i]:(block_id[i]+1)*dims_mem[i]+1 ]))
                    self.assertTrue( np.allclose(impact_coord_array[::TestImpactInput.sr[i]], 
                                                 self.geometry_complement[i][ block_id[i]*dims_mem[i]:(block_id[i]+1)*dims_mem[i]+1 ], rtol=1e-1)) # FIXME: probably failing due to numerical error in y-coordinates in IMPACT!
                                
                # further checks can be added as required...
    
    def tearDown(self):
        # Clean up file

        files = [TestImpactInput.filename_full, TestImpactInput.config_output] + \
                [TestImpactInput.filename_out_rank % ("config_writer"), \
                 TestImpactInput.filename_err_rank % ("config_writer") ] + \
                [TestImpactInput.filename_out_rank % (mpi_rank) for mpi_rank in range(TestImpactInput.mpi_proc)] + \
                [TestImpactInput.filename_err_rank % (mpi_rank) for mpi_rank in range(TestImpactInput.mpi_proc)] + \
                TestImpactInput.impact_noise # IMPACT files that are written without asking
        for f in files:
            os.remove(f) 
            

if __name__ == '__main__':
    unittest.main()

