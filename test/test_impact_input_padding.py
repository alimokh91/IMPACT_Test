import unittest
import subprocess as sp
import os
import numpy as np
import logging
from mr_io import HPCPredictMRI
from test_common import spatial_hyperslab_dims


class TestImpactInputPadding(unittest.TestCase): # FIXME: coordinates test...
    
    # number of Fortran MPI processes
    mpi_proc = 4
    
    # Filenames
    filename_prefix = "mr_io_test_impact_input"
    filename_full = filename_prefix + ".h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
    impact_noise = ["test_beta.txt", "test_multigrid_properties.txt"] # IMPACT files that are written without asking

    # python python/mr_io_impact_config.py --mri mr_io_test_space_time.h5 --sr 2 2 2 --tr 10 --config python/config.txt.j2 --output python/config.txt --np 8
    num_vox = (7, 13, 11)   # should be divisible by domain decomposition computed in this test
    domain_origin = (1., 4., 5.)
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

        velocity_mean = np.random.rand(*TestImpactInputPadding.num_vox,11,3)
        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInputPadding, velocity_mean)
        
        num_ext_vox = [((TestImpactInputPadding.num_vox[i] + block_dims[i] -1) // block_dims[i]) \
                       * block_dims[i] for i in range(3)]
        num_pad_vox_lhs = [(num_ext_vox[i]-TestImpactInputPadding.num_vox[i])//2 for i in range(3)]
        self.num_vox_per_proc = [num_ext_vox[i]//block_dims[i] for i in range(3)]
        
        geometry = [TestImpactInputPadding.domain_origin[i] + \
                    np.linspace(0, TestImpactInputPadding.domain_length[i],
                                2*num_ext_vox[i]+1)[1+2*num_pad_vox_lhs[i]:1+2*(num_pad_vox_lhs[i]+TestImpactInputPadding.num_vox[i]):2] for i in range(3)]
        geometry_complement = [TestImpactInputPadding.domain_origin[i] + \
                               np.linspace(0, TestImpactInputPadding.domain_length[i],
                                           2*num_ext_vox[i]+1)[::2] for i in range(3)]
        intensity = np.random.rand(*TestImpactInputPadding.num_vox,11)
        velocity_cov = np.random.rand(*TestImpactInputPadding.num_vox,11,3,3)

        self.mri = HPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov)
        self.geometry_complement = geometry_complement
        
        print(block_dims)
        print(dims_mem)
        print(dims_mem_boundary)

    def test_communicator(self):

        # Write HDF5 from Python
        self.mri.write_hdf5(TestImpactInputPadding.filename_full)

        ## Write configuration file for impact
        # python python/mr_io_impact_config.py --mri mr_io_test_space_time.h5 --sr 2 2 2 --tr 10 --config python/config.txt.j2 --output python/config.txt --np 8
        config_command = "python python/mr_io_impact_config.py --mri %s --sr %d %d %d --tr %d --config %s --output %s --np %d  1> %s 2> %s" % \
                               (TestImpactInputPadding.filename_full,
                                TestImpactInputPadding.sr[0],
                                TestImpactInputPadding.sr[1],
                                TestImpactInputPadding.sr[2],
                                TestImpactInputPadding.tr,
                                TestImpactInputPadding.config_template,
                                TestImpactInputPadding.config_output,
                                TestImpactInputPadding.mpi_proc,
                                TestImpactInputPadding.filename_out_rank % ("config_writer"),
                                TestImpactInputPadding.filename_err_rank % ("config_writer"))
        print(config_command)
        config_run = sp.run(["bash","-c",config_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=False)
                               
        config_run_stdout = config_run.stdout.decode("utf-8")
        config_run_stderr = config_run.stderr.decode("utf-8")        
        print(config_run_stderr)
        print(config_run_stdout)


        block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(TestImpactInputPadding, self.mri.velocity_mean)

#         fort_command = "xterm -geometry 73x31+$(( 100 + 600*(${PMI_RANK}/%d/%d) ))+$(( 100 + 1500*((${PMI_RANK}/%d) %% %d) + 600*(${PMI_RANK} %% %d) )) -e gdb fortran/mr_io_test_impact_input %s  1> %s 2> %s" % \
#                                (block_dims[1],block_dims[2],
#                                 block_dims[2],
#                                 block_dims[1],
#                                 block_dims[2],
#                                 TestImpactInputPadding.config_output,
#                                 TestImpactInputPadding.filename_out_rank % ("${PMI_RANK}"),
#                                 TestImpactInputPadding.filename_err_rank % ("${PMI_RANK}"))
        fort_command = "fortran/mr_io_test_impact_input %s %s  1> %s 2> %s" % \
                               (TestImpactInputPadding.config_output,
                                TestImpactInputPadding.filename_full,
                                TestImpactInputPadding.filename_out_rank % ("${PMI_RANK}"),
                                TestImpactInputPadding.filename_err_rank % ("${PMI_RANK}"))
        
        print(fort_command)
        fort = sp.run(["mpiexec","-np", "%d" % (TestImpactInputPadding.mpi_proc), \
                       "bash","-c",fort_command],
                                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
        
        
        # Check for equality. NOTE: Possibly regexpr for parsing would make this more elegant
        fort_stdout = fort.stdout.decode("utf-8")
        fort_stderr = fort.stderr.decode("utf-8")        
        print(fort_stderr)
        print(fort_stdout)

        #import pdb; pdb.set_trace()

        for mpi_rank in range(TestImpactInputPadding.mpi_proc):
            block_id = ( mpi_rank//(block_dims[1]*block_dims[2]),
                         (mpi_rank//block_dims[2]) % block_dims[1],
                         mpi_rank % block_dims[2] )
            
            filename_out = TestImpactInputPadding.filename_out_rank % (mpi_rank)
            
            with open(filename_out,'r') as out:
                out_lines = out.readlines()

                # import pdb; pdb.set_trace()                
                
                out_begin = out_lines.index(" ***** mr_io_test_impact_input output *****\n")                
                
                # validate coordinates
                for i in range(3):
                    impact_num_uniaxial_pressure_grid_points = np.fromstring(out_lines[out_begin+2*i+1][:-1], dtype=int, sep=' ')[0]
                    impact_coord_array = np.fromstring(out_lines[out_begin+2*i+2][:-1], dtype=float, sep=' ')
                    self.assertTrue( impact_num_uniaxial_pressure_grid_points == TestImpactInputPadding.sr[i]*(len(self.geometry_complement[i])-1)/int(block_dims[i])+1)
                    self.assertTrue( len(impact_coord_array) == impact_num_uniaxial_pressure_grid_points )
                    if not np.allclose(impact_coord_array[::TestImpactInputPadding.sr[i]], 
                                       self.geometry_complement[i][ block_id[i]*self.num_vox_per_proc[i]:(block_id[i]+1)*self.num_vox_per_proc[i]+1 ], rtol=1e-14):
                        logging.warning("Numerical error above double precision in {}-block along {}-th dimension:"\
                                        "\nIMPACT:    {}\nreference: {}".format(block_id, i, \
                                     impact_coord_array[::TestImpactInputPadding.sr[i]], \
                                     self.geometry_complement[i][ block_id[i]*self.num_vox_per_proc[i]:(block_id[i]+1)*self.num_vox_per_proc[i]+1 ]))
                    self.assertTrue( np.allclose(impact_coord_array[::TestImpactInputPadding.sr[i]], 
                                             self.geometry_complement[i][ block_id[i]*self.num_vox_per_proc[i]:(block_id[i]+1)*self.num_vox_per_proc[i]+1 ], rtol=1e-14)) # FIXME: probably failing due to numerical error in y-coordinates in IMPACT!
                                
                # further checks can be added as required...
    
    def tearDown(self):
        # Clean up file

        files = [TestImpactInputPadding.filename_full, TestImpactInputPadding.config_output] + \
                [TestImpactInputPadding.filename_out_rank % ("config_writer"), \
                 TestImpactInputPadding.filename_err_rank % ("config_writer") ] + \
                [TestImpactInputPadding.filename_out_rank % (mpi_rank) for mpi_rank in range(TestImpactInputPadding.mpi_proc)] + \
                [TestImpactInputPadding.filename_err_rank % (mpi_rank) for mpi_rank in range(TestImpactInputPadding.mpi_proc)] + \
                 TestImpactInputPadding.impact_noise # IMPACT files that are written without asking
        for f in files:
            os.remove(f) 
            

if __name__ == '__main__':
    unittest.main()

