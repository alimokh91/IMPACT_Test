import subprocess as sp

def start_impact_fortran(test_inst):
    test_cls = type(test_inst)

    ## Write configuration file for impact
    config_command = "python python/mr_io_impact_config.py --input-mri %s --output-mri %s --sr %d %d %d --padding %f %f %f --tr %d --config %s --output %s --np %d  1> %s 2> %s" % \
                           (test_cls.filename_mri_in,
                            test_cls.filename_mri_out,
                            test_cls.sr[0],
                            test_cls.sr[1],
                            test_cls.sr[2],
                            test_cls.padding[0],
                            test_cls.padding[1],
                            test_cls.padding[2],
                            test_cls.tr,
                            test_cls.config_template,
                            test_cls.config_output,
                            test_cls.mpi_proc,
                            test_cls.filename_out_rank % ("config_writer"),
                            test_cls.filename_err_rank % ("config_writer"))
                           
    print("IMPACT config generator command:")
    print(config_command)
    config_run = sp.run( 
         ["bash","-c",config_command],
                            stdout=sp.PIPE, stderr=sp.PIPE, check=True)


    print("Shell command returned out/err:")
    print(config_run.stdout.decode("utf-8"))
    print(config_run.stderr.decode("utf-8"))
         
    with open(test_cls.filename_out_rank % ("config_writer"), 'r') as f:
        print("IMPACT config generator returned out:")
        print(f.read())

    with open(test_cls.filename_err_rank % ("config_writer"), 'r') as f:
        print("IMPACT config generator returned err:")
        print(f.read())

#         fort_command = "xterm -geometry 73x31+$(( 100 + 600*(${PMI_RANK}/%d/%d) ))+$(( 100 + 1500*((${PMI_RANK}/%d) %% %d) + 600*(${PMI_RANK} %% %d) )) -e gdb fortran/test/mr_io_test_impact_input %s  1> %s 2> %s" % \
#                                (block_dims[1],block_dims[2],
#                                 block_dims[2],
#                                 block_dims[1],
#                                 block_dims[2],
#                                 test_cls.config_output,
#                                 test_cls.filename_out_rank % ("${PMI_RANK}"),
#                                 test_cls.filename_err_rank % ("${PMI_RANK}"))
    # FIXME: the command line parameters here are currently unused!
    fort_command = "%s %s 1> %s 2> %s" % \
                           (test_cls.filename_exec,
                            test_cls.config_output,
                            test_cls.filename_out_rank % ("${PMI_RANK}"),
                            test_cls.filename_err_rank % ("${PMI_RANK}"))
    
    print(fort_command)
    fort = sp.run(["mpiexec","-np", "%d" % (test_cls.mpi_proc), \
                   "bash","-c",fort_command],
                            stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    
    
    print("Shell command returned out/err:")
    print(fort.stdout.decode("utf-8"))
    print(fort.stderr.decode("utf-8"))
         
    for mpi_rank in range(test_cls.mpi_proc):
        with open(test_cls.filename_err_rank % (mpi_rank), 'r') as err:
            print("Fortran command for rank %d returned err:" % (mpi_rank))
            print(err.read())


