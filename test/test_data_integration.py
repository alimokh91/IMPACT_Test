import unittest
import subprocess as sp
import os
import logging
import numpy as np
import shutil

from mr_io import FlowMRI, SpaceTimeMRI
from test_common import spatial_hyperslab_dims_test, validate_array, validate_file_path, validate_refined_mri_coordinates, validate_replicated_mri_vector_array
from test_impact_common import start_impact_fortran

def convert_expt_to_hdf5(test_inst):
    test_cls = type(test_inst)
    conversion_command = "python  mri_datasource/convert_bern_expt_to_hpc_predict.py --input %s --output %s 1> %s 2> %s" \
        % (test_cls.data_source_json_filename, 
           test_cls.filename_mri_in,
           test_cls.filename_out_rank % ("data_conversion"),
           test_cls.filename_err_rank % ("data_conversion"))
    
    conversion_run= sp.run( 
        ["bash","-c",conversion_command],
                stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    
    print("Bern experimental data to hpc-predict-io HDF5 conversion command:")
    print(conversion_command)
    
    print("Shell command returned out/err:")
    print(conversion_run.stdout.decode("utf-8"))
    print(conversion_run.stderr.decode("utf-8"))
    
    with open(test_cls.filename_out_rank % ("data_conversion"), 'r') as f:
        print("Bern experimental data to hpc-predict-io HDF5 conversion returned out:")
        print(f.read())
    
    with open(test_cls.filename_err_rank % ("data_conversion"), 'r') as f:
        print("Bern experimental data to hpc-predict-io HDF5 conversion returned err:")
        print(f.read())
        
       
def remove_test_files(test_inst):
    test_cls = type(test_inst)
    # Clean up file
    files = [test_cls.filename_mri_in, test_cls.config_output] + \
            [test_cls.filename_out_rank % ("config_writer"), \
             test_cls.filename_err_rank % ("config_writer") ] + \
            [test_cls.filename_out_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)] + \
            [test_cls.filename_err_rank % (mpi_rank) for mpi_rank in range(test_cls.mpi_proc)] + \
            test_cls.impact_noise # IMPACT files that are written without asking
    for f in files:
        os.remove(f) 
    

class TestBerneExptIntegration(unittest.TestCase): # FIXME: coordinates test...
     
    # number of Fortran MPI processes
    mpi_proc = 4 #was 7
     
    # Filenames
    filename_prefix = os.environ['FORTRAN_TEST_BINARY_PATH'] + "mr_io_test_impact_mri"
 
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out_rank = filename_prefix + "_%s.out"
    filename_err_rank = filename_prefix + "_%s.err"
    impact_noise = ["test_beta.txt", "test_multigrid_properties.txt"] # IMPACT files that are written without asking

    data_filename_prefix = "mri_datasource/bern_data_experiments_hpc_predict/"
    data_filename_mri = data_filename_prefix + "bern_experimental_dataset_flow_mri.h5"
 
    data_source_filename_prefix = "mri_datasource/bern_data_experiments_source/"
    data_source_json_filename = data_source_filename_prefix + "bern_exp_metadata.json"

    num_vox = (-1, -1, -1)   # should be divisible by domain decomposition computed in this test
    sr = [2, 1, 2]
    tr = 2
    
    config_template = 'python/config.txt.j2'
    config_output = 'config.txt'
 
    padding = (0.5, 0.4, 0.7)
 
    def setUp(self):
        test_cls = type(self)
        if not os.path.exists(test_cls.data_source_json_filename):
            raise RuntimeError("The JSON file %s describing the dataset does not exist. Aborting..." % test_cls.data_source_json_filename)
        
                 
    def test_communicator(self):
        test_cls = type(self)
                    
        # Either ...
        convert_expt_to_hdf5(self)
        # ...or...
#         shutil.copy2(test_cls.data_filename_mri, test_cls.filename_mri_in)

        self.mri = FlowMRI.read_hdf5(test_cls.filename_mri_in)
        
        test_cls.num_vox = tuple(self.mri.intensity.shape[:3])
         
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
         
        num_pad_vox = [int(np.ceil(2.*test_cls.padding[i]*test_cls.num_vox[i])) for i in range(3)]
         
        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] -1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]        
        
        # Make sure number of pressure grid points per process is odd (number of (extended refined) MRI voxels must be even)
        for i in range(3):
            if ( (test_cls.sr[i]*num_ext_vox[i] // self.mpi_cart_dims[i]) % 2 != 0): # This will never be triggered by the above criterion
                test_cls.sr[i] += 1
        
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i]-test_cls.num_vox[i])//2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i]-test_cls.num_vox[i]+1)//2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i]//self.mpi_cart_dims[i] for i in range(3)]
         
        print("MPI cartesian dims: {}\nMRI voxels: {}\nDesired padding voxels: {}\nExtended MRI voxels: {}\nComputed padding voxels: {} {}".format(\
              self.mpi_cart_dims, test_cls.num_vox, num_pad_vox, num_ext_vox, test_cls.num_pad_vox_lhs, test_cls.num_pad_vox_rhs))
        
        start_impact_fortran(self)
 
        out_mri = SpaceTimeMRI.read_hdf5(type(self).filename_mri_out)
 
        validate_refined_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.voxel_feature)
 
    def tearDown(self):
#         pass
        remove_test_files(self)
        os.remove(type(self).filename_mri_out) 

        
        
if __name__ == '__main__':
    unittest.main()

