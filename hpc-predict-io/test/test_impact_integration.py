import unittest
import os
import logging
import numpy as np

from mr_io import SegmentedFlowMRI, SpaceTimeMRI
from mr_io_domain_decomp import mpi_cart_rank
from test_common import filename_mri_in, filename_mri_out, filename_out, filename_err, remove_test_files, \
                        spatial_hyperslab_dims_test, \
                        validate_array, \
                        validate_file_path, \
                        validate_refined_mri_coordinates, \
                        validate_replicated_mri_vector_array, \
                        read_out, print_out, print_err

    
def validate_impact_coordinates(test_inst):
    test_cls = type(test_inst)
    for mpi_rank in range(test_cls.mpi_proc):
        mpi_rank_cart = mpi_cart_rank(mpi_rank, test_inst.mpi_cart_dims)

        out_lines = read_out(filename_out(test_cls), mpi_rank)
        print_err(filename_err(test_cls), mpi_rank)

        out_begin = out_lines.index("***** mr_io_test_impact_input output *****")

        num_vox_per_proc = [(len(test_inst.geometry_complement[i])-1)//test_inst.mpi_cart_dims[i] for i in range(3)]
        # validate coordinates
        for i in range(3):
            impact_num_uniaxial_pressure_grid_points = np.fromstring(out_lines[out_begin+2*i+1], dtype=int, sep=' ')[0]
            impact_coord_array = np.fromstring(out_lines[out_begin+2*i+2], dtype=float, sep=' ')
            test_inst.assertTrue( impact_num_uniaxial_pressure_grid_points == test_cls.sr[i]*(len(test_inst.geometry_complement[i])-1)//int(test_inst.mpi_cart_dims[i])+1 )
            test_inst.assertTrue( len(impact_coord_array) == impact_num_uniaxial_pressure_grid_points )
            if not np.allclose(impact_coord_array[::test_cls.sr[i]],
                               test_inst.geometry_complement[i][mpi_rank_cart[i]*num_vox_per_proc[i]:(mpi_rank_cart[i]+1)*num_vox_per_proc[i]+1], rtol=1e-14):
                logging.warning("Numerical error above double precision in {}-block along {}-th dimension:"\
                                "\nIMPACT:    {}\nreference: {}".format(mpi_rank_cart, i, \
                                impact_coord_array[::test_cls.sr[i]], \
                                test_inst.geometry_complement[i][mpi_rank_cart[i]*num_vox_per_proc[i]:(mpi_rank_cart[i]+1)*num_vox_per_proc[i]+1]))

            validate_array(test_inst,
                           impact_coord_array[::test_cls.sr[i]],
                           test_inst.geometry_complement[i][mpi_rank_cart[i]*num_vox_per_proc[i]:(mpi_rank_cart[i]+1)*num_vox_per_proc[i]+1])

        # Validate IMPACT extended MRI grid parameters
        validate_array(test_inst, np.fromstring(out_lines[out_begin+7 ], dtype=int, sep=' '), np.array(test_cls.num_vox_per_proc))  # kalman_num_data_voxels_per_process
        validate_array(test_inst, np.fromstring(out_lines[out_begin+8], dtype=int, sep=' '), np.array(test_cls.num_pad_vox_lhs))  # kalman_num_padding_data_voxels_lhs
        validate_array(test_inst, np.fromstring(out_lines[out_begin+9], dtype=int, sep=' '), np.array(test_cls.num_pad_vox_rhs))  # kalman_num_padding_data_voxels_rhs

        out_metainfo_begin = out_begin + 10

        validate_file_path(test_inst, out_lines[out_metainfo_begin + 0].strip(), filename_mri_in(test_inst)) # kalman_mri_input_file_path
        validate_file_path(test_inst, out_lines[out_metainfo_begin + 1].strip(), filename_mri_out(test_inst)) # kalman_mri_output_file_path
        validate_array(test_inst, np.fromstring(out_lines[out_metainfo_begin + 2], dtype=int, sep=' '), np.array([test_inst.tr])) # kalman_num_time_refinements
        validate_array(test_inst, np.fromstring(out_lines[out_metainfo_begin + 3], dtype=int, sep=' '), np.array(test_inst.sr)) # kalman_num_spatial_refinements
        validate_array(test_inst, np.fromstring(out_lines[out_metainfo_begin + 4], dtype=float, sep=' '), np.array([test_inst.mri.time_heart_cycle_period])) # kalman_mri_input_attr_t_heart_cycle_period


class TestImpactInput(unittest.TestCase):
       
    fortran_exec = 'mr_io_test_impact_input'
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE']) # 2**2
       
    num_vox = (2**4, 2**5, 2**3)   # should be divisible by domain decomposition computed in this test
    domain_origin = (1., 4., 5.)
    domain_length = (3., 1.5, 2.) 
    sr = [2, 2, 2]
    tr = 10

    padding = (0., 0., 0.)
   
    def setUp(self):
        test_cls = type(self)
        # Initialize the MRI data
        time = np.linspace(0.,1.,11)
        time_heart_cycle_period = 1.
        intensity = np.random.rand(*TestImpactInput.num_vox,11)
        segmentation_prob = np.random.rand(*TestImpactInput.num_vox,11)
        velocity_mean = np.random.rand(*TestImpactInput.num_vox,11,3)
        velocity_cov = np.random.rand(*TestImpactInput.num_vox,11,3,3)
   
        geometry = [TestImpactInput.domain_origin[i] + \
                    np.linspace(0.,TestImpactInput.domain_length[i],
                                   2*TestImpactInput.num_vox[i]+1)[1:-1:2] for i in range(3)]
        self.geometry_complement = \
                   [TestImpactInput.domain_origin[i] + \
                    np.linspace(0.,TestImpactInput.domain_length[i],
                                   2*TestImpactInput.num_vox[i]+1)[::2] for i in range(3)]
   
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov, segmentation_prob)
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), self.mri.intensity)
   
        test_cls.num_pad_vox_lhs = [0]*3
        test_cls.num_pad_vox_rhs = [0]*3
        test_cls.num_vox_per_proc = [test_cls.num_vox[i]//self.mpi_cart_dims[i] for i in range(3)]
   
        print("MPI cartesian dims: {}\nMRI voxels: {}\nDesired padding voxels: {}\nExtended MRI voxels: {}\nComputed padding voxels: {} {}".format(\
              self.mpi_cart_dims, test_cls.num_vox, [0]*3, test_cls.num_vox, test_cls.num_pad_vox_lhs, test_cls.num_pad_vox_rhs))
                   
    def test_communicator(self):
        self.mri.write_hdf5(filename_mri_in(self))
        validate_impact_coordinates(self)
        # further checks can be added as required... (MRI reader testing not done here as already tested elsewhere)
          
    def tearDown(self):
        remove_test_files(self)
              
  
class TestImpactInputPadding(unittest.TestCase):

    fortran_exec = 'mr_io_test_impact_input'
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE'])  # FIXME: 7 !!

    num_vox = (95, 65, 55)   # should be divisible by domain decomposition computed in this test
    domain_origin = (1., 4., 5.)
    domain_length = (3., 1.5, 2.)
    sr = [2, 2, 2]
    tr = 10

    padding = (0.5, 0.4, 0.7)
  
    def setUp(self):
        test_cls = type(self)
          
        # Initialize the MRI data
        time = np.linspace(0.,1.,11)
        time_heart_cycle_period = 1.
        intensity = np.random.rand(*test_cls.num_vox,11)
        segmentation_prob = np.random.rand(*test_cls.num_vox,11)
        velocity_mean = np.random.rand(*test_cls.num_vox,11,3)
        velocity_cov = np.random.rand(*test_cls.num_vox,11,3,3)
  
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), intensity)
          
        num_pad_vox = [int(np.ceil(2.*test_cls.padding[i]*test_cls.num_vox[i])) for i in range(3)]
          
        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] -1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i]-test_cls.num_vox[i])//2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i]-test_cls.num_vox[i]+1)//2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i]//self.mpi_cart_dims[i] for i in range(3)]
          
        geometry = [test_cls.domain_origin[i] + \
                    np.linspace(0, test_cls.domain_length[i],
                                2*num_ext_vox[i]+1)[1+2*test_cls.num_pad_vox_lhs[i]:1+2*(test_cls.num_pad_vox_lhs[i]+test_cls.num_vox[i]):2] for i in range(3)]
        self.geometry_complement = \
                   [test_cls.domain_origin[i] + \
                               np.linspace(0, test_cls.domain_length[i],
                                           2*num_ext_vox[i]+1)[::2] for i in range(3)]
          
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov, segmentation_prob)
  
        print("MPI cartesian dims: {}\nMRI voxels: {}\nDesired padding voxels: {}\nExtended MRI voxels: {}\nComputed padding voxels: {} {}".format(\
              self.mpi_cart_dims, test_cls.num_vox, num_pad_vox, num_ext_vox, test_cls.num_pad_vox_lhs, test_cls.num_pad_vox_rhs))
                  
    def test_communicator(self):
        self.mri.write_hdf5(filename_mri_in(self))
        validate_impact_coordinates(self)
        # further checks can be added as required... (MRI reader testing not done here as already tested elsewhere)
  
    def tearDown(self):
        remove_test_files(self)
         
         
class TestImpactMRI(unittest.TestCase): # FIXME: coordinates test...

    fortran_exec = 'mr_io_test_impact_mri'
    # number of Fortran MPI processes
    mpi_proc = int(os.environ['HPC_PREDICT_IO_TEST_FORTRAN_MPI_SIZE'])  # 2**2

    num_vox = (95, 65, 55)   # should be divisible by domain decomposition computed in this test
    domain_origin = (1., 4., 5.)
    domain_length = (3., 1.5, 2.)
    sr = [2, 2, 2]
    tr = 2

    padding = (0.5, 0.4, 0.7)
 
    def setUp(self):
        test_cls = type(self)
         
        # Initialize the MRI data
        time = np.linspace(0.,1.,11)
        time_heart_cycle_period = 1.
        intensity = np.random.rand(*test_cls.num_vox,11)
        segmentation_prob = np.random.rand(*test_cls.num_vox,11)
        velocity_mean = np.random.rand(*test_cls.num_vox,11,3)
        velocity_cov = np.random.rand(*test_cls.num_vox,11,3,3)
 
        self.mpi_cart_dims = spatial_hyperslab_dims_test(type(self), intensity)
         
        num_pad_vox = [int(np.ceil(2.*test_cls.padding[i]*test_cls.num_vox[i])) for i in range(3)]
         
        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + self.mpi_cart_dims[i] -1) // self.mpi_cart_dims[i]) \
                       * self.mpi_cart_dims[i] for i in range(3)]
        test_cls.num_pad_vox_lhs = [(num_ext_vox[i]-test_cls.num_vox[i])//2 for i in range(3)]
        test_cls.num_pad_vox_rhs = [(num_ext_vox[i]-test_cls.num_vox[i]+1)//2 for i in range(3)]
        test_cls.num_vox_per_proc = [num_ext_vox[i]//self.mpi_cart_dims[i] for i in range(3)]
         
        geometry = [test_cls.domain_origin[i] + \
                    np.linspace(0, test_cls.domain_length[i],
                                2*num_ext_vox[i]+1)[1+2*test_cls.num_pad_vox_lhs[i]:1+2*(test_cls.num_pad_vox_lhs[i]+test_cls.num_vox[i]):2] for i in range(3)]
        self.geometry_complement = \
                   [test_cls.domain_origin[i] + \
                               np.linspace(0, test_cls.domain_length[i],
                                           2*num_ext_vox[i]+1)[::2] for i in range(3)]
         
        self.mri = SegmentedFlowMRI(geometry, time, time_heart_cycle_period, intensity, velocity_mean, velocity_cov, segmentation_prob)
 
        print("MPI cartesian dims: {}\nMRI voxels: {}\nDesired padding voxels: {}\nExtended MRI voxels: {}\nComputed padding voxels: {} {}".format(\
              self.mpi_cart_dims, test_cls.num_vox, num_pad_vox, num_ext_vox, test_cls.num_pad_vox_lhs, test_cls.num_pad_vox_rhs))
                 
    def test_communicator(self):
        self.mri.write_hdf5(filename_mri_in(self))

        for mpi_rank in range(type(self).mpi_proc):
            mpi_rank_cart = mpi_cart_rank(mpi_rank, self.mpi_cart_dims)
            print("OUT/ERR of Fortran cart rank {} ({}):".format(mpi_rank_cart, mpi_rank))
            print_out(filename_out(self), mpi_rank)
            print_err(filename_err(self), mpi_rank)

        out_mri = SpaceTimeMRI.read_hdf5(filename_mri_out(self))
 
        validate_refined_mri_coordinates(self, self.mri, out_mri)        
        validate_replicated_mri_vector_array(self, self.mri.velocity_mean, out_mri.vector_feature)
 
    def tearDown(self):
        remove_test_files(self)
        
        
if __name__ == '__main__':
    unittest.main()

