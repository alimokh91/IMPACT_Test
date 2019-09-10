import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, FlowMRI, SegmentedFlowMRI
from test_common import validate_group_name, validate_spatial_fort_array, validate_spacetime_fort_array, validate_spacetime_scalar_fort_array, validate_spacetime_vector_fort_array, validate_spacetime_matrix_fort_array

def write_hdf5_read_in_fortran(test_inst):
    test_cls = type(test_inst)
 
    # Write HDF5 from Python
    test_inst.mri.write_hdf5(test_cls.filename_mri)
    
    # Read HDF5 from Fortran
    # Run test
    fort = sp.run(["bash", "-c", "fortran/test/%s %s  1> %s 2> %s" %
                  (test_cls.filename_exec, 
                   test_cls.filename_mri,
                   test_cls.filename_out, 
                   test_cls.filename_err)], 
                  stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    # Debug test
#     fort = sp.run(["xterm","-e","gdb", "--args", *("fortran/test/%s %s" %
#                   (test_cls.filename_exec, test_cls.filename_mri)).split(" ")], 
#                   stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    print("Shell command returned out/err:")
    print(fort.stdout.decode("utf-8"))
    print(fort.stderr.decode("utf-8"))
         
    with open(test_cls.filename_err, 'r') as err:
        print("Fortran command returned err:")
        print(err.read())

def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri, test_cls.filename_err, test_cls.filename_out]
    for f in files:
        os.remove(f) 

        
class TestSpatialMRI(unittest.TestCase):
  
    # Filenames
    filename_prefix = "mr_io_test_reader"
     
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"
 
    mri_group_name = "spatial-mri"
 
    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(83,63,93)
        self.mri = SpatialMRI(voxel_feature)
      
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self) 
             
        with open(type(self).filename_out,'r') as out:
            out_lines = [l.strip(' \n') for l in out.readlines()]
            validate_group_name(self, out_lines[0])
 
            validate_spatial_fort_array(self, self.mri.voxel_feature, out_lines[1:3])
         
    def tearDown(self):
        remove_test_files(self)
  
   
   
class TestSpaceTimeMRI(unittest.TestCase):
      
    # Filenames
    filename_prefix = "mr_io_test_reader_space_time"
      
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"
  
    mri_group_name = "space-time-mri"
     
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        voxel_feature = np.random.rand(23,13,19,17,21)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self) 
              
        with open(type(self).filename_out,'r') as out:
            out_lines = [l.strip(' \n') for l in out.readlines()]
            validate_group_name(self, out_lines[0])
             
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])
 
            validate_spacetime_vector_fort_array(self, self.mri.voxel_feature, out_lines[9:11])
         
    def tearDown(self):
        # Clean up file
        remove_test_files(self)

  
class TestFlowMRI(unittest.TestCase):
  
    # Filenames
    filename_prefix = "mr_io_test_reader_flow"
     
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"
 
    mri_group_name = "flow-mri"

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity = np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        self.mri = FlowMRI(geometry, time, intensity, velocity_mean, velocity_cov)
  
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self) 
  
        with open(type(self).filename_out,'r') as out:
            out_lines = [l.strip(' \n') for l in out.readlines()]
            validate_group_name(self, out_lines[0])
            
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])

            validate_spacetime_scalar_fort_array(self, self.mri.intensity, out_lines[9:11])  
            validate_spacetime_vector_fort_array(self, self.mri.velocity_mean, out_lines[11:13])  
            validate_spacetime_matrix_fort_array(self, self.mri.velocity_cov, out_lines[13:15])  
  
    def tearDown(self):
        # Clean up file
        remove_test_files(self)
 
  
class TestSegmentedFlowMRI(unittest.TestCase):
 
    # Filenames
    filename_prefix = "mr_io_test_reader_segmented_flow"
      
    filename_exec = filename_prefix
    filename_mri = filename_prefix + ".h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"
  
    mri_group_name = "segmented-flow-mri" 
  
    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity = np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        segmentation_prob = np.random.rand(23,13,19,17)
        self.mri = SegmentedFlowMRI(geometry, time, intensity, velocity_mean, 
                                          velocity_cov, segmentation_prob)
  
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_read_in_fortran(self) 
   
        with open(type(self).filename_out,'r') as out:
            out_lines = [l.strip(' \n') for l in out.readlines()]
            validate_group_name(self, out_lines[0])
             
            validate_spatial_fort_array(self, self.mri.geometry[0], out_lines[1:3])
            validate_spatial_fort_array(self, self.mri.geometry[1], out_lines[3:5])
            validate_spatial_fort_array(self, self.mri.geometry[2], out_lines[5:7])
            validate_spatial_fort_array(self, self.mri.time, out_lines[7:9])
 
            validate_spacetime_scalar_fort_array(self, self.mri.intensity, out_lines[9:11])  
            validate_spacetime_vector_fort_array(self, self.mri.velocity_mean, out_lines[11:13])  
            validate_spacetime_matrix_fort_array(self, self.mri.velocity_cov, out_lines[13:15])  
            validate_spacetime_scalar_fort_array(self, self.mri.segmentation_prob, out_lines[15:17])  
         
    def tearDown(self):
        # Clean up file
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

