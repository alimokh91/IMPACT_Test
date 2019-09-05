import unittest
import subprocess as sp
import os
import numpy as np
from mr_io import SpatialMRI, SpaceTimeMRI, HPCPredictMRI, SegmentedHPCPredictMRI


def write_hdf5_exec_fortran(test_inst):
    test_cls = type(test_inst)
 
    # Write HDF5 from Python
    test_inst.mri.write_hdf5(test_cls.filename_mri_in)
    
    # Read HDF5 from Fortran
    # Run test
    fort = sp.run(["bash", "-c", "fortran/%s %s %s 1> %s 2> %s" %
                  (test_cls.filename_exec, 
                   test_cls.filename_mri_in, 
                   test_cls.filename_mri_out,
                   test_cls.filename_out, 
                   test_cls.filename_err)], 
                  stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    # Debug test
#     fort = sp.run(["xterm","-e","gdb", "--args", *("fortran/%s %s %s" %
#                   (test_cls.filename_exec, test_cls.filename_mri_in, test_cls.filename_mri_out)).split(" ")], 
#                   stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    print("Shell command returned out/err:")
    print(fort.stdout.decode("utf-8"))
    print(fort.stderr.decode("utf-8"))
         
    with open(test_cls.filename_err, 'r') as err:
        print("Fortran command returned err:")
        print(err.read())

def validate_group_name(test_inst, fort_group_name):
    test_inst.assertEqual(fort_group_name, type(test_inst).mri_group_name)        
    
def validate_array(test_inst, array, out_array):
    test_inst.assertTrue(np.array_equal(out_array.shape, array.shape))        
    test_inst.assertTrue(np.allclose(out_array, array, rtol=1e-14))

def remove_test_files(test_inst):
    # Clean up file
    test_cls = type(test_inst)
    files = [test_cls.filename_mri_in, test_cls.filename_mri_out, test_cls.filename_err, test_cls.filename_out]
    for f in files:
        os.remove(f) 


class TestSpatialMRIBidirectional(unittest.TestCase):

    # Filenames
    filename_prefix = "mr_io_test_reader_writer"
    
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"

    mri_group_name = "spatial-mri"

    def setUp(self):
        # Initialize the MRI data
        voxel_feature = np.random.rand(71,61,97)
        self.mri = SpatialMRI(voxel_feature)
    
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
            
        out_mri = SpatialMRI.read_hdf5(type(self).filename_mri_out)
        
        validate_group_name(self, out_mri.group_name)
        validate_array(self, self.mri.voxel_feature, out_mri.voxel_feature)        
        
    def tearDown(self):
        remove_test_files(self)


class TestSpaceTimeMRIBidirectional(unittest.TestCase):
 
    # Filenames
    filename_prefix = "mr_io_test_reader_writer_space_time"
    
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"

    mri_group_name = "space-time-mri"

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(13)
        geometry = [np.random.rand(26), np.random.rand(11), np.random.rand(7)]
        voxel_feature = np.random.rand(26,11,7,13,19)
        self.mri = SpaceTimeMRI(geometry, time, voxel_feature)
 
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
            
        out_mri = SpaceTimeMRI.read_hdf5(type(self).filename_mri_out)

        validate_group_name(self, out_mri.group_name)
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.voxel_feature, out_mri.voxel_feature)

    def tearDown(self):
        remove_test_files(self)
 
 
class TestHPCPredictMRIBidirectional(unittest.TestCase):
 
    # Filenames
    filename_prefix = "mr_io_test_reader_writer_hpc_predict"
    
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"

    mri_group_name = "hpc-predict-mri"

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity= np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        self.mri = HPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov)
 
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
            
        out_mri = HPCPredictMRI.read_hdf5(type(self).filename_mri_out)

        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
 
    def tearDown(self):
        remove_test_files(self)
         
 
class TestSegmentedHPCPredictMRIBidirectional(unittest.TestCase):
 
    # Filenames
    filename_prefix = "mr_io_test_reader_writer_segmented_hpc_predict"
    
    filename_exec = filename_prefix
    filename_mri_in = filename_prefix + "_in.h5"
    filename_mri_out = filename_prefix + "_out.h5"
    filename_out = filename_prefix + ".out"
    filename_err = filename_prefix + ".err"

    mri_group_name = "segmented-hpc-predict-mri"

    def setUp(self):
        # Initialize the MRI data
        time = np.random.rand(17)
        geometry = [np.random.rand(23), np.random.rand(13), np.random.rand(19)]
        intensity= np.random.rand(23,13,19,17)
        velocity_mean = np.random.rand(23,13,19,17,3)
        velocity_cov = np.random.rand(23,13,19,17,3,5)
        segmentation_prob = np.random.rand(23,13,19,17)
        self.mri = SegmentedHPCPredictMRI(geometry, time, intensity, velocity_mean, velocity_cov, segmentation_prob)
 
    def test_communicator(self):
        # Write HDF5 from Python and read HDF5 from Fortran   
        write_hdf5_exec_fortran(self) 
            
        out_mri = SegmentedHPCPredictMRI.read_hdf5(type(self).filename_mri_out)
     
        for i in range(3):
            validate_array(self, self.mri.geometry[i], out_mri.geometry[i])        
        validate_array(self, self.mri.time, out_mri.time)
        validate_array(self, self.mri.intensity, out_mri.intensity)
        validate_array(self, self.mri.velocity_mean, out_mri.velocity_mean)
        validate_array(self, self.mri.velocity_cov, out_mri.velocity_cov)
        validate_array(self, self.mri.segmentation_prob, out_mri.segmentation_prob)
 
    def tearDown(self):
        remove_test_files(self)


if __name__ == '__main__':
    unittest.main()

