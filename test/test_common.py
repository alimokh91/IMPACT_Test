import numpy as np
from mr_io_domain_decomp import spatial_hyperslab_dims, spatial_hyperslab_loc, DomainPadding

def validate_group_name(test_inst, fort_group_name):
    test_inst.assertEqual(fort_group_name, type(test_inst).mri_group_name)        

# For bidirectional test cases only    
def validate_array(test_inst, array, out_array):
    test_inst.assertTrue(np.array_equal(out_array.shape, array.shape))        
    test_inst.assertTrue(np.allclose(out_array, array, rtol=1e-14))

# For arrays read from Fortran streams
def validate_spatial_fort_array(test_inst, array, out_lines):
    fort_dims = np.fromstring(out_lines[0], dtype=int, sep=' ') 
    fort_array= np.fromstring(out_lines[1], dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose()

    test_inst.assertTrue(np.array_equal(fort_dims, array.shape)) #[0])
    test_inst.assertTrue(np.allclose(fort_array, array, rtol=1e-14))


def validate_spacetime_fort_array(test_inst, array, out_lines, transpose_dims):
    fort_dims = np.fromstring(out_lines[0], dtype=int, sep=' ') 
    fort_array= np.fromstring(out_lines[1], dtype=float, sep=' ').reshape(np.flip(fort_dims)).transpose(transpose_dims)

    test_inst.assertTrue(np.array_equal(fort_dims[-3:], array.shape[:3]))
    test_inst.assertTrue(np.array_equal(fort_dims[-4], array.shape[3]))
    test_inst.assertTrue(np.array_equal(fort_dims[:-4], array.shape[4:]))
    test_inst.assertTrue(np.allclose(fort_array, array, rtol=1e-14))

def validate_spacetime_scalar_fort_array(test_inst, array, out_lines):
    validate_spacetime_fort_array(test_inst, array, out_lines, transpose_dims=(2,1,0,3))

def validate_spacetime_vector_fort_array(test_inst, array, out_lines):
    validate_spacetime_fort_array(test_inst, array, out_lines, transpose_dims=(2,1,0,3,4))

def validate_spacetime_matrix_fort_array(test_inst, array, out_lines):
    validate_spacetime_fort_array(test_inst, array, out_lines, transpose_dims=(2,1,0,3,5,4))

# Spatial domain decomposition layout
def spatial_hyperslab_dims_test(cls, voxel_feature): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
    return spatial_hyperslab_dims(mpi_size=cls.mpi_proc, voxel_feature_shape=voxel_feature.shape[:3])


def spatial_hyperslab_loc_test(cls, mpi_rank, mpi_cart_dims, voxel_feature):
    domain_padding = DomainPadding(pad_vox_lhs=cls.num_pad_vox_lhs, 
                                   pad_vox_rhs=cls.num_pad_vox_rhs) \
                     if hasattr(cls, "num_pad_vox_lhs") else None
    return spatial_hyperslab_loc(mpi_rank, mpi_cart_dims, voxel_feature.shape, domain_padding=domain_padding)
