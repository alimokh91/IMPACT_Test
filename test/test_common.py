import numpy as np

def validate_group_name(test_inst, fort_group_name):
    test_inst.assertEqual(fort_group_name, type(test_inst).mri_group_name)        

# For bidirectional test cases only    
def validate_array(test_inst, array, out_array):
    test_inst.assertTrue(np.array_equal(out_array.shape, array.shape))        
    test_inst.assertTrue(np.allclose(out_array, array, rtol=1e-14))

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

def mpi_cart_rank(mpi_rank, mpi_cart_dims):
    return (  mpi_rank // (mpi_cart_dims[1] * mpi_cart_dims[2]),
             (mpi_rank // mpi_cart_dims[2]) % mpi_cart_dims[1],
              mpi_rank % mpi_cart_dims[2] )

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
    dims_mem_boundary = np.array([coord_shape % dims_mem[i] \
                                  if coord_shape % dims_mem[i] != 0 
                                  else dims_mem[i] 
                                  for i,coord_shape in enumerate(voxel_feature.shape[:3])]) 
    return block_dims, dims_mem, dims_mem_boundary #FIXME : discard the last two return values, only compute dims...

def spatial_hyperslab_dims_new(cls, voxel_feature): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
    #import pdb; pdb.set_trace()
    mpi_size = cls.mpi_proc
    dims_mem = np.array(voxel_feature.shape[:3])
    block_dims = np.array((1,1,1))
    while (mpi_size > 1):
        refinement_axis = np.argmax(dims_mem)
        block_dims[refinement_axis] *= 2
        dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
        mpi_size //= 2
    return block_dims

def spatial_hyperslab_loc(cls, mpi_rank, mpi_cart_dims, voxel_feature): 
    #block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(cls, voxel_feature)
    with_padding = hasattr(cls, "num_pad_vox_lhs")
    num_vox_ext = np.array([ voxel_feature.shape[i] + 
                            (cls.num_pad_vox_lhs[i] + cls.num_pad_vox_rhs[i] if with_padding else 0) for i in range(3)])
    dims_mem = np.array([(num_vox_ext[i] + mpi_cart_dims[i] - 1)//mpi_cart_dims[i] for i in range(3)])
    dims_mem_boundary = np.array([coord_shape % dims_mem[i] \
                                  if coord_shape % dims_mem[i] != 0 
                                  else dims_mem[i] 
                                  for i,coord_shape in enumerate(num_vox_ext)])
    
    mpi_rank_cart = mpi_cart_rank(mpi_rank, mpi_cart_dims)
    hyperslab_offset = np.array([dims_mem[i]*mpi_rank_cart[i] for i in range(3)] + [0]*(len(voxel_feature.shape)-3))
    hyperslab_shape = np.array([dims_mem[i] if mpi_rank_cart[i] + 1 < mpi_cart_dims[i] else dims_mem_boundary[i] for i in range(3)] + list(voxel_feature.shape[3:]))
    if with_padding:
        hyperslab_shape[:3] -= np.maximum(hyperslab_offset[:3]+hyperslab_shape[:3]-(num_vox_ext-cls.num_pad_vox_rhs), [0,0,0])
        hyperslab_shape[:3] -= np.maximum(cls.num_pad_vox_lhs - hyperslab_offset[:3], [0,0,0])
        if any(hyperslab_shape[:3] <= 0):
            hyperslab_shape[:3] = [0, 0, 0]
        hyperslab_offset[:3] = np.maximum(hyperslab_offset[:3] - cls.num_pad_vox_lhs, [0,0,0])
    
    return hyperslab_offset, hyperslab_shape
