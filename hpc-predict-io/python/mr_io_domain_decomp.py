import numpy as np


def mpi_cart_rank(mpi_rank, mpi_cart_dims):
    return (  mpi_rank // (mpi_cart_dims[1] * mpi_cart_dims[2]),
             (mpi_rank // mpi_cart_dims[2]) % mpi_cart_dims[1],
              mpi_rank % mpi_cart_dims[2] )


def spatial_hyperslab_dims(mpi_size, voxel_feature_shape): 
    
    # Compute prime factorization
    prime_factors = []
    current_factor = 2
    remainder = mpi_size
    while remainder >= current_factor:
        while remainder % current_factor == 0:
            remainder /= current_factor
            prime_factors.append(current_factor)
        current_factor += 1

    # Compute domain decomposition
    dims_mem = np.array(voxel_feature_shape)
    block_dims = np.array((1,1,1))
    for prime_factor in reversed(prime_factors):
        refinement_axis = np.argmax(dims_mem)
        block_dims[refinement_axis] *= prime_factor
        dims_mem[refinement_axis] = (dims_mem[refinement_axis] + prime_factor - 1)// prime_factor    
    
#     while (mpi_size > 1):
#         refinement_axis = np.argmax(dims_mem)
#         block_dims[refinement_axis] *= 2
#         dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
#         mpi_size //= 2

    return block_dims


class DomainPadding:
    def __init__(self, pad_vox_lhs, pad_vox_rhs):
        if (len(pad_vox_lhs) != 3):
            raise ValueError("List of number of lhs padding voxels to MRI grid must have length 3 (has %d instead)." \
                         % (len(pad_vox_lhs)))
        if (len(pad_vox_rhs) != 3):
            raise ValueError("List of number of rhs padding voxels to MRI grid must have length 3 (has %d instead)." \
                         % (len(pad_vox_rhs)))
        self.pad_vox_lhs = pad_vox_lhs
        self.pad_vox_rhs = pad_vox_rhs


def spatial_hyperslab_loc(mpi_rank, mpi_cart_dims, voxel_feature_shape, domain_padding):
    with_padding = isinstance(domain_padding, DomainPadding)
    
    # Compute extended domain and subdomain geometry
    num_vox_ext = np.array([ voxel_feature_shape[i] + 
                            (domain_padding.pad_vox_lhs[i] + domain_padding.pad_vox_rhs[i] if with_padding else 0) for i in range(3)])
    dims_mem = np.array([(num_vox_ext[i] + mpi_cart_dims[i] - 1)//mpi_cart_dims[i] for i in range(3)])
    dims_mem_boundary = np.array([coord_shape % dims_mem[i] \
                                  if coord_shape % dims_mem[i] != 0 
                                  else dims_mem[i] 
                                  for i,coord_shape in enumerate(num_vox_ext)])
    # Compute hyperslab location for left-most position of MRI
    mpi_rank_cart = mpi_cart_rank(mpi_rank, mpi_cart_dims)
    hyperslab_offset = np.array([dims_mem[i]*mpi_rank_cart[i] for i in range(3)] + [0]*(len(voxel_feature_shape)-3))
    hyperslab_shape = np.array([dims_mem[i] if mpi_rank_cart[i] + 1 < mpi_cart_dims[i] \
                                else dims_mem_boundary[i] for i in range(3)] + list(voxel_feature_shape[3:]))
    # Center MRI around padding
    if with_padding:
        hyperslab_shape[:3] -= np.maximum(hyperslab_offset[:3]+hyperslab_shape[:3]-(num_vox_ext-domain_padding.pad_vox_rhs), [0,0,0])
        hyperslab_shape[:3] -= np.maximum(domain_padding.pad_vox_lhs - hyperslab_offset[:3], [0,0,0])
        if any(hyperslab_shape[:3] <= 0):
            hyperslab_shape[:3] = [0, 0, 0]
        hyperslab_offset[:3] = np.maximum(hyperslab_offset[:3] - domain_padding.pad_vox_lhs, [0,0,0])
    
    return hyperslab_offset, hyperslab_shape