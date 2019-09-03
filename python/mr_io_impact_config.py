from mr_io import HPCPredictMRI
from jinja2 import Environment, FileSystemLoader
import numpy as np
import argparse
import logging

def spatial_hyperslab_dims(mpi_size, voxel_feature_shape): 
    #import pdb; pdb.set_trace()
    dims_mem = np.array(voxel_feature_shape)
    block_dims = np.array((1,1,1))
    while (mpi_size > 1):
        refinement_axis = np.argmax(dims_mem)
        block_dims[refinement_axis] *= 2
        dims_mem[refinement_axis] = (dims_mem[refinement_axis] + 2 - 1)// 2
        mpi_size //= 2
    dims_mem_boundary = np.array([ coord_shape % dims_mem[i] if coord_shape % dims_mem[i] != 0 else dims_mem[i] for i,coord_shape in enumerate(voxel_feature_shape) ]) 
    return block_dims, dims_mem, dims_mem_boundary

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Generate IMPACT input configuration file from MRI HDF5-message.')
    parser.add_argument('--mri', type=str,
                        help='hpc-predict-io file containting MRI')
    parser.add_argument('--sr', type=int, nargs=3,
                    help='Number of spatial refinements along each dimension compared to MRI grid')
    parser.add_argument('--tr', type=int,
                    help='Number of temporal refinements compared to MRI grid')
    parser.add_argument('--config', type=str, default='config.txt.j2',
                    help='Jinja2 template for configuration file in IMPACT')
    parser.add_argument('--output', type=str,  default='config.txt',
                    help='Output configuration file for IMPACT')
    parser.add_argument('--np', type=int, help='Number of MPI processes')
    args = parser.parse_args()

    mri = HPCPredictMRI.read_hdf5(args.mri)
    
    block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(args.np, [len(axis) for axis in mri.geometry])
    
    # This will fail
    #assert((dims_mem == dims_mem_boundary).all())
    
    template_args = dict()
    for i in range(3):
#         if (len(mri.geometry[i]) % block_dims[i] != 0):
#             logging.error("Number of MRI-voxels (%d) along %d-th dimension not divisible by # processes along this direction (%d). " \
#                           "This will result in a pressure grid decomposition that is nonconforming with MRI-voxel decomposition."
#                          % (len(mri.geometry[i]), i, block_dims[i]))
#             raise ValueError("Number of MRI-voxels (%d) along %d-th dimension must be divisible by number of processes along this direction (%d)." \
#                          % (len(mri.geometry[i]), i, block_dims[i]))

#         # sr define number of pressure grid "voxels" (bounded by pressure grid coordinates) per MRI voxel. 
#         # Due to IMPACT's limitations this has to be adjusted
#         num_refined_mri_grid_voxels = args.sr[i]*len(mri.geometry[i])
#         if (num_refined_mri_grid_voxels % block_dims[i] != 0): # This will never be triggered by the above criterion
#             logging.warn("Number of \"pressure grid voxels\" (%d) along %d-th dimension not divisible by # processes along this direction (%d) - " \
#                          "will be rounded up in order for IMPACT to start correctly." \
#                          % (num_refined_mri_grid_voxels, i, block_dims[i]))
#         # compute rounded up number of local "pressure grid voxels" (does not change num_refined_mri_grid_voxels if divisible by block_dims[i])
#         N_i_minus_one = (num_refined_mri_grid_voxels+block_dims[i]-1)//block_dims[i]
#         if N_i_minus_one % 2 == 0: # valid number of pressure grid voxels (# pressure grid points - 1 )
#             num_pressure_grid_points = N_i_minus_one*block_dims[i]+1 # the +1 at the end is to convert number of voxels to number of grid points
#         else: # need to add one pressure grid point per block to make per-process number odd
#             num_pressure_grid_points = (N_i_minus_one+1)*block_dims[i]+1 # the +1 at the end is to convert number of voxels to number of grid points

        num_mri_voxels = len(mri.geometry[i])
        
        # Extend MRI voxel grid to allow for uniform decomposition among MPI processes
        num_ext_mri_voxels = ((num_mri_voxels + block_dims[i] -1) // block_dims[i]) * block_dims[i]
        # TODO: With extra added padding (relative to MRI voxel grid size)
        #num_ext_mri_voxels = (( np.ceil((padding+1.0)*num_mri_voxels) + block_dims[i] -1) // block_dims[i]) * block_dims[i]
                
        num_padding_voxels = num_ext_mri_voxels - num_mri_voxels 
        num_padding_voxels_lhs = num_padding_voxels//2
        num_padding_voxels_rhs = num_padding_voxels - num_padding_voxels_lhs
 
        # Make sure number of pressure grid points per process is odd (number of (extended refined) MRI voxels must be even)
        if ( (args.sr[i]*num_ext_mri_voxels // block_dims[i]) % 2 != 0): # This will never be triggered by the above criterion
            args.sr[i] += 1
        num_refined_ext_mri_voxels = args.sr[i]*num_ext_mri_voxels
        
#         if ( (num_refined_ext_mri_voxels // block_dims[i]) % 2 != 0): # This will never be triggered by the above criterion
#             num_refined_ext_mri_voxels += block_dims[i]    
#             logging.error("Number of pressure grid voxels (or points when counting the boundary as 0.5) per process (%d) along %d-th dimension has to be even as required by IMPACT." \
#                          % (num_refined_ext_mri_voxels // block_dims[i], i))
#             raise ValueError
        
        num_pressure_grid_points = num_refined_ext_mri_voxels+1 # the +1 at the end is to convert number of voxels to number of grid points
                
        template_args['M%d' % (i+1)]  = num_pressure_grid_points # Number of pressure grid points in each direction
        # second line converts MRI point grid extent to MRI voxel grid extent,
        # third line, which is commented out, would convert pressure voxel grid extent to pressure point grid extent
        mri_voxel_width = (mri.geometry[i][-1]-mri.geometry[i][0])/(len(mri.geometry[i])-1.)
        template_args['L%d' % (i+1)]  = mri_voxel_width*num_ext_mri_voxels
        template_args['NB%d' % (i+1)] = block_dims[i] # Number of processes along this dimension
        template_args['y%d_origin' % (i+1)] = -1.*(mri.geometry[i][0] - mri_voxel_width/2. - mri_voxel_width*num_padding_voxels_lhs)


        template_args['kalman_num_data_voxels_per_process_%d' % (i+1)]  = num_ext_mri_voxels//block_dims[i]
        template_args['kalman_num_padding_data_voxels_lhs_%d' % (i+1)]  = num_padding_voxels_lhs
        template_args['kalman_num_padding_data_voxels_rhs_%d' % (i+1)]  = num_padding_voxels_rhs
    

    template_args['time_start'] = mri.time[0]   # 0.
    template_args['time_end']   = mri.time[-1]  # 5000.
    
    template_args['n_timesteps'] = args.tr*len(mri.time) #'100000000'
    
    # Kalman filter variables (probably subject to change) # FIXME: These fields are probably all outdated
    template_args['dtime_out_scal'] = '0.0' # Delta time for saving data from DNS simulation 
                                            #for subsequent Kalman-filtered sim.
    template_args['dtime_out_kalm'] = '0.2' # Delta time for Kalman-filtered sim.
    template_args['vel_initcond_file_yes'] = 'F' # for DNS set to true, for Kalman-filtered sim. to false
    
    # Load config.txt termplate and instantiate variables
    env = Environment(loader=FileSystemLoader(searchpath="./"))
    template = env.get_template(args.config)
    with open(args.output,'w') as config_file:
        config_file.write(template.render(**template_args))
    
    
if __name__ == '__main__':
    main()
    
