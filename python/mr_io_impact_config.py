from mr_io import SpaceTimeMRI
from jinja2 import Environment, FileSystemLoader
import numpy as np
import argparse

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

    mri = SpaceTimeMRI.read_hdf5(args.mri)
    
    block_dims, dims_mem, dims_mem_boundary = spatial_hyperslab_dims(args.np, [len(axis) for axis in mri.geometry])
    
    assert((dims_mem == dims_mem_boundary).all())
    
    template_args = dict()
    for i in range(3):
        template_args['M%d' % (i+1)]  = args.sr[i]*len(mri.geometry[i]) # Number of pressure grid points in each direction
        template_args['L%d' % (i+1)]  = mri.geometry[i][-1]-mri.geometry[i][0] # position of last pressure grid point (first is at 0.)
        template_args['NB%d' % (i+1)] = block_dims[i] # Number of processes along each dimension
    
    template_args['time_start'] = mri.time[0] # 0.
    template_args['time_end']   = mri.time[-1]  # 5000.
    
    template_args['n_timesteps'] = args.tr*len(mri.time) #'100000000'
    
    # Kalman filter variables (probably subject to change) # FIXME TBD
    template_args['dtime_out_scal'] = '0.0' # Delta time for saving data from DNS simulation 
                                            #for subsequent Kalman-filtered sim.
    template_args['dtime_out_kalm'] = '0.2' # Delta time for Kalman-filtered sim.
    template_args['vel_initcond_file_yes'] = 'F' # for DNS set to true, for Kalman-filtered sim. to false
    
    env = Environment(loader=FileSystemLoader(searchpath="./"))
    template = env.get_template(args.config)
    with open(args.output,'w') as config_file:
        config_file.write(template.render(**template_args))
    
    
if __name__ == '__main__':
    main()
    