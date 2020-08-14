import os
import numpy as np
from mr_io_domain_decomp import spatial_hyperslab_dims,\
                                spatial_hyperslab_loc,\
                                DomainPadding
from mr_io_locking_utils import LockedFileReader
import inspect
import re


def get_test_filename(cls, filename_ending, relative_to_cache_folder=False):
    relative_path = '_'.join([cls.__module__, cls.__name__]) + filename_ending
    if not relative_to_cache_folder:
        return os.path.join(os.environ['CI_CACHE_FOLDER'], relative_path)
    else:
        return relative_path


def filename_mri(cls, relative_to_cache_folder=False):
    if not inspect.isclass(cls):
        cls = type(cls)
    return get_test_filename(cls, '.h5', relative_to_cache_folder)


def filename_mri_in(cls, relative_to_cache_folder=False):
    if not inspect.isclass(cls):
        cls = type(cls)
    return get_test_filename(cls, '_in.h5', relative_to_cache_folder)


def filename_mri_out(cls, relative_to_cache_folder=False):
    if not inspect.isclass(cls):
        cls = type(cls)
    return get_test_filename(cls, '_out.h5', relative_to_cache_folder)


def filename_out(cls, relative_to_cache_folder=False):
    if not inspect.isclass(cls):
        cls = type(cls)
    if cls.__module__.startswith('test_sequential'):
        return get_test_filename(cls, '.out', relative_to_cache_folder)
    elif cls.__module__.startswith('test_parallel'):
        return get_test_filename(cls, '_%s.out', relative_to_cache_folder)
    else:
        raise ValueError("{} should be a class from either test_sequential or test_parallel module".format(cls))


def filename_err(cls, relative_to_cache_folder=False):
    if not inspect.isclass(cls):
        cls = type(cls)
    if cls.__module__.startswith('test_sequential'):
        return get_test_filename(cls, '.err', relative_to_cache_folder)
    elif cls.__module__.startswith('test_parallel'):
        return get_test_filename(cls, '_%s.err', relative_to_cache_folder)
    else:
        raise ValueError("{} should be a class from either test_sequential or test_parallel module".format(cls))


def remove_test_files(cls):
    if not inspect.isclass(cls):
        cls = type(cls)
    if cls.__module__.startswith('test_sequential'):
        if not cls.__name__.endswith('Bidirectional'):
            files = [filename_mri(cls), filename_out(cls), filename_err(cls)]
        else:
            files = [filename_mri_in(cls), filename_mri_out(cls), filename_out(cls), filename_err(cls)]
    elif cls.__module__.startswith('test_parallel'):
        if not cls.__name__.endswith('Bidirectional'):
            files = [filename_mri(cls)] + \
                    [filename_out(cls) % (mpi_rank) for mpi_rank in range(cls.mpi_proc)] + \
                    [filename_err(cls) % (mpi_rank) for mpi_rank in range(cls.mpi_proc)]
        else:
            files = [filename_mri_in(cls), filename_mri_out(cls)] + \
                    [filename_out(cls) % (mpi_rank) for mpi_rank in range(cls.mpi_proc)] + \
                    [filename_err(cls) % (mpi_rank) for mpi_rank in range(cls.mpi_proc)]
    else:
        raise ValueError("{} should be a class from either test_sequential or test_parallel module".format(cls))

    # Clean up files
    for f in files:
        os.remove(f)


def mri_group_name(cls):
    if not inspect.isclass(cls):
        cls = type(cls)
    mri_type = re.match("Test(.*MRI).*", cls.__name__).groups()[0]
    return (re.sub(r'([A-Z])', r'-\1', mri_type[:-3])[1:] + '-MRI').lower()


def validate_group_name(test_inst, fort_group_name):
    test_inst.assertEqual(fort_group_name, mri_group_name(test_inst))


def validate_file_path(test_inst, file_path, out_file_path):
    test_inst.assertEqual(file_path, out_file_path)


# For bidirectional test cases only
def validate_array(test_inst, array, out_array):
    test_inst.assertTrue(np.array_equal(out_array.shape, array.shape))        
    test_inst.assertTrue(np.allclose(out_array, array, rtol=1e-14))


def validate_array_approx(test_inst, array, out_array):
    test_inst.assertTrue(np.array_equal(out_array.shape, array.shape))        
    test_inst.assertTrue(np.allclose(out_array, array, rtol=1e-12))


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


def validate_replicated_mri_coordinates(test_inst, in_mri, out_mri):
    test_cls = type(test_inst)

    for i in range(3):
        for j in range(test_cls.sr[i]):
            validate_array_approx(test_inst, in_mri.geometry[i], out_mri.geometry[i][j::test_cls.sr[i]])
    
    for j in range(test_cls.tr):
        validate_array_approx(test_inst, in_mri.time, out_mri.time[j::test_cls.tr])


def validate_refined_mri_coordinates(test_inst, in_mri, out_mri):
    test_cls = type(test_inst)

    def refined_edge_grid(grid, num_refinements):
        edge_width = (grid[-1] - grid[0])/float(len(grid)-1)
        begin = grid[0] - 0.5*edge_width
        end = grid[-1] + 0.5*edge_width
        return np.linspace(begin,end,num_refinements*len(grid)*2+1)[1::2]
    
    for i in range(3):
        validate_array_approx(test_inst, refined_edge_grid(in_mri.geometry[i], test_cls.sr[i]), out_mri.geometry[i])        
    
    for j in range(test_cls.tr):
        validate_array_approx(test_inst, in_mri.time, out_mri.time[j::test_cls.tr])


def validate_replicated_mri_vector_array(test_inst, in_array, out_array):
    test_cls = type(test_inst)
    for ix in range(test_cls.sr[0]):
        for iy in range(test_cls.sr[1]):
            for iz in range(test_cls.sr[2]):
                for it in range(test_cls.tr):
                    for iv in range(3):
                        # refined_cell_index = iv + 3*(iz + test_cls.sr[2]*(iy + test_cls.sr[1]*ix))
                        validate_array_approx(
                            test_inst, 
                            in_array[:,:,:,:,iv], # + refined_cell_index 
                            out_array[ix::test_cls.sr[0], 
                                      iy::test_cls.sr[1], 
                                      iz::test_cls.sr[2], 
                                      it::test_cls.tr,iv])


# Spatial domain decomposition layout
def spatial_hyperslab_dims_test(cls, voxel_feature): # NOTE: This is not the hyperslab_shape of particular block at (i,j,k), which has to account for boundaries
    return spatial_hyperslab_dims(mpi_size=cls.mpi_proc, voxel_feature_shape=voxel_feature.shape[:3])


def spatial_hyperslab_loc_test(cls, mpi_rank, mpi_cart_dims, voxel_feature):
    domain_padding = DomainPadding(pad_vox_lhs=cls.num_pad_vox_lhs, 
                                   pad_vox_rhs=cls.num_pad_vox_rhs) \
                     if hasattr(cls, "num_pad_vox_lhs") else None
    return spatial_hyperslab_loc(mpi_rank, mpi_cart_dims, voxel_feature.shape, domain_padding=domain_padding)


def read_out(filename_out, rank=None):
    if rank is not None:
        filename_out = filename_out % rank
    out_reader_lock = LockedFileReader(filename_out)
    out_reader_lock.open()
    out_content = out_reader_lock.file().read()
    out_lines = [l.strip(' \n') for l in out_content.split('\n')]
    print("STDOUT of Fortran rank {} ({}):".format(rank if rank is not None else "", filename_out))
    print('\n'.join(l if len(l) < 160 else l[:125]+" ... "+l[-30:] for l in out_lines))
    out_reader_lock.close()
    return out_lines


def print_err(filename_err, rank=None):
    if rank is not None:
        filename_err = filename_err % rank
    if os.stat(filename_err).st_size > 0:
        err_reader_lock = LockedFileReader(filename_err)
        err_reader_lock.open()
        err_content = err_reader_lock.file().read()
        if len(err_content) > 0:
            print("STDERR of Fortran rank {} ({}):".format(rank if rank is not None else "0", filename_err))
            print(err_content)
        err_reader_lock.close()
    else:
        print("STDERR of Fortran rank {} empty ({}).".format(rank if rank is not None else "0", filename_err))
