import os
import argparse
import numpy as np
from mr_io_domain_decomp import spatial_hyperslab_dims
from test_common import fortran_strace, \
                        fortran_exec, \
                        filename_mri, \
                        filename_mri_in, \
                        filename_mri_out, \
                        fortran_out, \
                        fortran_err
# from test_parallel_mr_io_container import TestFlowMRIPadded,\
#                                           TestSegmentedFlowMRIPadded
# from test_parallel_mr_io_bidirectional_container import TestFlowMRIPaddedBidirectional,\
#                                                         TestSegmentedFlowMRIPaddedBidirectional,\
#                                                         TestFlowMRIPaddedToSpaceTimeBidirectional,\
#                                                         TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional


def get_fortran_args(test_cls):
    mpi_cart_dims = spatial_hyperslab_dims(test_cls.mpi_proc, test_cls.spatial_feature_shape)
    fort_args = " {:d} {:d} {:d} ".format(*mpi_cart_dims)

    if not test_cls.__name__.endswith('Bidirectional'):
        fort_args += " %s " % \
                       (filename_mri(test_cls))
    else:
        fort_args += " %s %s " % \
                    (filename_mri_in(test_cls),
                     filename_mri_out(test_cls))

    # if test_cls in [TestFlowMRIPadded,
    #                 TestSegmentedFlowMRIPadded,
    #                 TestFlowMRIPaddedBidirectional,
    #                 TestSegmentedFlowMRIPaddedBidirectional,
    #                 TestFlowMRIPaddedToSpaceTimeBidirectional,
    #                 TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional]:
    if "Padded" in test_cls.__name__: # add LHS & RHS padding voxels
        test_cls.num_vox = test_cls.spatial_feature_shape

        num_pad_vox = [int(np.ceil(2. * test_cls.padding[i] * test_cls.num_vox[i])) for i in range(3)]

        num_ext_vox = [((test_cls.num_vox[i] + num_pad_vox[i] + mpi_cart_dims[i] - 1) // mpi_cart_dims[i]) \
                       * mpi_cart_dims[i] for i in range(3)]
        num_pad_vox_lhs = [(num_ext_vox[i] - test_cls.num_vox[i]) // 2 for i in range(3)]
        num_pad_vox_rhs = [(num_ext_vox[i] - test_cls.num_vox[i] + 1) // 2 for i in range(3)]
        num_vox_per_proc = [num_ext_vox[i] // mpi_cart_dims[i] for i in range(3)]

        fort_args += " %d %d %d %d %d %d " % (*num_pad_vox_lhs, *num_pad_vox_rhs)

    # if test_cls in [TestFlowMRIPaddedToSpaceTimeBidirectional,
    #                 TestSegmentedFlowMRIPaddedToSpaceTimeBidirectional]:
    if "ToSpaceTime" in test_cls.__name__: # add temporal and spatial refinements
        fort_args += " %d %d %d %d " % (test_cls.tr, *test_cls.sr)

    return fort_args


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate Fortran MPI test command.')
    parser.add_argument('--test', type=str, required=True,
                        help='Unit test (TestModule.TestClass)')
    parser.add_argument('--type', type=str, default="all",
                        help='Parameter type [all, exec, args, out, err, strace]')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    module_name, test_cls_name = args.test.split('.')

    if not module_name.startswith('test_parallel'):
        raise ValueError("{} should be a class from test_parallel module".format(args.test))

    import importlib
    module = importlib.import_module(module_name)
    test_cls = getattr(module, test_cls_name)

    cmd = fortran_exec(test_cls) + \
        " " + get_fortran_args(test_cls) + \
        " 1> " + fortran_out(test_cls) + \
        " 2> " + fortran_err(test_cls)

    if args.type == 'all':
        print(cmd)
    elif args.type == 'strace':
        print(fortran_strace(test_cls) + cmd)
    else:
        raise ValueError("type parameter must be one of [all, strace]")
