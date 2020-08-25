import argparse
from test_common import fortran_strace, \
                        fortran_exec, \
                        filename_mri, \
                        filename_mri_in, \
                        filename_mri_out, \
                        fortran_out, \
                        fortran_err


def get_fortran_args(test_cls):

    if not test_cls.__name__.endswith('Bidirectional'):
        fort_args = " %s " % \
                   (filename_mri(test_cls))
    else:
        fort_args = "%s %s" % \
                    (filename_mri_in(test_cls),
                     filename_mri_out(test_cls))
    return fort_args


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate Fortran MPI test command.')
    parser.add_argument('--test', type=str, required=True,
                        help='Unit test (TestModule.TestClass)')
    parser.add_argument('--type', type=str, default="all",
                        help='Parameter type [all, strace]')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    module_name, test_cls_name = args.test.split('.')

    if not module_name.startswith('test_sequential'):
        raise ValueError("{} should be a class from test_sequential module".format(args.test))

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
        raise ValueError("type parameter must be one of [all,strace]")
