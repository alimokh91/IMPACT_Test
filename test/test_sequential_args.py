import os
import argparse
from test_common import filename_mri, \
                        filename_mri_in, \
                        filename_mri_out, \
                        filename_out, \
                        filename_err


def get_fortran_exec(test_cls):
    return test_cls.fortran_exec


def get_fortran_args(test_cls):

    if not test_cls.__name__.endswith('Bidirectional'):
        fort_args = " %s " % \
                   (filename_mri(test_cls))
    else:
        fort_args = "%s %s" % \
                    (filename_mri_in(test_cls),
                     filename_mri_out(test_cls))
    return fort_args


def get_fortran_out(test_cls):
    return filename_out(test_cls)


def get_fortran_err(test_cls):
    return filename_err(test_cls)


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate Fortran MPI test command.')
    parser.add_argument('--test', type=str, required=True,
                    help='Unit test (TestModule.TestClass)')
    parser.add_argument('--type', type=str, default="all",
                    help='Parameter type [all, exec, args, out, err]')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    module_name, test_cls_name = args.test.split('.')

    if not module_name.startswith('test_sequential'):
        raise ValueError("{} should be a class from test_sequential module".format(args.test))

    import importlib
    module = importlib.import_module(module_name)
    test_cls = getattr(module, test_cls_name)

    if args.type == "all":
        print(get_fortran_exec(test_cls) + \
              " " + get_fortran_args(test_cls) + \
              " 1> " + get_fortran_out(test_cls) + \
              " 2> " + get_fortran_err(test_cls))
    elif args.type == "exec":
        print(get_fortran_exec(test_cls))
    elif args.type == "args":
        print(get_fortran_args(test_cls))
    elif args.type == "out":
        print(get_fortran_out(test_cls))
    elif args.type == "err":
        print(get_fortran_err(test_cls))
    else:
        raise ValueError("type parameter must be one of [all, exec, args, out, err]")
