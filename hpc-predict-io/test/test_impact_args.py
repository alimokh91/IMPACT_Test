import argparse
from test_common import fortran_strace, \
                        fortran_exec, \
                        filename_mri_in, \
                        filename_mri_out, \
                        filename_config, \
                        filename_out, \
                        filename_err, \
                        fortran_out, \
                        fortran_err


def config_writer_out(test_cls):
    return filename_out(test_cls) % ("config_writer")


def config_writer_err(test_cls):
    return filename_err(test_cls) % ("config_writer")


def get_config_args(test_cls):
    return " python3 -m mr_io_impact_config " + \
           " --input-mri " + filename_mri_in(test_cls) + \
           " --output-mri " + filename_mri_out(test_cls) + \
           " --sr %d %d %d --padding %f %f %f --tr %d " % (*test_cls.sr, *test_cls.padding, test_cls.tr) + \
           " --np %d " % test_cls.mpi_proc + \
           " --output " + filename_config(test_cls) + \
           " 1> " + config_writer_out(test_cls) + \
           " 2> " + config_writer_err(test_cls)


def get_impact_args(test_cls):
    return fortran_exec(test_cls) + \
        " " + filename_config(test_cls) + \
        " 1> " + fortran_out(test_cls) + \
        " 2> " + fortran_err(test_cls)


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate Fortran MPI test command.')
    parser.add_argument('--test', type=str, required=True,
                    help='Unit test (TestModule.TestClass)')
    parser.add_argument('--type', type=str, default="all",
                    help='Parameter type [all, strace]')
    parser.add_argument('--component', type=str,
                    help='Component to generate command for [config, impact]')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    module_name, test_cls_name = args.test.split('.')

    if not module_name.startswith('test_impact'):
        raise ValueError("{} should be a class from test_impact module".format(args.test))

    import importlib
    module = importlib.import_module(module_name)
    test_cls = getattr(module, test_cls_name)

    if args.component == 'config':
        cmd = get_config_args(test_cls)
    elif args.component == 'impact':
        cmd = get_impact_args(test_cls)
    else:
        raise ValueError("component parameter must be one of [config, impact]")

    if args.type == 'all':
        print(cmd)
    elif args.type == 'strace':
        print(fortran_strace(test_cls) + cmd)
    else:
        raise ValueError("type parameter must be one of [all, strace]")
