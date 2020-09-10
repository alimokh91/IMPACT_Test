import argparse
import yaml
import importlib

# Parse data input and output directories
def parse_args():
    parser = argparse.ArgumentParser(description='Validate pipeline from data interface specifications')
    parser.add_argument('--pipeline', type=str, required=True, nargs='+',
                    help='Data interfaces of individual applications running in sequence in the pipeline')
    return parser.parse_args()

args = parse_args()

data_interfaces = []
for (data_interface_filename, app_name) in [name.split(':') for name in args.pipeline]:
    with open(data_interface_filename) as data_interface_file:
        data_interfaces.append( (app_name, yaml.load(data_interface_file, Loader=yaml.FullLoader)) )

for writer, reader in zip(data_interfaces[:-1], data_interfaces[1:]):
    (writer_app_name, writer_yml) = writer
    (reader_app_name, reader_yml) = reader
   
    print("Validating {} -> {}".format(writer_yml['app'], reader_yml['app']))

    for name, python_type in reader_yml[reader_app_name]['input'].items():
        assert name in writer_yml[writer_app_name]['output']
        assert python_type == writer_yml[writer_app_name]['output'][name]

print("All internal pipeline interfaces validated & compatible!")

