import argparse
import yaml
import importlib

# Parse data input and output directories
def parse_args():
    parser = argparse.ArgumentParser(description='Validate data interface specification')
    parser.add_argument('--app', type=str, required=True, 
                    help='Data interface of an individual application running in the pipeline')
    return parser.parse_args()

args = parse_args()

(data_interface_filename, app_name) = args.app.split(':')
with open(data_interface_filename) as data_interface_file:
    app_yml = yaml.load(data_interface_file, Loader=yaml.FullLoader)

def validate_app_interface(interface_yml, interface_name):
    print("Validating interface {}".format(interface_name))
    for name, python_type in interface_yml.items():
        module_name = python_type.split('.')[0]
        class_name = python_type[len(module_name)+1:]
        module = importlib.import_module(module_name) # TODO in a subprocess
        if hasattr(module, class_name):
            print("    {}: found class {}".format(name, getattr(module, class_name)))
        else:
            print("    {}: class {} not found!".format(name, python_type))
            raise ValueError("Class {} for {} does not exist in {}!".format(python_type, name, interface_name))

print("Validating data interface of {} app (@{})".format(app_yml['app'], data_interface_filename))
validate_app_interface(app_yml[app_name]['input'],  "{}:{}[{}]".format(app_yml['app'], app_name, 'input' ))
validate_app_interface(app_yml[app_name]['output'], "{}:{}[{}]".format(app_yml['app'], app_name, 'output')) 
print("Data interface of {} is valid!".format(args.app))

