import argparse
import yaml
import importlib

# Parse data input and output directories
def parse_args():
    parser = argparse.ArgumentParser(description='Validate pipeline from data interface specifications')
    parser.add_argument('--pipeline', type=str, required=True, nargs='+',
                    help='Data interfaces of individual applications running in sequence in the pipeline')
    parser.add_argument('--diagram', type=str, 
                    help='Jinja2 template of pipeline diagram to render')
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

if not args.diagram is None:
    import os
    from jinja2 import FileSystemLoader, Environment    

    assert os.path.exists(args.diagram) 
    assert args.diagram.endswith('.j2') 
    assert not os.path.exists(args.diagram[:-3])

    print("Writing rendered diagram template to {}".format(args.diagram[:-3]))

    template_args = dict()
    for (interface_filename, (interface_app_name, interface_yml)) in zip([name.split(':')[0] for name in args.pipeline], data_interfaces):
        assert interface_filename.endswith('_data_interface.yml')
        template_args[interface_filename[:-len('_data_interface.yml')] + '_app_name'] = interface_yml['app'] # TODO: key from submodule name or similar
        template_args[interface_filename[:-len('_data_interface.yml')] + '_input'] = '/'.join(['.'.join(fqn.split('.')[1:]) for fqn in interface_yml[interface_app_name]['input'].values()])
        template_args[interface_filename[:-len('_data_interface.yml')] + '_output'] = '/'.join(['.'.join(fqn.split('.')[1:]) for fqn in interface_yml[interface_app_name]['output'].values()])


    print("Template arguments: {}".format(template_args))

    env = Environment(loader=FileSystemLoader(searchpath=os.path.dirname(args.diagram)))
    template = env.get_template(os.path.basename(args.diagram))
    with open(args.diagram[:-3],'w') as diagram_file:
        diagram_file.write(template.render(**template_args))

