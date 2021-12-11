
import argparse
from lxml import etree
from lxml.isoschematron import Schematron

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--schema', default='flow_mri.sch',
                    help='Schematron schema file path')
parser.add_argument('--h5dump_headers_xml', default='flow_CStest_Volunteer_R4.h5.xml',
                    help='h5dump -H --xml <MRI>.h5 output file')
args = parser.parse_args()

# https://github.com/openSUSE/schvalidator may be useful for nicer output
sch = Schematron(etree.parse(args.schema_file),
                 store_report=True,
                 error_finder=Schematron.ASSERTS_AND_REPORTS)
sch.validate(etree.parse(args.h5dump_headers_xml))
print(str(sch.validation_report).replace('\xa0',' '))