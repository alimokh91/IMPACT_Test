#!/bin/bash

h5dump --xml -H "$1" | sed 's/\/\/hdfgroup.org/\/\/support.hdfgroup.org/g' | sed 's/www.hdfgroup.org/support.hdfgroup.org/g' > "$1.xml"

# For now just assume it's a flow-mri (could use the group name in choosing the Schematron schema file)
python h5_validate_schematron.py --schema flow_mri.sch --h5dump_headers_xml "$1.xml"
