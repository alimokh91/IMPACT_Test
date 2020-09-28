#!/bin/bash

set -euo pipefail

echo "Removing existing diagram and replacing it by new one"
rm ../doc/diagrams/hpc-predict-app-logic-architecture-combined-template.svg || true

set -x
# $(id -u ${USER}):$(id -g ${USER})
docker run -it --rm -v $(realpath data_interfaces):/src/hpc-predict/hpc-predict-io/data_interfaces -v $(realpath ../doc/diagrams):/src/hpc-predict/doc/diagrams --entrypoint /bin/bash cscs-ci/hpc-predict/io/deploy -c 'source /src/hpc-predict/hpc-predict-io/python/venv/bin/activate && cd /src/hpc-predict/hpc-predict-io/data_interfaces  && PYTHONPATH=/src/hpc-predict/hpc-predict-io/python  python3 validate_pipeline.py --pipeline ../../flownet/data_interface/flownet_data_interface.yml:inference ../../segmenter/data_interface/cnn_segmenter_data_interface.yml:inference ../../IMPACT/data_interface/impact_data_interface.yml:data-assimilation --diagram ../../doc/diagrams/hpc-predict-app-logic-architecture-combined-template.svg.j2'
# inkscape could be run from a container
inkscape ../doc/diagrams/hpc-predict-app-logic-architecture-combined-template.svg --export-type=png --export-background=WHITE
set +x

