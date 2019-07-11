#!/bin/bash

set -euxo pipefail

IMPACT_DEBUG_INFO_DIR=/tmp/impact_debug_info
cat ${IMPACT_DEBUG_INFO_DIR}/xterm_processes | sort | tee ${IMPACT_DEBUG_INFO_DIR}/xterm_processes
