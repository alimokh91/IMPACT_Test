IMPACT_DEBUG_INFO_DIR=/tmp/impact_debug_info
flock -x -w 5 ${IMPACT_DEBUG_INFO_DIR}/xterm_processes echo $$ >> ${IMPACT_DEBUG_INFO_DIR}/xterm_processes

