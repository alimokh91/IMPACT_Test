IMPACT_DEBUG_INFO_DIR=/tmp/impact_debug_info

echo $$ >> ${IMPACT_DEBUG_INFO_DIR}/xterm_processes/$$

sleep 0.5

unset xterm_processes
unset mpi_processes
#readarray xterm_processes < ${IMPACT_DEBUG_INFO_DIR}/xterm_processes

xterm_processes=()
while IFS=' ' read -r value; do
    xterm_processes+=( $(echo ${value} | tr -d '\n') )
done < <(cat ${IMPACT_DEBUG_INFO_DIR}/xterm_processes/*)
#done < ${IMPACT_DEBUG_INFO_DIR}/xterm_processes

declare -A mpi_processes

while IFS=' ' read -r key value; do
    mpi_processes[$key]=$(echo ${value} | tr -d '\n')
done < <(cat ${IMPACT_DEBUG_INFO_DIR}/mpi_processes/*)
#done < ${IMPACT_DEBUG_INFO_DIR}/mpi_processes

for i in "${!mpi_processes[@]}"; do printf "\"%s\":\"%s\"\n" "$i" "${mpi_processes[$i]}";done

mpi_process=-1

for i in "${!xterm_processes[@]}"; do
   echo "Comparing \"${xterm_processes[$i]}\" to \"$$\"..."
   if [[ "${xterm_processes[$i]}" = "$$" ]]; then
       mpi_process=${mpi_processes["${i}"]};
   fi
done

if [[ mpi_process = -1 ]]; then
    echo "No corresponding MPI process to attach to found for this xterm... exiting."
    exit -1
fi

gdb attach ${mpi_process}
