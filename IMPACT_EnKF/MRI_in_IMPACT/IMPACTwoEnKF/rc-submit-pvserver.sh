#!/bin/bash -l

# tested for ParaView 5.8

# usage
#echo "Usage : %1:Session name ($1)"
#echo "        %2:Job Wall Time ($2)"
#echo "        %3:server-num-nodes ($3)"
#echo "        %4:server-num-tasks-per-node ($4)"
#echo "        %5:server-port ($5)"
#echo "        %6:login node ($6)"
#echo "        %7:Version number ("GNU-5.8" or "GNU-5.8-OSMesa") ($7)"
#echo "        %8:Queue's name ("normal" or "debug") ($8)"
#echo "        %9:Memory per Node ("standard" or "high")($9)"

# Create a temporary filename to write our launch script into
TEMP_FILE=`mktemp`
echo "Temporary FileName is :" $TEMP_FILE

nservers=$[$3 * $4]

# Create a job script
echo "#!/bin/bash"                                    >> $TEMP_FILE
echo "#SBATCH --job-name=$1"                          >> $TEMP_FILE
echo "#SBATCH --nodes=$3"                             >> $TEMP_FILE
echo "#SBATCH --ntasks-per-node=$4"                   >> $TEMP_FILE
echo "#SBATCH --ntasks=$nservers"                     >> $TEMP_FILE
echo "#SBATCH --time=$2"                              >> $TEMP_FILE
echo "#SBATCH --partition=$8"                         >> $TEMP_FILE
echo "#SBATCH --account=c12"       >> $TEMP_FILE

if [ "$7" = "GNU-5.9" ]; then
echo "#SBATCH --constraint=gpu"                    >> $TEMP_FILE
echo "module load daint-gpu"                       >> $TEMP_FILE
echo "module load /scratch/snx3000/jfavre/daint/modules/all/oidn/1.2.4-CrayGNU-20.08" >> $TEMP_FILE 
echo "module load /scratch/snx3000/jfavre/daint/modules/all/ospray/2.4.0-CrayGNU-20.08" >> $TEMP_FILE
echo "module load /scratch/snx3000/jfavre/daint/modules/all/ParaView/5.9.0-RC2-CrayGNU-20.08-EGL-python3" >> $TEMP_FILE
echo "export NVINDEX_PVPLUGIN_HOME=/apps/common/UES/easybuild/sources/p/ParaView" >> $TEMP_FILE

elif [ "$7" = "GNU-5.9-OSMesa" ]; then
if [ "$9" = "high" ]; then
echo "#SBATCH --mem=122G"                    >> $TEMP_FILE
fi
echo "#SBATCH --cpus-per-task=72"                  >> $TEMP_FILE
echo "#SBATCH --ntasks-per-core=2"                 >> $TEMP_FILE
echo "#SBATCH --hint=multithread"                       >> $TEMP_FILE
echo "#SBATCH --constraint=mc"                    >> $TEMP_FILE
echo "module load daint-mc"                       >> $TEMP_FILE
echo "export GALLIUM_DRIVER=swr"                       >> $TEMP_FILE
echo "module load /scratch/snx3000/jfavre/daint/modules/all/oidn/1.2.4-CrayGNU-20.08" >> $TEMP_FILE
echo "module load /scratch/snx3000/jfavre/daint/modules/all/ospray/2.4.0-CrayGNU-20.08" >> $TEMP_FILE
echo "module load /scratch/snx3000/jfavre/daint/modules/all/ParaView/5.9.0-RC2-CrayGNU-20.08-OSMesa-python3" >> $TEMP_FILE

elif [ "$7" = "GNU-5.8" ]; then
echo "#SBATCH --constraint=gpu"                       >> $TEMP_FILE
echo "module load daint-gpu"                          >> $TEMP_FILE
echo "module load ParaView/5.8.1-CrayGNU-20.08-EGL-python3"  >> $TEMP_FILE

elif [ "$7" = "GNU-5.8-OSMesa" ]; then
echo "#SBATCH --constraint=mc"                    >> $TEMP_FILE
echo "module load daint-mc"                                      >> $TEMP_FILE
echo "module load ParaView/5.8.1-CrayGNU-20.08-OSMesa-python3"  >> $TEMP_FILE

if [ "$9" = "high" ]; then
echo "#SBATCH --mem=122G"                    >> $TEMP_FILE
fi
fi

echo ""                                               >> $TEMP_FILE
echo ""                                               >> $TEMP_FILE
echo "srun -n $nservers -N $3 --cpu_bind=sockets pvserver -rc -ch=$6 -sp=$5" >> $TEMP_FILE
cat $TEMP_FILE

# submit the job
sbatch $TEMP_FILE

# wipe the temp file
rm $TEMP_FILE
