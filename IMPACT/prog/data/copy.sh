#!/bin/bash
echo "copy hdf5 files script"
for ((r=0;r <= 23; r++))#r in {0..23}
    do
    s=`expr $r + 96`
    if [ $r -lt 10 ]; then
        NUM_IN_STR="0$r"
    else
        NUM_IN_STR="$r"
    fi
    if [ $s -lt 10 ]; then
        NUM_IN_ST="0$s"
    else
        NUM_IN_ST="$s"
    fi
    echo "running repetition: $NUM_IN_STR and $NUM_IN_ST"
    cp bern_experimental_dataset_flow_mri_single_rep_${NUM_IN_STR}_dw2.h5 bern_experimental_dataset_flow_mri_single_rep_${NUM_IN_ST}_dw2.h5
    done
