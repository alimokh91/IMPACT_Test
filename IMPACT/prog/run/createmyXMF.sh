#!/bin/bash

for i in $(seq -f "%08g" 1 10000)
do 
  sed s/00000000/$i/g IMPACT_xdmf3d_00000000.xmf > IMPACT_xdmf3d_$i.xmf
done
