#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

#rm -r output
#mkdir output
export OMP_NUM_THREADS=8

periods=1

pi=$(echo "scale=16; 4*a(1)" | bc -l)
#t_max=$(echo "scale=16; $periods*2.0*$pi /sqrt(0.2*9.80665)" | bc -l)
#echo $t_max
t_max=6.0

i_min=14
i_max=14
order=1
output=output/output_${order}_unlimited_$i_max
mkdir $output


bin/samoa_${order}_unlimited -tmax $t_max -dmin $i_min -dmax $i_max -tout 0.05  -output_dir $output -drytolerance 0.000001 -dry_dg_guard 0.001 -coastheight 0.0001 -xmlpointoutput -sections 1 -courant 0.45






