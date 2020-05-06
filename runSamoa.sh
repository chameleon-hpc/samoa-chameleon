#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"


#mkdir output
export OMP_NUM_THREADS=4

periods=1

pi=$(echo "scale=16; 4*a(1)" | bc -l)
#t_max=$(echo "scale=16; $periods*2.0*$pi /sqrt(0.2*9.80665)" | bc -l)
#echo $t_max
t_max=6.0

i_min=8
i_max=8
order=4
output=output/output_${order}_unlimited_$i_max
rm -r $output
mkdir $output


bin/samoa_${order}_unlimited -tmax $t_max -dmin $i_min -dmax $i_max -tout 0.1  -output_dir $output -drytolerance 0.000001 -dry_dg_guard -0.001 -coastheight -0.0001 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-8 -courant 0.05
