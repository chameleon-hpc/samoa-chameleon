#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"


#mkdir output
export OMP_NUM_THREADS=1

periods=1

pi=$(echo "scale=16; 4*a(1)" | bc -l)
#t_max=$(echo "scale=16; $periods*2.0*$pi /sqrt(0.2*9.80665)" | bc -l)
#echo $t_max
t_max=1000.0

i_min=10
i_max=12
order=2
output=output/output_${order}_all_$i_max

rm -r $output
mkdir $output



#bin/samoa_${order}_${limiter} -tmax $t_max -dmin $i_min -dmax $i_max -tout 0.01 -output_dir $output -drytolerance 0.001 -dry_dg_guard 0.03 -coastheight 0.025 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-12 -limiter_buffer 0.0001d0 -courant 0.05

mpiexec -np 4 bin/samoa_${order}_all -tmax $t_max -dmin $i_min -dmax $i_max -tout 0.1  -output_dir $output -drytolerance 0.000001 -dry_dg_guard -0.001 -coastheight -0.0001 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-8 -courant 0.05 -limiter_buffer 0.00001d0

