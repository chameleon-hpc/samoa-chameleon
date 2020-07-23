#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"


#mkdir output
export OMP_NUM_THREADS=1

periods=1

pi=$(echo "scale=16; 4*a(1)" | bc -l)
#t_max=$(echo "scale=16; $periods*2.0*$pi /sqrt(0.2*9.80665)" | bc -l)
#echo $t_max
t_max=1.0

i_min=10
i_max=10
order=4
output=output/output_${order}_all_$i_max

rm -rf ${output}_1
rm -rf ${output}_2

mkdir ${output}_1
mkdir ${output}_2


bin/samoa_${order}_all_True -nmax 10 -dmin $i_min -dmax $i_max -tout 0.01 -output_dir ${output}_1 -drytolerance 0.000001 -dry_dg_guard 0.01 -coastheight -0.025 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-12 -limiter_buffer 0.1d0 -courant 0.05
bin/samoa_${order}_all_False -nmax 10 -dmin $i_min -dmax $i_max -tout 0.01 -output_dir ${output}_2 -drytolerance 0.000001 -dry_dg_guard 0.01 -coastheight -0.025 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-12 -limiter_buffer 0.1d0 -courant 0.05

#bin/samoa_${order}_all -tmax $t_max -dmin $i_min -dmax $i_max -drytolerance 0.000001 -dry_dg_guard 0.01 -coastheight -0.025 -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-12 -limiter_buffer 0.1d0 -courant 0.05

#gdb --args bin/samoa_${order}_all -tmax $t_max -dmin $i_min -dmax $i_max -tout 0  -output_dir $output -drytolerance 0.00000001 -dry_dg_guard 0.1 -coastheight -0.0001 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-8 -courant 0.15 -limiter_buffer 0.0001d0

#bin/samoa_${order}_all -tmax $t_max -dmin $i_min -dmax $i_max -tout 0 -output_dir $output -drytolerance 0.001 -dry_dg_guard 1010.0 -coastheight 1000.0 -xmlpointoutput -sections 1 -max_picard_iterations 20 -max_picard_error 10.0d-12 -limiter_buffer 100000.0d0 -courant 0.05 -fbath bath.nc -fdispl disp.nc -domain "-2800000 -1950000 1400000 2250000"

