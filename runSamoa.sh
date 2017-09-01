#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

rm -r output
mkdir output
export OMP_NUM_THREADS=8

i_min=10
i_max=10
order=1

mkdir output/output_${order}_unlimited_$i_max

bin/samoa_${order}_unlimited -tmax 4.483 -dmin $i_min -dmax $i_max -tout 5.0  -output_dir output/output_${order}_unlimited_$i_max -drytolerance 0.0001 -coastheight 0.001 -xmlpointoutput -sections 1






