#!/bin/bash

#rm -r output
#mkdir output
export OMP_NUM_THREADS=8

for i in `seq 10 10`;
do

mkdir output/osc_$i

bin/samoa_2_unlimited_osc -tmax 7.0 -dmin $i -dmax $i -tout 0.025 -stestpoints 0.5 0.5 -output_dir output/osc_$i -drytolerance 0.00005 -xmloutput

done