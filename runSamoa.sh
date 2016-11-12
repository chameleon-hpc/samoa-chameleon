#!/bin/bash

rm -r output
mkdir output
export OMP_NUM_THREADS=8

for i in `seq 4 20`;
do

mkdir output/output_2_unlimited_$i

bin/samoa_2_unlimited -tmax 0.05 -dmin $i -dmax $i -tout 0.025 -stestpoints 0.0 0.5 -output_dir output/output_2_unlimited_$i -drytolerance 0.00000000001 -xmloutput
done