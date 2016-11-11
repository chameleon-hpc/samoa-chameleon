#!/bin/bash

rm -r output
mkdir output
export OMP_NUM_THREADS=8

bin/samoa_swe_nompi_noasagi_gnu -tmax 24 -dmin 11 -dmax 11 -xmloutput -tout 0.005 -drytolerance 0.00005 -stestpoints 0.5 0.0
