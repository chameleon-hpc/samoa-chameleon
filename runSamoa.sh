#!/bin/bash

rm -r output
mkdir output
export OMP_NUM_THREADS=8
#bin/samoa_swe_nompi_noasagi_llf_gnu_debug -tsteps 10000 -dmin 10 -dmax 10 -xmloutput -tout 0.005
#bin/samoa_swe_nompi_noasagi_fwave_gnu_debug -tsteps 10000 -dmin 10 -dmax 10 -xmloutput -tout 0.0005
#bin/samoa_swe_nompi_noasagi_hlle_gnu_debug -tsteps 10000 -dmin 12 -dmax 12 -xmloutput -tout 0.005 -drytolerance 1.0e-10
#bin/samoa_swe_nompi_noasagi_hlle_gnu_debug -tsteps 10000 -dmin 10 -dmax 10 -xmloutput -tout 0.005 -drytolerance 0.01
bin/samoa_swe_nompi_noasagi_hlle_gnu_debug -tsteps 2 -dmin 10 -dmax 10 -xmloutput -tout 0.05 -drytolerance 0.05
