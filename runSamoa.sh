#!/bin/bash


export OMP_NUM_THREADS=1
#bin/samoa_swe_nompi_noasagi_fwave_gnu_debug -tsteps 10000 -dmin 8 -dmax 8 -xmloutput -tout 0.0005

bin/samoa_swe_nompi_noasagi_hlle_gnu_debug -tsteps 10000 -dmin 8 -dmax 8 -xmloutput -tout 0.0005

#bin/samoa_swe_nompi_gnu_debug -tsteps 1000000 -dmin 8 -dmax 8 

