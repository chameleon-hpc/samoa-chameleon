#!/bin/sh

# Load the tracing library (choose C/Fortran)
export LD_PRELOAD=${EXTRAE_HOME}/lib/libompitrace.so    # C

# Run the program
$*

