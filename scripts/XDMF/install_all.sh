#!/bin/bash

function print_help()
{
    echo "Specify MPI support (mpi | nompi) as first parameter"
    echo "Specify compiler choice (gnu | intel) as second parameter"
    echo "Specify installation directory as third parameter"
}

if [ "$1" = "mpi" ] || [ "$1" = "nompi" ]; then
    MPI=$1
else
    echo "Unknown MPI option $1"
    print_help
    exit 1
fi
if [ "$2" = "gnu" ] || [ "$2" = "intel" ]; then
    COMP=$2
else
    echo "Unknown compiler option $2"
    print_help
    exit 1
fi
if [ ! -z "$3" ]; then
    LIBS_R=$3
else
    echo "Missing install directory"
    print_help
    exit 1
fi

./install_lib.sh fox $MPI $COMP $LIBS_R
./install_lib.sh hdf5 $MPI $COMP $LIBS_R
./install_lib.sh netcdf_c $MPI $COMP $LIBS_R
./install_lib.sh netcdf_cxx $MPI $COMP $LIBS_R
./install_lib.sh asagi $MPI $COMP $LIBS_R

