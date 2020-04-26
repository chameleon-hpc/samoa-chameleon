#!/bin/bash

function print_help()
{
    echo "Specify library name (hdf5 | fox) as first parameter"
    echo "Specify MPI support (mpi | nompi) as second parameter"
    echo "Specify compiler choice (gnu | intel) as third parameter"
    echo "Specify installation directory as fourth parameter"
}

MPI="nompi"
if [ ! -z "$2" ]; then
    MPI=$2
fi
COMP="gnu"
if [ ! -z "$3" ]; then
    COMP=$3
fi
if [ ! -z "$4" ]; then
    LIBS_R="$4/$COMP"
else
    echo "Missing install directory"
    print_help
    exit 1
fi

export GIT_PROTOCOL=${GIT_PROTOCOL:='https'}

if [ "$MPI" = "mpi" ]; then
    LIBS_R="$LIBS_R/parallel"
    if [ "$COMP" = "gnu" ]; then
        export CC=mpicc
        export MPICC=mpicc
        export CXX=mpic++
        export MPICXX=mpic++
        export FC=mpifort
        export F90=mpifort
        export F9X=mpifort
        export MPIFC=mpifort
        export MPIF90=mpifort
        export MPIF9X=mpifort
        export MPICH_F90=gfortran
        export OMPI_F90=gfortran
        export I_MPI_F90=gfortran
        export MPICH_F77=gfortran
        export OMPI_F77=gfortran
        export I_MPI_F77=gfortran
        export MPICH_FC=gfortran
        export OMPI_FC=gfortran
        export I_MPI_FC=gfortran
        export MPICH_F9X=gfortran
        export OMPI_F9X=gfortran
        export I_MPI_F9X=gfortran
        export MPICH_CC=gcc
        export OMPI_CC=gcc
        export I_MPI_CC=gcc
        export MPICH_CXX=g++
        export OMPI_CXX=g++
        export I_MPI_CXX=g++
    elif [ "$COMP" = "intel" ]; then
        export CC=mpiicc
        export MPICC=mpiicc
        export CXX=mpiicpc
        export MPICXX=mpiicpc
        export FC=mpiifort
        export F90=mpiifort
        export F9X=mpiifort
        export MPIFC=mpiifort
        export MPIF90=mpiifort
        export MPIF9X=mpiifort
        export MPICH_F90=ifort
        export OMPI_F90=ifort
        export I_MPI_F90=ifort
        export MPICH_F77=ifort
        export OMPI_F77=ifort
        export I_MPI_F77=ifort
        export MPICH_FC=ifort
        export OMPI_FC=ifort
        export I_MPI_FC=ifort
        export MPICH_F9X=ifort
        export OMPI_F9X=ifort
        export I_MPI_F9X=ifort
        export MPICH_CC=icc
        export OMPI_CC=icc
        export I_MPI_CC=icc
        export MPICH_CXX=icpc
        export OMPI_CXX=icpc
        export I_MPI_CXX=icpc
    else
        echo "Unknown compiler choice $COMP"
        print_help
        exit 1
    fi
elif [ "$MPI" = "nompi" ]; then
    LIBS_R="$LIBS_R/serial"
    if [ "$COMP" = "gnu" ]; then
        export CC=gcc
        export CXX=g++
        export FC=gfortran
    elif [ "$COMP" = "intel" ]; then
        export CC=icc
        export CXX=icpc
        export FC=ifort
    else
        echo "Unknown compiler choice $COMP"
        print_help
        exit 1
    fi
else
    echo "Unknown MPI option $MPI"
    print_help
    exit 1
fi

if [ -f "$1.sh" ]; then
    if ! mkdir -p $LIBS_R; then
        echo "Missing permission to write to installation directory $LIBS_R"
        print_help
        exit 1
    fi
    LIBS=`cd "$LIBS_R"; pwd`
    echo "Compiling and installing $1 (MPI: $MPI) using $COMP compiler to $LIBS"
    ./$1.sh $MPI $COMP $LIBS
else
    echo "Unknown library option $1"
    print_help
    exit 1
fi
