#!/bin/bash

HDF5_URL=${HDF5_URL:="$GIT_PROTOCOL://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git"}

FLAGS="--prefix=$3 --enable-fortran --enable-optimization=high --enable-trace --enable-build-mode=production --disable-tests --disable-tools --enable-hl"
if [ "$1" = "mpi" ]; then
    FLAGS="$FLAGS --enable-parallel"
fi

mkdir -p libsrc
cd libsrc
git clone $HDF5_URL
cd hdf5
git checkout hdf5-1_10_4
git reset --hard

echo "$3 HDF5: MPI: $1, Compiler: $2"
mkdir -p "build_$1_$2"
cd "build_$1_$2"
../configure $FLAGS && make -j 8 install