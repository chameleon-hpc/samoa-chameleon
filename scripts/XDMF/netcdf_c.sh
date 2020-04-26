#!/bin/bash

FLAGS="--prefix=$3 --disable-doxygen --disable-filter-testing --disable-dap --disable-examples --disable-testsets"
if [ "$1" = "mpi" ]; then
    FLAGS="$FLAGS --enable-parallel4"
else
    FLAGS="$FLAGS --disable-parallel4"
fi

mkdir -p libsrc
cd libsrc
git clone "$GIT_PROTOCOL://github.com/Unidata/netcdf-c.git"
cd netcdf-c
git checkout v4.7.1
git reset --hard

PATH=$PATH:$3/bin

echo "$3 NetCDF C: MPI: $1, Compiler: $2"
mkdir -p "build_$1_$2"
cd "build_$1_$2"
CPPFLAGS=-I$3/include LDFLAGS=-L$3/lib ../configure $FLAGS && make -j 8 install