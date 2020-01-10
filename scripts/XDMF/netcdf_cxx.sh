#!/bin/bash

FLAGS="--prefix=$3 --disable-doxygen"

mkdir -p libsrc
cd libsrc
git clone "$GIT_PROTOCOL://github.com/Unidata/netcdf-cxx4.git"
cd netcdf-cxx4
git checkout v4.3.1
git reset --hard

PATH=$PATH:$3/bin

echo "$3 NetCDF C++: MPI: $1, Compiler: $2"
mkdir -p "build_$1_$2"
cd "build_$1_$2"
CPPFLAGS=-I$3/include LDFLAGS=-L$3/lib ../configure $FLAGS && make -j 8 install