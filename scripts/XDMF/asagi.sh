#!/bin/bash

#todo korrekte netcdf libs

ASAGI_NUMA=${ASAGI_NUMA:='-DNONUMA=1'}

FLAGS="-DCMAKE_INSTALL_PREFIX=$3 -DNETCDF_LIBRARIES=$3/lib -DNETCDF_INCLUDES=$3/include $ASAGI_NUMA"
if [ "$1" = "nompi" ]; then
    FLAGS="$FLAGS -DNOMPI=1"
fi

mkdir -p libsrc
cd libsrc
git clone "$GIT_PROTOCOL://github.com/StarGate01/ASAGI.git"
cd ASAGI
git checkout 3e572ee2c6f8f816dd551418e5c6d9c8900bf7a7
git reset --hard
sed -i "s/https/$GIT_PROTOCOL/g" .gitmodules
git submodule update --init --recursive

PATH=$PATH:$3/bin

echo "$3 ASAGI: MPI: $1, Compiler: $2"
mkdir -p "build_$1_$2"
cd "build_$1_$2"
CPPFLAGS=-I$3/include LDFLAGS=-L$3/lib cmake .. $FLAGS && make -j 8 install