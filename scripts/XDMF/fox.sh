#!/bin/bash

FLAGS="-DCMAKE_INSTALL_PREFIX=$3 -DBUILD_SHARED_LIBS=1 -DFoX_ENABLE_DOM=1 -DFoX_ENABLE_EXAMPLES=0 -DFoX_ENABLE_WCML=0 -DFoX_ENABLE_WKML=0 -DFoX_ENABLE_WXML=0 -DFoX_SUPPRESS_WARNINGS=1"

mkdir -p libsrc
cd libsrc
git clone "$GIT_PROTOCOL://github.com/StarGate01/fox.git"
cd fox
git checkout 3199302381316608fb11daad9bf93ae2166bfe13
git reset --hard

echo "$3 FoX: MPI: $1, Compiler: $2"
mkdir -p "build_$1_$2"
cd "build_$1_$2"
cmake .. $FLAGS && make -j 8 install