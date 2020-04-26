#!/bin/bash

LIBS_R="~/local"
if [ ! -z "$1" ]; then
    LIBS_R=$1
fi

./install_all.sh mpi intel $LIBS_R