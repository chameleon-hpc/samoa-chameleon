#!/bin/bash

source "runs/$1.sh"

cd $(dirname $0)/../../../../../../

echo 'Spawning process'
LC_NUMERIC='en_US.UTF-8' OMP_NUM_THREADS=$THREADS mpiexec -n $RANKS bin/samoa_flash_gnu_${BINARY}_$2 $PARAMS -xdmfoutput -output_dir $3
echo 'Process terminated'