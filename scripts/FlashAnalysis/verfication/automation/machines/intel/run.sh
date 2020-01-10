#!/bin/bash

INTELRC=${INTELRC:="~/.intelrc"}

source "runs/$1.sh"
cd $(dirname $0)/../../../../../../

echo 'Spawning process'
source $INTELRC
LC_NUMERIC='en_US.UTF-8' OMP_NUM_THREADS=$THREADS mpiexec -n $RANKS bin/samoa_flash_intel_${BINARY}_$2 $PARAMS -xdmfoutput -output_dir $3
echo 'Process terminated'