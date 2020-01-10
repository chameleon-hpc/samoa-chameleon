#!/bin/bash

source "runs/$1.sh"

cd $(dirname $0)/../../../../../docker

echo 'Spawning container'
docker-compose run --rm --name "gnu_${BINARY}_$2_run_$1" samoa_gnu /bin/bash -c "LC_NUMERIC='en_US.UTF-8' OMP_NUM_THREADS=$THREADS mpiexec -n $RANKS bin/samoa_flash_gnu_${BINARY}_$2 $PARAMS -xdmfoutput -output_dir $3"
echo 'Container terminated'