#!/bin/bash

XDMFDIR='/app/samoa_xdmf_libs/gnu/parallel'
ASAGIDIR='/app/samoa_xdmf_libs/gnu/parallel'
NETCDFDIR='/app/samoa_xdmf_libs/gnu/parallel'
SCONSPARAMS="$COMMON_SCONSPARAMS swe_scenario=$1 limiter=$2 compiler=gnu xdmf_fox_dir=$XDMFDIR xdmf_hdf5_dir=$XDMFDIR asagi_dir=$ASAGIDIR netcdf_dir=$NETCDFDIR exe=samoa_flash_gnu_$1_$2 $BUILDFLAGS"

cd $(dirname $0)/../../../../../docker

echo 'Spawning container'
docker-compose run --rm --name "gnu_$1_$2_build" samoa_gnu /bin/bash -c "scons -j $SCONSTHREADS $SCONSPARAMS"
echo 'Container terminated'