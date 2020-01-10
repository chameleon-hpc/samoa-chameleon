#!/bin/bash

XDMFDIR='/app/samoa_xdmf_libs/intel/parallel'
ASAGIDIR='/app/samoa_xdmf_libs/intel/parallel'
NETCDFDIR='/app/samoa_xdmf_libs/intel/parallel'
SCONSPARAMS="$COMMON_SCONSPARAMS swe_scenario=$1 limiter=$2 compiler=intel xdmf_fox_dir=$XDMFDIR xdmf_hdf5_dir=$XDMFDIR asagi_dir=$ASAGIDIR netcdf_dir=$NETCDFDIR exe=samoa_flash_intel_$1_$2 $BUILDFLAGS"

DOCKER="$(dirname $0)/../../../../../docker"
cd $DOCKER

echo 'Spawning container'
docker-compose run --rm --name "intel_$1_$2_build" samoa_intel /bin/bash -c "source /app/intel/intelrc.sh && scons -j $SCONSTHREADS $SCONSPARAMS"
echo 'Container terminated'