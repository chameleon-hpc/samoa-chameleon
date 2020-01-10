#!/bin/bash

XDMFDIR=${XDMFDIR:='/opt/samoa_xdmf_libs/gnu/parallel'}
ASAGIDIR=${ASAGIDIR:='/opt/samoa_xdmf_libs/gnu/parallel'}
NETCDFDIR=${NETCDFDIR:='/opt/samoa_xdmf_libs/gnu/parallel'}
SCONSPARAMS="$COMMON_SCONSPARAMS swe_scenario=$1 limiter=$2 compiler=gnu xdmf_fox_dir=$XDMFDIR xdmf_hdf5_dir=$XDMFDIR asagi_dir=$ASAGIDIR netcdf_dir=$NETCDFDIR exe=samoa_flash_gnu_$1_$2 $BUILDFLAGS"

cd $(dirname $0)/../../../../../../

echo 'Spawning process'
scons -j $SCONSTHREADS $SCONSPARAMS
echo 'Process terminated'