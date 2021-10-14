#!/bin/bash

SAMOA_DIR=/home/ps659535/chameleon/samoa_chameleon
SAMOA_PATCH_ORDER=7

module use /home/ps659535/.modules
module purge
module load DEVELOP
module load intel/19.0
module load intelmpi/2018
module load cmake/3.6.0
module load chameleon     
module load python/3.6.0
module load intelitac/2019

CUR_DIR=$(pwd)
cd ${SAMOA_DIR}
# standard builds
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=radial_dam_break flux_solver=aug_riemann assertions=on compiler=intel -j8
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=oscillating_lake flux_solver=aug_riemann assertions=on compiler=intel -j8
scons asagi=on  target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sampl flux_solver=aug_riemann assertions=on compiler=intel -j8

# standard build with packing/unpacking
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=radial_dam_break flux_solver=aug_riemann assertions=on compiler=intel chameleon=2 -j8 
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=oscillating_lake flux_solver=aug_riemann assertions=on compiler=intel chameleon=2  -j8 
scons asagi=on  target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel chameleon=2 -j8

# chameleon builds
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=radial_dam_break flux_solver=aug_riemann assertions=on compiler=intel chameleon=1 -j8
scons asagi=off target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=oscillating_lake flux_solver=aug_riemann assertions=on compiler=intel chameleon=1 -j8 
scons asagi=on target=release scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel chameleon=1 -j8

cd ${CUR_DIR}
