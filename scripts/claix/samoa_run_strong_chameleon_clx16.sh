#!/bin/bash

# normal run using packing/unpacking for fair comparison
unset CHAMELEON_TOOL_LIBRARIES


# use 23 OMP Threads

export OMP_NUM_THREADS=23

FILE_APPEND_TMP=${FILE_APPEND}

# CCP
#export LB_SETTINGS=1
#export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

#export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=0
#export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_0"
#source samoa_run_chameleon.sh

#export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=2
#export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_2"
#source samoa_run_chameleon.sh

#export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=23
#export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_23"
#source samoa_run_chameleon.sh


# no CCP
export LB_SETTINGS=1000000
export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

#export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=0
#export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_0"
#source samoa_run_chameleon.sh

#export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=2
#export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_2"
#source samoa_run_chameleon.sh

export MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION=23
export FILE_APPEND="${FILE_APPEND_TMP}_abs_before_migration_23"
source samoa_run_chameleon.sh

