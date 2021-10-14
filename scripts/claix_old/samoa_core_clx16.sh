#!/bin/bash

export SAMOA_BIN=/home/ps659535/chameleon/samoa_chameleon/bin
#export SAMOA_OUTPUT_DIR=/work/jk869269/repos/chameleon/samoa_xmloutput
export SAMOA_OUTPUT_DIR=/home/ps659535/chameleon/samoa_chameleon/output

export CUR_DATE_STR="$(date +"%Y%m%d_%H%M%S")"
export CUR_OUTP_STR="${CUR_DATE_STR}_output"
export CUR_TIMES_STR="${CUR_DATE_STR}_times"

#ulimit -s unlimited
#ulimit -c unlimited

export I_MPI_PIN=1
#export I_MPI_DEBUG=5
export I_MPI_PIN_DOMAIN=auto
export I_MPI_TMI_NBITS_RANK=16 #set env var for enabling larger tag size on OPA with IntelMPI and psm2
export I_MPI_FABRICS="shm:tmi" #select different fabric because of weird IntelMPI bug that truncates messages or switches order... don't know

export OMP_PLACES=cores
#export OMP_PROC_BIND=spread
export OMP_PROC_BIND=close
#export KMP_AFFINITY=verbose

export SAMOA_PARAMS="-output_dir ${SAMOA_OUTPUT_DIR} -dmin ${CUR_DMIN} -dmax ${CUR_DMAX} -sections ${NUM_SECTIONS} ${TOHOKU_PARAMS} ${SIM_LIMIT} "
export PREF_VAR_EXPORT="-genvlist "
export ENVS_FOR_EXPORT="PATH,CPLUS_INCLUDE_PATH,C_INCLUDE_PATH,CPATH,INCLUDE,LD_LIBRARY_PATH,LIBRARY_PATH,I_MPI_PIN,I_MPI_DEBUG,I_MPI_PIN_DOMAIN,I_MPI_TMI_NBITS_RANK,OMP_NUM_THREADS,OMP_PLACES,OMP_PROC_BIND,I_MPI_FABRICS,MIN_ABS_LOAD_IMBALANCE_BEFORE_MIGRATION,CHAMELEON_TOOL_LIBRARIES"

# reduce number of threads by one due to additional communication thread
#export OMP_NUM_THREADS=$((${OMP_NUM_THREADS}-1))



