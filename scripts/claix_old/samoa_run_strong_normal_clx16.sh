#!/bin/bash

# normal run using packing/unpacking for fair comparison
unset CHAMELEON_TOOL_LIBRARIES

# use 24 OMP Threads

#export OMP_NUM_THREADS=24

# CCP
#export LB_SETTINGS=1
#export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

#${CUR_MPI_CMD} \
#${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT} \
#${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_packing \
#${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_packing_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log

# no CCP
#export LB_SETTINGS=1000000
#export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

#${CUR_MPI_CMD} \
#${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT} \
#${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_packing \
#${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_packing_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log


# use 23 OMP Threads

export OMP_NUM_THREADS=23

# CCP
#export LB_SETTINGS=1
#export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

#${CUR_MPI_CMD} \
#${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT} \
#${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_packing \
#${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_packing_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log

# no CCP
export LB_SETTINGS=1000000
export SAMOA_PARAMS="${SAMOA_PARAMS} -lbfreq $LB_SETTINGS"

${CUR_MPI_CMD} \
${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT} \
${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_packing \
${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_packing_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log
