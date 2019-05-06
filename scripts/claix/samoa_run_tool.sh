#!/bin/bash

# run chameleon with samoa specific tool
#export CHAMELEON_TOOL_LIBRARIES=/home/ps659535/chameleon/chameleon-lib/examples/samoa_tool/samoa_tool.so
export CHAMELEON_TOOL_LIBRARIES=/home/ps659535/chameleon/chameleon-lib/examples/tool_loadbalancing_to_minimum/tool.so

${CUR_MPI_CMD} \
${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT},CHAMELEON_TOOL_LIBRARIES \
${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_chameleon \
${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_tool_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log

