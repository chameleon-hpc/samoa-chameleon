#!/usr/local_rwth/bin/zsh

# normal run using packing/unpacking for fair comparison
unset CHAMELEON_TOOL_LIBRARIES

${CUR_MPI_CMD} \
${PREF_VAR_EXPORT}${ENVS_FOR_EXPORT} \
${SAMOA_VTUNE_PREFIX} no_numa_balancing ${SAMOA_BIN}/samoa_swe${SWE_SCENARIO}${ASAGI_NAME_EXT}_packing \
${SAMOA_PARAMS} &> ${CUR_OUTP_STR}_packing_${OMP_NUM_THREADS}_dmax_${CUR_DMAX}_lbfreq_${LB_SETTINGS}${FILE_APPEND}.log
