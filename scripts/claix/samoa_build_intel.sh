#!/usr/local_rwth/bin/zsh

# parameters that can be set from outside
export SAMOA_DIR=${SAMOA_DIR:-"/path/to/samoa/dir"}
export ASAGI_DIR=${ASAGI_DIR:-"/path/to/asagi/dir"}
export SAMOA_PATCH_ORDER=${SAMOA_PATCH_ORDER:-7}
export ENABLE_TRACING_ITAC=${ENABLE_TRACING_ITAC:-0}
export ENABLE_TRACING_EXTRAE=${ENABLE_TRACING_EXTRAE:-0}
export COMPILE_WS=${COMPILE_WS:-0}
export COMPILE_TASKING=${COMPILE_TASKING:-0}
export SWE_SCENARIO=${SWE_SCENARIO:-"oscillating_lake"}
# export SWE_SCENARIO=${SWE_SCENARIO:-"radial_dam_break"}

echo ">>> Building with the following config:"
echo "SAMOA_DIR: ${SAMOA_DIR}"
echo "ASAGI_DIR: ${ASAGI_DIR}"
echo "SAMOA_PATCH_ORDER: ${SAMOA_PATCH_ORDER}"
echo "ENABLE_TRACING_ITAC: ${ENABLE_TRACING_ITAC}"
echo "ENABLE_TRACING_EXTRAE: ${ENABLE_TRACING_EXTRAE}"
echo "COMPILE_WS: ${COMPILE_WS}"
echo "COMPILE_TASKING: ${COMPILE_TASKING}"
echo "SWE_SCENARIO: ${SWE_SCENARIO}"

# remember current working directory
CUR_DIR=$(pwd)
echo "CUR_DIR=${CUR_DIR}"

if [ "${ENABLE_TRACING_ITAC}" = "1" ]; then
    module load intelitac
fi
if [ "${ENABLE_TRACING_EXTRAE}" = "1" ]; then
    module load DEV-TOOLS extrae papi
    # TODO: include current directory as it might contain a fixed Fortran extrae module file
fi

# change to samoa dir
cd ${SAMOA_DIR}

# The following versions of samoa can be created here:
# - Regular worksharing based (section indices are calculated for each thread)
# - Tasking version that improves load balance within a process
# - TODO: Chameleon version for task-based balancing between processes

# ========== synthetic scenario: oscillating_lake ==========
if [ "${COMPILE_WS}" = "1" ]; then
    scons asagi=off target=release openmp=notasks trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=${SWE_SCENARIO} flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi
if [ "${COMPILE_TASKING}" = "1" ]; then
    scons asagi=off target=release openmp=tasks   trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} data_refinement=sample swe_scenario=${SWE_SCENARIO} flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi

# ========== ASAGI version ==========
# TODO:

# switch back to working directory
cd ${CUR_DIR}
