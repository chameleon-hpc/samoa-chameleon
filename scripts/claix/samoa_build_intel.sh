#!/usr/local_rwth/bin/zsh

# parameters that can be set from outside
export CHAMELEON_DIR=${CHAMELEON_DIR:-"/path/to/chameleon/dir"}
export SAMOA_DIR=${SAMOA_DIR:-"/path/to/samoa/dir"}
export ASAGI_DIR=${ASAGI_DIR:-"/path/to/asagi/dir"}
export USE_ASAGI=${USE_ASAGI:-1}
export SWE_SCENARIO=${SWE_SCENARIO:-"oscillating_lake"}
# export SWE_SCENARIO=${SWE_SCENARIO:-"radial_dam_break"}
export SAMOA_PATCH_ORDER=${SAMOA_PATCH_ORDER:-7}
export ENABLE_TRACING_ITAC=${ENABLE_TRACING_ITAC:-0}
export ENABLE_TRACING_EXTRAE=${ENABLE_TRACING_EXTRAE:-0}
export COMPILE_WS=${COMPILE_WS:-0}
export COMPILE_TASKING=${COMPILE_TASKING:-0}
export COMPILE_CHAMELEON=${COMPILE_CHAMELEON:-0}
export COMPILE_PACKING=${COMPILE_PACKING:-0}


echo ">>> Building with the following config:"
echo "CHAMELEON_DIR: ${CHAMELEON_DIR}"
echo "SAMOA_DIR: ${SAMOA_DIR}"
echo "ASAGI_DIR: ${ASAGI_DIR}"
echo "USE_ASAGI: ${USE_ASAGI}"
echo "SWE_SCENARIO: ${SWE_SCENARIO}"
echo "SAMOA_PATCH_ORDER: ${SAMOA_PATCH_ORDER}"
echo "ENABLE_TRACING_ITAC: ${ENABLE_TRACING_ITAC}"
echo "ENABLE_TRACING_EXTRAE: ${ENABLE_TRACING_EXTRAE}"
echo "COMPILE_WS: ${COMPILE_WS}"
echo "COMPILE_TASKING: ${COMPILE_TASKING}"
echo "COMPILE_CHAMELEON: ${COMPILE_CHAMELEON}"
echo "COMPILE_PACKING: ${COMPILE_PACKING}"

# remember current working directory
CUR_DIR=$(pwd)
echo "CUR_DIR=${CUR_DIR}"

if [ "${ENABLE_TRACING_ITAC}" = "1" ]; then
    module load intelitac
fi
if [ "${ENABLE_TRACING_EXTRAE}" = "1" ]; then
    module load DEV-TOOLS extrae papi
    # TODO: Extrae builds might contain a broken Fortran extrae module file.
    #       It is possible to separately build this module and provide the folder path
    #       where the fixed file resides using `export EXTRAE_ADD_INCLUDE=<path>`
fi

# change to samoa dir
cd ${SAMOA_DIR}

# The following versions of samoa can be created here:
# - Regular worksharing based (section indices are calculated for each thread)
# - Tasking version that improves load balance within a process
# - Chameleon version for task-based balancing between processes

# ========== synthetic scenario: oscillating_lake ==========
if [ "${COMPILE_WS}" = "1" ]; then
    scons target=release openmp=notasks \
        asagi=${USE_ASAGI} asagi_dir=${ASAGI_DIR} \
        trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} \
        scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} swe_scenario=${SWE_SCENARIO} \
        data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi
if [ "${COMPILE_TASKING}" = "1" ]; then
    scons target=release openmp=tasks \
        asagi=${USE_ASAGI} asagi_dir=${ASAGI_DIR} \
        trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} \
        scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} swe_scenario=${SWE_SCENARIO} \
        data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi
if [ "${COMPILE_PACKING}" = "1" ]; then
    scons target=release chameleon=2 chameleon_dir=${CHAMELEON_DIR} \
        asagi=${USE_ASAGI} asagi_dir=${ASAGI_DIR} \
        trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} \
        scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} swe_scenario=${SWE_SCENARIO} \
        data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi
if [ "${COMPILE_CHAMELEON}" = "1" ]; then
    scons target=release chameleon=1 chameleon_dir=${CHAMELEON_DIR} \
        asagi=${USE_ASAGI} asagi_dir=${ASAGI_DIR} \
        trace_extrae=${ENABLE_TRACING_EXTRAE} trace_itac=${ENABLE_TRACING_ITAC} \
        scenario=swe swe_patch_order=${SAMOA_PATCH_ORDER} swe_scenario=${SWE_SCENARIO} \
        data_refinement=sample flux_solver=aug_riemann assertions=on compiler=intel mpi=intel -j8
fi

# switch back to working directory
cd ${CUR_DIR}
