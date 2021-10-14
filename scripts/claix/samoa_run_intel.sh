#!/usr/local_rwth/bin/zsh
#SBATCH --time=01:00:00
#SBATCH --partition=c18m
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --job-name=samoa-osc
#SBATCH --output=samoa-osc.%J.txt

# parameters that can be set from outside
export SAMOA_DIR=${SAMOA_DIR:-"/path/to/samoa/dir"}
export OUTPUT_DIR=${OUTPUT_DIR:-"/path/to/samoa/dir"}
export ENABLE_TRACING_ITAC=${ENABLE_TRACING_ITAC:-0}
export ENABLE_TRACING_EXTRAE=${ENABLE_TRACING_EXTRAE:-0}
export ENABLE_PROFILING=${ENABLE_PROFILING:-0}
export RUN_WS=${RUN_WS:-0}
export RUN_TASKING=${RUN_TASKING:-0}
export SWE_SCENARIO=${SWE_SCENARIO:-"oscillating_lake"}
# export SWE_SCENARIO=${SWE_SCENARIO:-"radial_dam_break"}
export NUM_RANKS=${NUM_RANKS:-2}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-11}

echo ">>> Running with the following config:"
echo "SAMOA_DIR: ${SAMOA_DIR}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "ENABLE_TRACING_ITAC: ${ENABLE_TRACING_ITAC}"
echo "ENABLE_TRACING_EXTRAE: ${ENABLE_TRACING_EXTRAE}"
echo "RUN_WS: ${RUN_WS}"
echo "RUN_TASKING: ${RUN_TASKING}"
echo "SWE_SCENARIO: ${SWE_SCENARIO}"
echo "NUM_RANKS: ${NUM_RANKS}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"

ulimit -c unlimited

if [ "${ENABLE_TRACING_ITAC}" = "1" ]; then
    module load intelitac
fi
if [ "${ENABLE_TRACING_EXTRAE}" = "1" ]; then
    module load DEV-TOOLS extrae papi
    export EXTRAE_CONFIG_FILE=$(pwd)/extrae.xml
    #export CMD_EXTRAE_PREFIX="$(pwd)/trace.sh"
    export LIBRARY_PATH="$LIBRARY_PATH:${EXTRAE_HOME}/lib"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${EXTRAE_HOME}/lib"
fi
if [[ "${ENABLE_PROFILING}" == "1" ]]; then
    module load intelvtune
fi

# create temporary output directory
CUR_DIR=$(pwd)
mkdir -p ${OUTPUT_DIR}

# ========== execution settings ==========
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# choose depending on batch or local environment
if [ -z "${SLURM_JOBID}" ]
then
    export I_MPI_PIN=1 
    export I_MPI_PIN_DOMAIN=auto
    export CUR_MPI_CMD="mpiexec.hydra -n ${NUM_RANKS} -genvlist VT_LOGFILE_PREFIX"
else
    export CUR_MPI_CMD="${MPIEXEC} ${FLAGS_MPI_BATCH} "
fi

# ========== samoa specific settings ==========
export DMIN=15
export DMAX=23
export NUM_SECTIONS=16
export NUM_STEPS=100
export SIM_TIME_SEC=10
export LB_FREQ="-lbfreq 1000000"
export SIM_LIMIT="-nmax ${NUM_STEPS}"
#export SIM_LIMIT="-tmax ${SIM_TIME_SEC}"
export OUTPUT_PARAMS=""
#export OUTPUT_PARAMS="-nout 10 -xmloutput"

# final parameters
export SAMOA_PARAMS="-output_dir ${OUTPUT_DIR} -dmin ${DMIN} -dmax ${DMAX} -sections ${NUM_SECTIONS} ${SIM_LIMIT} ${OUTPUT_PARAMS} ${LB_FREQ}"

# ========== application runs: oscillating lake (worksharing) ==========
if [[ "${RUN_WS}" == "1" ]]; then
    if [[ "${ENABLE_TRACING_ITAC}" == "1" ]]; then
        mkdir -p ${CUR_DIR}/trace_osc_ws
        export VT_LOGFILE_PREFIX=${CUR_DIR}/trace_osc_ws
    fi
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        TMP_EXTRAE="${CUR_DIR}/trace_extrae_osc_ws"
        mkdir -p ${TMP_EXTRAE}
        cd ${TMP_EXTRAE}
    fi
    if [[ "${ENABLE_PROFILING}" == "1" ]]; then
        export CMD_VTUNE_PREFIX="amplxe-cl –collect hotspots –r ${CUR_DIR}/profile_osc_ws_${OMP_NUM_THREADS}t -trace-mpi -- "
    fi
    ${CUR_MPI_CMD} ${CMD_VTUNE_PREFIX} ${CMD_EXTRAE_PREFIX} ${SAMOA_DIR}/bin/samoa_swe_${SWE_SCENARIO}_notasks_intel_noasagi ${SAMOA_PARAMS} &> ${CUR_DIR}/log_samoa_run_osc_ws.log
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        mpi2prv -f TRACE.mpits -v -paraver
    fi
fi

# ========== application runs: oscillating lake (tasking) ==========
if [[ "${RUN_TASKING}" == "1" ]]; then
    if [[ "${ENABLE_TRACING_ITAC}" == "1" ]]; then
        mkdir -p ${CUR_DIR}/trace_osc_tasking
        export VT_LOGFILE_PREFIX=${CUR_DIR}/trace_osc_tasking
    fi
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        TMP_EXTRAE="${CUR_DIR}/trace_extrae_osc_tasking"
        mkdir -p ${TMP_EXTRAE}
        cd ${TMP_EXTRAE}
    fi
    if [[ "${ENABLE_PROFILING}" == "1" ]]; then
        export CMD_VTUNE_PREFIX="amplxe-cl –collect hotspots –r ${CUR_DIR}/profile_osc_tasking_${OMP_NUM_THREADS}t -trace-mpi -- "
    fi
    ${CUR_MPI_CMD} ${CMD_VTUNE_PREFIX} ${CMD_EXTRAE_PREFIX} ${SAMOA_DIR}/bin/samoa_swe_${SWE_SCENARIO}_intel_noasagi ${SAMOA_PARAMS} &> ${CUR_DIR}/log_samoa_run_osc_tasks.log
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        mpi2prv -f TRACE.mpits -v -paraver
    fi
fi
