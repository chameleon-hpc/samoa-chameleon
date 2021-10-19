#!/usr/local_rwth/bin/zsh
#SBATCH --time=01:00:00
#SBATCH --partition=c18m
#SBATCH --exclusive
#SBATCH --job-name=samoa-job
#SBATCH --output=samoa-job.%J.out
#SBATCH --error=samoa-job.%J.err

# parameters that can be set from outside
export OUTPUT_DIR=${OUTPUT_DIR:-"/path/to/samoa/dir"}
export SAMOA_DIR=${SAMOA_DIR:-"/path/to/samoa/dir"}
export ASAGI_DATA_DIR=${ASAGI_DATA_DIR:-"/path/to/asagi/data/dir"}
export USE_ASAGI=${USE_ASAGI:-1}
export ENABLE_TRACING_ITAC=${ENABLE_TRACING_ITAC:-0}
export ENABLE_TRACING_EXTRAE=${ENABLE_TRACING_EXTRAE:-0}
export ENABLE_PROFILING=${ENABLE_PROFILING:-0}
export RUN_WS=${RUN_WS:-0}
export RUN_TASKING=${RUN_TASKING:-0}
export RUN_CHAMELEON=${RUN_CHAMELEON:-0}
export RUN_PACKING=${RUN_PACKING:-0}
export SWE_SCENARIO=${SWE_SCENARIO:-"oscillating_lake"}
# export SWE_SCENARIO=${SWE_SCENARIO:-"radial_dam_break"}
export NUM_RANKS=${NUM_RANKS:-2}
export OMP_NUM_THREADS_VAR=${OMP_NUM_THREADS_VAR:-11}

echo ">>> Running with the following config:"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "SAMOA_DIR: ${SAMOA_DIR}"
echo "ASAGI_DATA_DIR: ${ASAGI_DATA_DIR}"
echo "USE_ASAGI: ${USE_ASAGI}"
echo "ENABLE_TRACING_ITAC: ${ENABLE_TRACING_ITAC}"
echo "ENABLE_TRACING_EXTRAE: ${ENABLE_TRACING_EXTRAE}"
echo "ENABLE_PROFILING: ${ENABLE_PROFILING}"
echo "RUN_WS: ${RUN_WS}"
echo "RUN_TASKING: ${RUN_TASKING}"
echo "RUN_CHAMELEON: ${RUN_CHAMELEON}"
echo "RUN_PACKING: ${RUN_PACKING}"
echo "SWE_SCENARIO: ${SWE_SCENARIO}"
echo "NUM_RANKS: ${NUM_RANKS}"
echo "OMP_NUM_THREADS_VAR: ${OMP_NUM_THREADS_VAR}"

ulimit -c unlimited

if [ "${ENABLE_TRACING_ITAC}" = "1" ]; then
    module load intelitac
fi
if [ "${ENABLE_TRACING_EXTRAE}" = "1" ]; then
    module load DEV-TOOLS extrae papi
    export EXTRAE_CONFIG_FILE=$(pwd)/extrae.xml
    export LIBRARY_PATH="$LIBRARY_PATH:${EXTRAE_HOME}/lib"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${EXTRAE_HOME}/lib"
fi
if [[ "${ENABLE_PROFILING}" == "1" ]]; then
    module load intelvtune
fi
if [[ "${RUN_CHAMELEON}" == "1" ]]; then
    # Append linux env vars with Chameleon include and lib folder (here: realized using an environment module)
    source ~/.zshrc
    module load chameleon
fi
if [[ "${USE_ASAGI}" == "0" ]]; then
    export ASAGI_NAME_EXT="_noasagi"
    export SWE_NAME_EXT="${SWE_SCENARIO}_"
fi

# create temporary output directory
CUR_DIR=$(pwd)
mkdir -p ${OUTPUT_DIR}

# ========== execution settings ==========
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=${OMP_NUM_THREADS_VAR}

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
export SIM_TIME_SEC=600
export LB_FREQ="-lbfreq 1000000"
export ASAGI_PARAMS="-fbath ${ASAGI_DATA_DIR}/bath.nc -fdispl ${ASAGI_DATA_DIR}/displ.nc"
# export SIM_LIMIT="-nmax ${NUM_STEPS}"
export SIM_LIMIT="-tmax ${SIM_TIME_SEC}"
export OUTPUT_PARAMS=""
#export OUTPUT_PARAMS="-nout 10 -xmloutput"

# final parameters
export SAMOA_PARAMS="-output_dir ${OUTPUT_DIR} -dmin ${DMIN} -dmax ${DMAX} -sections ${NUM_SECTIONS} ${ASAGI_PARAMS} ${SIM_LIMIT} ${OUTPUT_PARAMS} ${LB_FREQ}"

run_experiment() {
    if [[ "${ENABLE_TRACING_ITAC}" == "1" ]]; then
        mkdir -p ${CUR_DIR}/trace_itac_${NAME_EXTENSION}
        export VT_LOGFILE_PREFIX=${CUR_DIR}/trace_itac_${NAME_EXTENSION}
    fi
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        TMP_EXTRAE="${CUR_DIR}/trace_extrae_${NAME_EXTENSION}"
        mkdir -p ${TMP_EXTRAE}
        cd ${TMP_EXTRAE}
    fi
    if [[ "${ENABLE_PROFILING}" == "1" ]]; then
        export CMD_VTUNE_PREFIX="amplxe-cl –collect hotspots –r ${CUR_DIR}/profile_${NAME_EXTENSION}_${OMP_NUM_THREADS}t -trace-mpi -- "
    fi
    ${CUR_MPI_CMD} ${CMD_VTUNE_PREFIX} ${CMD_EXTRAE_PREFIX} ${SAMOA_DIR}/bin/${EXE_NAME} ${SAMOA_PARAMS} &> ${CUR_DIR}/log_samoa_run_${NAME_EXTENSION}.log
    if [[ "${ENABLE_TRACING_EXTRAE}" == "1" ]]; then
        mpi2prv -f TRACE.mpits -v -paraver
    fi
}

# ========== application runs: oscillating lake (work-sharing) ==========
if [[ "${RUN_WS}" == "1" ]]; then
    echo "Running work-sharing version"
    export EXE_NAME=samoa_swe_${SWE_NAME_EXT}notasks_intel${ASAGI_NAME_EXT}
    export NAME_EXTENSION=ws
    run_experiment
fi

# ========== application runs: oscillating lake (tasking) ==========
if [[ "${RUN_TASKING}" == "1" ]]; then
    echo "Running tasking version"
    export EXE_NAME=samoa_swe_${SWE_NAME_EXT}intel${ASAGI_NAME_EXT}
    export NAME_EXTENSION=tasking
    run_experiment
fi

# ========== application runs: oscillating lake (packing = OpenMP tasks) ==========
if [[ "${RUN_PACKING}" == "1" ]]; then
    echo "Running packing version"
    export EXE_NAME=samoa_swe_${SWE_NAME_EXT}intel${ASAGI_NAME_EXT}_packing
    export NAME_EXTENSION=packing
    run_experiment
fi

# ========== application runs: oscillating lake (chameleon) ==========
if [[ "${RUN_CHAMELEON}" == "1" ]]; then
    echo "Running Chameleon version"
    export EXE_NAME=samoa_swe_${SWE_NAME_EXT}intel${ASAGI_NAME_EXT}_chameleon
    export NAME_EXTENSION=chameleon
    run_experiment
fi
