#!/usr/local_rwth/bin/zsh
#BSUB -J "samoa_chameleon"
#BSUB -W 05:00
#BSUB -m c24m128
#BSUB -n 4
#BSUB -x
#BSUB -R "span[ptile=1]"
#BSUB -a intelmpi
#BSUB -P jara0001
#BSUB -M 120000
##BSUB -u samfass@in.tum.de
##BSUB -B
##BSUB -N

#bhist -l ${LSB_JOBID} > ${LSB_JOBNAME}.submission.${LSB_JOBID}

source ~/.bash_profile
module use /home/ps659535/.modules

module load chameleon
module unload openmpi
module load intelmpi/2018.3
module switch intel intel/19.0


# ===== Tohoku Scenario =====
export TOHOKU_PARAMS="-fbath /home/ps659535/chameleon/samoa_chameleon/data/tohoku_static/bath.nc -fdispl /home/ps659535/chameleon/samoa_chameleon/data/tohoku_static/displ.nc"
export ASAGI_NAME_EXT=
export SWE_SCENARIO=

# ===== Other Scenarios =====
#export TOHOKU_PARAMS=
#export ASAGI_NAME_EXT=_noasagi
#export SWE_SCENARIO=_oscillating_lake
#export SWE_SCENARIO=_radial_dam_break

# ===== Simulation Settings =====
export NUM_SECTIONS=16
export CUR_DMIN=15
export CUR_DMAX=23
export NUM_STEPS=150
export SIM_TIME_SEC=3600

#export SIM_LIMIT="-nmax ${NUM_STEPS}"
export SIM_LIMIT="-tmax ${SIM_TIME_SEC}"

#export LB_SETTINGS= 					# w load balancing
export LB_SETTINGS="100" 	# w/o load balancing

# ===== MPI / OpenMP =====
export OMP_NUM_THREADS=24
#export CUR_MPI_CMD="mpiexec.hydra -n 2 "
export CUR_MPI_CMD="$MPIEXEC $FLAGS_MPI_BATCH"

# ===== Run / Compare =====
source samoa_core_clx16.sh

source samoa_run_normal.sh
source samoa_run_chameleon.sh

export CUR_DMAX=25
export OMP_NUM_THREADS=24

source samoa_core_clx16.sh
source samoa_run_normal.sh
source samoa_run_chameleon.sh
