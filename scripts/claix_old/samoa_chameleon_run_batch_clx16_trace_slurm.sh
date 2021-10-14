#!/bin/bash
#SBATCH --job-name=samoa_chameleon
#SBATCH --output=output_samoa_chameleon.%J.txt
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --account=jara0001
# Broadwell Nodes
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --partition=c16m

module use /home/ps659535/.modules
module purge
module load DEVELOP
module load intel/19.0
module load intelmpi/2018
module load cmake/3.6.0
module load chameleon/chameleon.test
module load python/3.6.0
module load intelitac
#module load intelvtune/XE2019-u02

#module load intelvtune

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
export CUR_DMAX=25
export NUM_STEPS=300
#export SIM_TIME_SEC=200

export SIM_LIMIT="-nmax ${NUM_STEPS}"
#export SIM_LIMIT="-tmax ${SIM_TIME_SEC}"

#export LB_SETTINGS= 					# w load balancing
export LB_SETTINGS="100000" 	# w/o load balancing

# ===== MPI / OpenMP =====
export OMP_NUM_THREADS=23
#export CUR_MPI_CMD="mpiexec.hydra -n 2 "
export CUR_MPI_CMD="mpiexec -n 4"

# ===== Run / Compare =====

#export SAMOA_VTUNE_PREFIX="amplxe-cl --collect hotspots -r vtune_trace "
export VT_LOGFILE_NAME=samoa_swe_23.stf
source samoa_core_clx16.sh

#source samoa_run_normal.sh
source samoa_run_chameleon.sh
#source samoa_run_tool.sh
