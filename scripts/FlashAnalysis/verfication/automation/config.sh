#!/bin/bash

export SCONSTHREADS=${SCONSTHREADS:=8}
export THREADS=${THREADS:=4}
export RANKS=${RANKS:=1}

export COMMON_SCONSPARAMS='scenario=flash flash_order=1 mpi=default target=release data_refinement=sample xdmf=true'

export LIMITERS='BJ_edge BJ_vertex'

declare -A INDICES
INDICES[resting_lake]=1
INDICES[linear_beach]=2
INDICES[longwave_basin]=3
INDICES[oscillating_lake]=4
INDICES[okushiri]=5
INDICES[conical_island]=6

echo "Configuration: Scons build threads = $SCONSTHREADS, OpenMP threads = $THREADS, MPI ranks = $RANKS"
