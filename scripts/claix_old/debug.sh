#!/usr/bin/env zsh
#BSUB -J terminal
#BSUB -W 03:00
#BSUB -m c24m128
#BSUB -n 4
#BSUB -R "span[ptile=1]"
#BSUB -x
#BSUB -a intelmpi
#BSUB -P jara0001
#BSUB -M 120000
#BSUB -XF

module use ~/.modules
module load chameleon
module unload intel
module load intel/19.0
module switch openmpi intelmpi/2018.3

OMP_PLACES=cores 
OMP_PROC_BIND=spread 
KMP_AFFINITY=verbose
I_MPI_PIN=1 
I_MPI_PIN_DOMAIN=auto 
I_MPI_DEBUG=5

xterm

