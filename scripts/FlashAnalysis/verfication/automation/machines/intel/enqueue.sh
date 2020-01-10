#!/bin/bash

TIME=${TIME:='00:30:00'}

source "runs/$1.sh"

echo 'Generating slurm config'
NAME="${BINARY}_$2_run_$1"
SCONFIG="$CURRENT/slurm_$NAME.sh"
sed "s/JOBNAME/$NAME/g; s/LIMITER/$2/g; s/RANKS/$RANKS/g; s/TIME/$TIME/g; s/SCENARIO/$BINARY/g; s/PARAMS/${PARAMS//\//\\/}/g; s/NTHREADS/$THREADS/g;" "$PROFILE/slurm.template" > $SCONFIG

echo 'Enqueueing slurm job'
sbatch $SCONFIG