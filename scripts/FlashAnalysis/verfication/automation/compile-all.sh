#!/bin/bash

source './machinesetup.sh'

echo 'Starting batch compilation'

for SCENARIO in builds/*; do
    SCENARIO=$(basename -- "$SCENARIO")
    SCENARIO=${SCENARIO%.*}
    ./compile-single.sh $1 $SCENARIO
done

echo 'Batch compilation ended'
