#!/bin/bash

source './machinesetup.sh'

echo 'Starting batch execution'

for RUN in runs/*; do
    RUN=$(basename -- "$RUN")
    RUN=${RUN%.*}
    ./run-single.sh $1 $RUN
done

echo 'Batch execution ended'
