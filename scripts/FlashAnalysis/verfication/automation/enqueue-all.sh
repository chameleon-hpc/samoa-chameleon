#!/bin/bash

source './machinesetup.sh'

echo 'Starting batch enqueueing'

for RUN in runs/*; do
    RUN=$(basename -- "$RUN")
    RUN=${RUN%.*}
    ./enqueue-single.sh $1 $RUN
done

echo 'Batch enqueueing ended'
