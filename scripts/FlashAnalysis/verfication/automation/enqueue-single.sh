#!/bin/bash

source './machinesetup.sh'

RUN=$(basename -- "$2")
RUN=${RUN%.*}

echo 'Starting enqueuing job'

for LIMITER in $LIMITERS; do
    echo "Enqueueing ${RUN}_$LIMITER"
    "$PROFILE/enqueue.sh" $RUN $LIMITER 2>&1 | tee "$CURRENT/${RUN}_${LIMITER}_enqueue.log"
done

echo 'Enqueueing job ended'
