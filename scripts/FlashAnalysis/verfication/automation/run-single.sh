#!/bin/bash

source './machinesetup.sh'

RUN=$(basename -- "$2")
RUN=${RUN%.*}

echo 'Starting execution job'

for LIMITER in $LIMITERS; do
    OUTPUT="output/verify/${RUN}_$LIMITER"
    OUTPUTREL="$(dirname $0)/../../../../$OUTPUT"
    rm -rf $OUTPUTREL
    mkdir -p $OUTPUTREL
    echo "Running ${RUN}_$LIMITER"
    "$PROFILE/run.sh" $RUN $LIMITER $OUTPUT 2>&1 | tee "$CURRENT/${RUN}_${LIMITER}_run.log"
done

echo 'Execution job ended'
