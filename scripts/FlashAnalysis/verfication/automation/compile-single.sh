#!/bin/bash

source './machinesetup.sh'

SCENARIO=$(basename -- "$2")
SCENARIO=${SCENARIO%.*}

echo 'Starting compilation job'

for LIMITER in $LIMITERS; do
    echo "Compiling ${SCENARIO}_$LIMITER"
    source "builds/$SCENARIO.sh"
    "$PROFILE/compile.sh" $SCENARIO $LIMITER 2>&1 | tee "$CURRENT/${SCENARIO}_${LIMITER}_build.log"
done

echo 'Compilation job ended'
