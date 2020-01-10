#!/bin/bash

source './config.sh'

export CURRENT=${CURRENT:="./logs/$(date +%Y-%m-%d_%H-%M-%S)_render"}
mkdir -p $CURRENT

RENDER=$(basename -- "$1")
RENDER=${RENDER%.*}

export DATASRC=${DATASRC:='../../../output/verify'}

echo 'Starting rendering job'

for LIMITER in $LIMITERS; do
    echo "Running ${RENDER}_$LIMITER"
    "./renders/$RENDER.sh" $LIMITER 2>&1 | tee "$CURRENT/${RENDER}_${LIMITER}_render.log"
done

echo 'Rendering job ended'
