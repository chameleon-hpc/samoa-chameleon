#!/bin/bash

source './config.sh'

export CURRENT=${CURRENT:="./logs/$(date +%Y-%m-%d_%H-%M-%S)_render"}
mkdir -p $CURRENT

echo 'Starting batch rendering'

for RENDER in renders/*; do
    RENDER=$(basename -- "$RENDER")
    RENDER=${RENDER%.*}
    ./render-single.sh $RENDER
done

echo 'Batch rendering ended'
