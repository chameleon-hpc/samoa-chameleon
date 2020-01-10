#!/bin/bash

echo 'Starting pack results job'

cd '../../../../output/verify'
mkdir -p '../verify-pack'

for RESULT in ./*; do
    RESULTE=$(basename -- "$RESULT")
    echo "Packing result $RESULTE"
    tar -zcvf "../verify-pack/$RESULTE.tar.gz" $RESULT
done

echo 'Pack results job ended'
