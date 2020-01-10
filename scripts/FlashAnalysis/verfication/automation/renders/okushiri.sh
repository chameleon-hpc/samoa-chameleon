#!/bin/bash

cd ..

OUTPUT="results/okushiri/$1"
mkdir -p $OUTPUT
./okushiri.py -i `find "$DATASRC/okushiri_$1" -type f -name "*.xmf"` -s surface -o $OUTPUT
./okushiri.py -i `find "$DATASRC/okushiri_$1" -type f -name "*.xmf"` -d "./reference/okushiri" -s gauges -o $OUTPUT -n 8
