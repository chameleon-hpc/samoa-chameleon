#!/bin/bash

cd ..

OUTPUT="results/resting_lake/$1"
OUTPUT2="results/resting_lake/$1"
mkdir -p $OUTPUT $OUTPUT2
./resting_lake.py -s island -i `find "$DATASRC/resting_lake_$1" -type f -name "*.xmf"` -o $OUTPUT -n 8
./resting_lake.py -s overlapping -i `find "$DATASRC/resting_lake2_$1" -type f -name "*.xmf"` -o $OUTPUT2 -n 8
