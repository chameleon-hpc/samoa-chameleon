#!/bin/bash

cd ..

OUTPUT="results/oscillating_lake/$1"
mkdir -p $OUTPUT
./oscillating_lake.py -i `find "$DATASRC/oscillating_lake_$1" -type f -name "*.xmf"` -o $OUTPUT
./oscillating_lake.py -s 0.001 -i `find "$DATASRC/oscillating_lake_$1" -type f -name "*.xmf"` -o $OUTPUT
./oscillating_lake.py -s "mass" -i `find "$DATASRC/oscillating_lake_$1" -type f -name "*.xmf"` -o $OUTPUT -n 8
./oscillating_lake.py -s "energy" -i `find "$DATASRC/oscillating_lake_$1" -type f -name "*.xmf"` -o $OUTPUT -n 8
./oscillating_lake.py -s "series_$1" -i $DATASRC -o $OUTPUT

