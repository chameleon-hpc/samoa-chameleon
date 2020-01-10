#!/bin/bash

cd ..

OUTPUT="results/longwave_basin/dw10e-2/$1"
OUTPUT2="results/longwave_basin/dw10e-8/$1"
mkdir -p $OUTPUT $OUTPUT2
./longwave_basin.py -i `find "$DATASRC/longwave_basin-dw10e-2-light_$1" -type f -name "*.xmf"` -o $OUTPUT
./longwave_basin.py -s 0.01 -i `find "$DATASRC/longwave_basin-dw10e-2-light_$1" -type f -name "*.xmf"` -o $OUTPUT
./longwave_basin.py -i `find "$DATASRC/longwave_basin-dw10e-8-light_$1" -type f -name "*.xmf"` -o $OUTPUT2
./longwave_basin.py -s 0.00000001 -i `find "$DATASRC/longwave_basin-dw10e-8-light_$1" -type f -name "*.xmf"` -o $OUTPUT2
