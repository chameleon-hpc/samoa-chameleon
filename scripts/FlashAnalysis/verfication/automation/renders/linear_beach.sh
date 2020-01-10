#!/bin/bash

cd ..

OUTPUT="results/linear_beach/$1"
mkdir -p $OUTPUT
./linear_beach.py -i `find "$DATASRC/linear_beach-light_$1" -type f -name "*.xmf"` -d "./reference/linear_beach" -o $OUTPUT
