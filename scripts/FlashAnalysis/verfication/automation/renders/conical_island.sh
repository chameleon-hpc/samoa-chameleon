#!/bin/bash

cd ..

OUTPUTA="results/conical_island/a/$1"
OUTPUTC="results/conical_island/c/$1"
mkdir -p $OUTPUTA $OUTPUTC
./conical_island.py -i `find "$DATASRC/conical_island_a_$1" -type f -name "*.xmf"` -s surface_a -o $OUTPUTA
./conical_island.py -i `find "$DATASRC/conical_island_a-light_$1" -type f -name "*.xmf"` -d "./reference/conical_island" -s gauges_a -o $OUTPUTA -n 8
./conical_island.py -i `find "$DATASRC/conical_island_c_$1" -type f -name "*.xmf"` -s surface_c -o $OUTPUTC
./conical_island.py -i `find "$DATASRC/conical_island_c-light_$1" -type f -name "*.xmf"` -d "./reference/conical_island" -s gauges_c -o $OUTPUTC -n 8
