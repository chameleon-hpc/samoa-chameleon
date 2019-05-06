#!/bin/bash

rm table.csv
echo "chameleon ranks threads lbfreq sections dmin dmax time" > table.csv

for f in *.log
do
  python parse_samoa_log.py $f >> table.csv 
done
