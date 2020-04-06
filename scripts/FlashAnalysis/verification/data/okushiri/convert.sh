#!/bin/bash

./csv2nc.py -i bathymetry.csv -o bathymetry.nc --scale=-1.0 --sample=2
./csv2nc.py -i displacement.csv -o displacement.nc --scale=1.0 --sample=0