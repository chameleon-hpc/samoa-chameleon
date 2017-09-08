#!/usr/bin/env bash
set -o errexit
set -o nounset
unset CDPATH

cores=$(grep -c '^processor' '/proc/cpuinfo') #Includes hyperthreading. Which should be exactly what we want.

# Compile dg for all orders.
for (( order=1; order <= 2; order+=1)) do
    scons config=my_conf.py swe_dg_order=${order} target=release -j ${cores} swe_scenario="linear_dam_break" exe="samoa_dg_${order}_unlimited_dam" 
    scons config=my_conf.py swe_dg_order=${order} target=release -j ${cores} swe_scenario="oscillating_lake" exe="samoa_dg_${order}_unlimited_lake" 
done
