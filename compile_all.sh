#!/usr/bin/env bash
set -o errexit
set -o nounset
unset CDPATH

cores=$(grep -c '^processor' '/proc/cpuinfo') #Includes hyperthreading. Which should be exactly what we want.

# Compile dg for all orders.
for (( order=1; order <= 4; order+=1)) do
    scons config=my_conf.py swe_dg_order=${order} -j ${cores}
done
