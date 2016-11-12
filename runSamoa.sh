#!/bin/bash

rm -r output
mkdir output
export OMP_NUM_THREADS=8

# n=1

# #x=`echo "scale=16;1/(${n}+1)"| bc | awk '{printf "%f", $0}'`
# #y=`echo "scale=16;1/(${n}+1)"| bc | awk '{printf "%f", $0}'`

# x=`echo "scale=16;1/(${n}+1)"| bc `
# y=`echo "scale=16;1/(${n}+1)"| bc `


# nodes="${x} ${y}"

# for i in `seq 1 ${n}`;
# do
#     for j in `seq 1 ${n}`;
#     do
# 	echo "scale=16;1.0/(${n}+1.0)*${i}"| bc
# #	x=`echo "scale=16;1.0/(${n}+1.0)*${i}"| bc | awk '{printf "%f", $0}'`
# #	y=`echo "scale=16;1.0/(${n}+1.0)*${j}"| bc | awk '{printf "%f", $0}'`
# 	x=`echo "scale=16;1.0/(${n}+1.0)*${i}"| bc `
# 	y=`echo "scale=16;1.0/(${n}+1.0)*${j}"| bc `

# 	nodes="${nodes},${x} ${y}"
#     done
# done

nodes='0.166667 0.166667,0.166667 0.333333,0.166667 0.500000,0.166667 0.666667,0.166667 0.833333,0.333333 0.166667,0.333333 0.333333,0.333333 0.500000,0.333333 0.666667,0.333333 0.833333,0.500000 0.166667,0.500000 0.333333,0.500000 0.500000,0.500000 0.666667,0.500000 0.833333,0.666667 0.166667,0.666667 0.333333,0.666667 0.500000,0.666667 0.666667,0.666667 0.833333,0.833333 0.166667,0.833333 0.333333,0.833333 0.500000,0.833333 0.666667,0.833333 0.833333'

for i in `seq 9 10`;
do

mkdir output/output_1_unlimited_$i
bin/samoa_1_unlimited -tmax 0.05 -dmin $i -dmax $i -tout 0.0 -output_dir output/output_1_unlimited_$i -drytolerance 0.0001 -stestpoints ${nodes} -xmloutput
done