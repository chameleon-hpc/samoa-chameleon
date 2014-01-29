# Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
# Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
# This program is licensed under the GPL, for details see the file LICENSE


#!/bin/bash
echo "Plotting results..."

cd $1
echo "#Component breakdown" > "darcy.plt"

for file in darcy*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")

	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $(($processes * $threads)) $sections" " >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "integrity: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "load balancing: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Transport" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Gradient" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Permeability" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Pressure Solver" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Time step phase time" $file | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	grep -E "r0.*Cells[ ]+:[ ]+[0-9]+" $file | grep -oE "Cells[ ]+:[ ]+[0-9]+" | grep -oE "[0-9]+" | tr "\n" " " | cat >> "darcy.plt"
	echo "" >> "darcy.plt"
done

echo "#Component breakdown" > "swe.plt"

for file in swe*.log; do
	flags=$(echo $file | grep -oE "(_no[a-zA-Z0-9]+)+")
	processes=$(echo $file | grep -oE "_p[0-9]+" | grep -oE "[0-9]+")
	threads=$(echo $file | grep -oE "_t[0-9]+" | grep -oE "[0-9]+")
	sections=$(echo $file | grep -oE "_s[0-9]+" | grep -oE "[0-9]+")
	
	processes=${processes:-1}
	threads=${threads:-1}
	sections=${sections:-1}

	echo -n $(($processes * $threads)) $sections" "  >> "swe.plt"
	    
	grep -E "r0.*Adaptions" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "integrity: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Adaptions" $file | grep -oE "load balancing: [0-9]*\.[0-9]+" | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Time steps" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Displace" $file | grep -oE "(travs|time): [0-9]*\.?[0-9]+" | grep -oE "[0-9]*\.?[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Time Step phase time" $file | grep -oE "[0-9]*\.[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	grep -E "r0.*Cells[ ]+:[ ]+[0-9]+" $file | grep -oE "Cells[ ]+:[ ]+[0-9]+" | grep -oE "[0-9]+" | tr "\n" " " | cat >> "swe.plt"
	echo ""  >> "swe.plt"
done

sort -t" " -n -k 1,1 -k 2,2 -k 3,3 darcy.plt -o darcy.plt
sort -t" " -n -k 1,1 -k 2,2 -k 3,3 swe.plt -o swe.plt

#gnuplot &> /dev/null << EOT
gnuplot << EOT

set terminal postscript enhanced color font ',25'
set xlabel "Cores"
set key below font ",20" spacing 0.8 width -3
set xtics rotate
set yrange [0:*]
set auto x

set style line 999 lt 2 lw 8 lc rgb "black"

set for [n=1:64] style line n lt 1 lw 2
set style line 1 lc rgb "cyan"
set style line 2 lc rgb "orange"
set style line 3 lc rgb "magenta"
set style line 4 lc rgb "red"
set style line 5 lc rgb "blue"
set style line 6 lc rgb "green"
set style line 7 lc rgb "brown"
set style line 8 lc rgb "purple"

set style data histogram
set style histogram rowstacked
set style fill solid border 0
set boxwidth 0.75

#*******
# Darcy
#*******

set title "Darcy component breakdown"
set ylabel "Sec. per core (wall clock time)"
set output '| ps2pdf - darcy_components.pdf'
	
plot "darcy.plt" u (\$8) ls 2 t "Conformity", \
    '' u (\$6):xtic(1) ls 1 t "Adaption", \
	'' u (\$10) ls 3 t "Load Balancing", \
	'' u (\$14) ls 5 t "Gradient", \
	'' u (\$12) ls 4 t "Transport", \
	'' u (\$16) ls 6 t "Permeability", \
	'' u (\$20) ls 7 t "Pressure Solver"

set title "Darcy component breakdown - normalized"
set ylabel "Sec. per element (CPU time)"
set output '| ps2pdf - darcy_components_norm.pdf'
	
plot "darcy.plt" u (10.0 * \$20/\$19 * \$1/\$23) ls 7 t "Pressure Solver", \
	'' u (\$14/\$13 * \$1/\$23) ls 5 t "Gradient", \
	'' u (\$12/\$11 * \$1/\$23) ls 4 t "Transport", \
	'' u (\$16/\$15 * \$1/\$23) ls 6 t "Permeability", \
	'' u (\$8/\$5 * \$1/\$23) ls 2 t "Conformity", \
    '' u (\$6/\$5 * \$1/\$23):xtic(1) ls 1 t "Adaption", \
	'' u (\$10/\$5 * \$1/\$23) ls 3 t "Load Balancing"

#*****
# SWE
#*****

set title "SWE component breakdown"
set ylabel "Sec. per core (wall clock time)"
set output '| ps2pdf - swe_components.pdf'

plot    "swe.plt" u (\$12) ls 4 t "Time step", \
	    '' u (\$14) ls 5 t "Displace", \
	    '' u (\$8) ls 2 t "Conformity", \
        '' u (\$6):xtic(1) ls 1 t "Adaption", \
	    '' u (\$10) ls 3 t "Load Balancing"

set title "SWE component breakdown - normalized"
set ylabel "Sec. per element (CPU time)"
set output '| ps2pdf - swe_components_norm.pdf'
	
plot "swe.plt" u (\$12/\$11 * \$1/\$17) ls 4 t "Time step", \
    '' u (\$8/\$5 * \$1/\$17) ls 2 t "Conformity", \
    '' u (\$6/\$5 * \$1/\$17):xtic(1) ls 1 t "Adaption", \
	'' u (\$10/\$5 * \$1/\$17) ls 3 t "Load Balancing"
#	'' u (\$14/\$13 * \$1/\$17) ls 5 t "Displace", \


EOT
