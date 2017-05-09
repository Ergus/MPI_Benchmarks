#!/bin/env gnuplot

set term postscript eps noenhanced color linewidth 0.5 dashlength 4

#set output '| epstopdf --filter --outfile=Time-Speedup.pdf'
set output "MPI_Matmult.eps"

set multiplot layout 1, 2 title "MPI_Matmult"
set grid
set xlabel "Nodes"
set ylabel "Time(ms)"
plot "data.csv" w linespoints pt 5 ps 0.5 lt 1 notitle

command="<awk '$1==1{var=$2}{print $1, var/$2}' data.csv"

set ylabel "Speedup"
plot sprintf(command) w linespoints pt 5 ps 0.5 lt 1 notitle