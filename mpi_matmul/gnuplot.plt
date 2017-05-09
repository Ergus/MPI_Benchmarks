#!/bin/env gnuplot

filter="<awk 'FNR==1{dim=$5; th=$7; n=$9}\
/algorithm_time/{time=$5/10^6} \
/Arrays match/{print dim, th, n, time; nextfile}' %s/*.out | \
sort -n -k3 | \
awk \'BEGIN{print \"Nodes\", \"time\"} \
{if(\$3!=N){ \
    if(NR>1) \
        {print N, mean/cont};\
    N=\$3; cont=0; mean=0}}\
{cont++;mean+=\$4}\
END{print N, mean/cont}'"

set term postscript eps noenhanced color linewidth 0.5 dashlength 4

#set output '| epstopdf --filter --outfile=Time-Speedup.pdf'
set output dir.".eps"

set multiplot layout 1, 2 title "MPI_Matmult"
set grid
set xlabel "Nodes"
set ylabel "Time(ms)"
plot sprintf(filter,dir) u 1:2 w linespoints pt 5 ps 0.5 lt 1 notitle

command="<awk '$1==1{var=$2}{print $1, var/$2}' data.csv"

command2=filter." | awk 'NR==1{next}$1==1{var=$2}{print $1, var/$2}'"

print sprintf(command2,dir)

set ylabel "Speedup"
plot sprintf(command2,dir) u 1:2 w linespoints pt 5 ps 0.5 lt 1 notitle