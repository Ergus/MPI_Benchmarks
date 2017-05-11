filter="awk 'FNR==1{dim=$5; th=$7; n=$9}\
            /algorithm_time/{time=$5/10^6} \
            /Arrays match/{print dim, th, n, time; nextfile}' %s/*.out \
            | sort -n -k3 > %s_general.csv "
averag="awk 'BEGIN{print \"Nodes\", \"time\"}\
            {if(\$3!=N){ \
               if(NR>1) \
                 {print N, mean/cont};\
               N=\$3; cont=0; mean=0}}\
            {cont++;mean+=$4}\
            END{print N, mean/cont}' %s_general.csv \
            | tee %s_reduced.csv "

speedup="awk 'NR==1{next}$1==1{var=$2}{print $1, var/$2}' %s_reduced.csv "

splitted(in)=word(system(sprintf("echo %s | sed -e \"s|_| |g\"",in)),2)

system(sprintf(filter,dir,dir))

filename=dir.".eps"

set style line 1 pt 5 ps 0.5 lt 1 lc rgb 'red'

set term postscript eps noenhanced color linewidth 0.5 dashlength 4

set output filename

set multiplot layout 1, 2 title dir
set grid
set xlabel "Nodes"
set ylabel "Time(ms)"
plot sprintf("<".averag, dir, dir) u 1:2 w lp ls 1 title splitted(dir)

set ylabel "Speedup"
plot sprintf("<".speedup,dir) u 1:2 w lp ls 1 title splitted(dir)

print "Graph saved to: ".filename