filter="awk 'FNR==1{dim=$5; th=$7; n=$9} \
            /algorithm_time/{time=$5/10^6} \
            /Arrays match/{print dim, th, n, time; nextfile}' %s/*.out \
            | sort -n -k3 > %s_general.csv "

averag="awk 'function avg(ar,n,mean,  r){ \
               for(i=0;i<n;++i){ \
                 r+=(ar[i]-mean)^2 } \
               return sqrt(r/(n*(n-1))) } \
            BEGIN{print \"Nodes\", \"Time\", \"Error\"} \
            {if($3!=N){ \
               if(cont>0) \
                 {print N, mean/=cont, avg(vals,cont,mean)}; \
               N=$3; cont=0; mean=0}} \
            {vals[cont++]=$4; mean+=$4} \
            END{print N, mean/=cont, avg(vals,cont,mean)}' %s_general.csv \
            | tee %s_reduced.csv "

speedup="awk 'NR==1{next}$1==1{var=$2;err=$3} \
              {print $1, var/$2, $2*err+var*$3}' %s_reduced.csv "

splitted(in)=word(system(sprintf("echo %s | sed -e \"s|_| |g\"",in)),2)

system(sprintf(filter,dir,dir))

set term postscript eps noenhanced color linewidth 0.5 dashlength 4
set output dir.".eps

# line styles
set style line  1 lt 1 pt 7 lw 4 ps 0.75
set style line  2 lt 1 pt 6 lw 0.75 ps 1

colors ="#0072bd #d95319 #edb120 #7e2f8e #77ac30 #4dbeee #a2142f"

set multiplot layout 1, 2 title dir
#set grid
set xlabel "Nodes"
set ylabel "Time(ms)"
plot sprintf("<".averag,dir,dir) u 1:2:3 w errorbars ls 1 notitle, \
     '' u 1:2 w lines ls 1 title splitted(dir)

set ylabel "Speedup"
plot sprintf("<".speedup,dir) u 1:2 w lp ls 1 title splitted(dir)

print "Processed folder: ".dir