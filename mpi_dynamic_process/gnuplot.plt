
#filter_spawn="awk '/^World/{print $2,$4,$7,$15,$19,$21}' mpi_dynamic_spawn_76445.out"
#filter_delete="awk '/^World/ && $2>$4 && $2==17 {print $2,$4,$7,$13,$17}' mpi_dynamic_delete_76446.out"

#mean="awk '{val[$3][cont[$3]++]=$5; sum[$3]+=$5}END{for(r in val) {avg[r]=sum[r]/cont[r]; for (v in val[r]){ std[r]+=(val[r][v]-avg[r])^2 } printf(\"%f %f\n\",avg[r],sqrt(std[r]/(cont[r]-1)))}}'"

set style line  1 lt 1 pt 7 lw 0.75 ps 1
set style line  2 lt 1 pt 6 lw 0.75 ps 1

colors ="#0072bd #d95319 #edb120 #7e2f8e #77ac30 #4dbeee #a2142f"
tit ="Split Total"
col ="13 17"
initial ="9 17"

set term postscript eps enhanced color linewidth 2 dashlength 4

set xlabel "Deleted Nodes"
set ylabel "Time (ns)"

mean_delete="<awk -v col=%s '/^World/ && $2>$4 && $2==%s {val[$7][cont[$7]++]=$col; sum[$7]+=$col}END{for(r in val) {avg[r]=sum[r]/cont[r]; for (v in val[r]){ std[r]+=(val[r][v]-avg[r])^2 } print r,avg[r],sqrt(std[r]/(cont[r]-1))}}' %s"

do for [it=1:2] {
  init=word(initial,it)
  set title sprintf("Split-Reduce initial %s processes",init)
  set output sprintf('|epstopdf --filter --outfile=Split_Reduce-%s.pdf',init)
  plot for [i=2:1:-1] sprintf(mean_delete,word(col,i),init,"mpi_dynamic_delete_76446.out") using 1:2 with filledcurve y1=0 lc rgb word(colors,i) title word(tit,i), \
       for [i=2:1:-1] sprintf(mean_delete,word(col,i),init,"mpi_dynamic_delete_76446.out") using 1:2:3 with errorbars pt 7 lc black notitle
}
