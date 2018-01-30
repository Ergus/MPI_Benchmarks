#!/usr/bin/gnuplot -c

if (ARGC < 2) {
	print sprintf("Usage: %s spawn_output shrink_output", ARG0)
	exit
}

# Basic awk filtering
spawn1="<awk '/^World/ && $2<$4 && $2==%s { print $7, $15 }' %s"
spawn2="<awk '/^World/ && $2<$4 && $7==%s { print $2, $15 }' %s"

shrink1="<awk '/^World/ && $2>$4 && $7==%s { print $2, $15 }' %s"
shrink2="<awk '/^World/ && $2>$4 && $2==%s { print $7, $15 }' %s"

# This is the average receives 2 columns 
mean=" | awk '{sum[$1]+=$2; var[$1][cont[$1]+=1]=$2} \
			END{ \
				for (s in var) { \
                	n = cont[s]; \
					m = sum[s]/n; \
					std = 0.0; \
					for (i in var[s]) { \
						std += (var[s][i]-m)^2; \
						}; \
                    e = sqrt(std / (n*(n-1)) ); \
                    print s, m, e; \
				} \
			}' "

########## Global sets ########################
set style line  1 lt 1 pt 7 lw 4 ps 0.75
set term postscript eps noenhanced color linewidth 0.5 dashlength 4
set grid

set key left

set ylabel "Time(ns)"

colors ="#0072bd #d95319 #edb120 #7e2f8e #77ac30 #4dbeee #a2142f"

########## Graph 1 (spawn size) ################

set output "| epstopdf --filter --outfile=spawn_size.pdf"

set xlabel "Spawn size"

vals="1 2 4 8 17 31 39"

plot for [i=1:words(vals)] sprintf(spawn1.mean, word(vals,i), ARG1) u 1:2:3 w errorbars ls 1 lc rgb word(colors,i) notitle, \
     for [i=1:words(vals)] sprintf(spawn1.mean, word(vals,i), ARG1) u 1:2 w lines ls 1 lc rgb word(colors,i) title "Initial Nodes: ".word(vals,i)

########## Graph 2 (spawn initial) ###############

set output "| epstopdf --filter --outfile=spawn_initial.pdf"

set xlabel "Initial Nodes"

vals="1 2 4 8 16 32 40"

plot for [i=1:words(vals)] sprintf(spawn2.mean, word(vals,i), ARG1) u 1:2:3 w errorbars ls 1 lc rgb word(colors,i) notitle, \
     for [i=1:words(vals)] sprintf(spawn2.mean, word(vals,i), ARG1) u 1:2 w lines ls 1 lc rgb word(colors,i) title "Spawn Size: ".word(vals,i)

########## Graph 3 (shrink size) #################

set output "| epstopdf --filter --outfile=shrink_size.pdf"
set key right bottom

set xlabel "Shrink size"

vals="2 3 5 9 17 33"

plot for [i=1:words(vals)] sprintf(shrink2.mean, word(vals,i), ARG2) u 1:2:3 w errorbars ls 1 lc rgb word(colors,i) notitle, \
     for [i=1:words(vals)] sprintf(shrink2.mean, word(vals,i), ARG2) u 1:2 w lines ls 1 lc rgb word(colors,i) title "Initial Nodes: ".word(vals,i)


########## Graph 4 (shrink initial) ##############

set output "| epstopdf --filter --outfile=shrink_initial.pdf"

set xlabel "Initial Nodes"

vals="1 2 4 8 16 32"

plot for [i=1:words(vals)] sprintf(shrink1.mean, word(vals,i), ARG2) u 1:2:3 w errorbars ls 1 lc rgb word(colors,i) notitle, \
     for [i=1:words(vals)] sprintf(shrink1.mean, word(vals,i), ARG2) u 1:2 w lines ls 1 lc rgb word(colors,i) title "Shrink Size: ".word(vals,i)

### Exit ###
