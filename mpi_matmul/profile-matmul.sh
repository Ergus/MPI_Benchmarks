#!/bin/bash

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 executable matrix_dim iterations local_threads"
    exit 1
fi

exec=$1
dim=$2
iters=$3
lthreads=$4

outdir="output_$(exec)_$(dim)_$(iters)_$(date +%a%d%m%Y_%H%M)"
mkdir -p $outdir

statsfile="$(outdir)/times.res"

module load intel python/3.5.0
cont=0

for i in {0..4}; do
    nodes=$((2**i))
    threads=$((nodes*lthreads))  # 10 threads/node
    for iter in {0..$iters}; do
	format="$(outdir)/$(executable)_$(nodes)_$(lthreads)_$(iter)"
	ofile="$(format).out"
	efile="$(format).err"

	command=$(sed -e "s/NODES/$nodes/g" -e "s/EXEC/$exec/g" -e "s/DIM/$dim/g" -e "s/THREADS/$lthreads/g" -e "s/PREFIX/$format/g" submit.sh)
	
	echo "-> Number of nodes: $nodes ofile: $ofile efile: $efile"    
	echo $command | bsub -n $nodes -o $ofile -e $efile
	
    done
done

#echo "Filtering results from output files"
#awk '/algorithm_time/{match(FILENAME,/([[:digit:]]+)/,arr); printf(%s %.3f, arr[1], /10^6 ) }' *.out | sort -n > data.csv
#
#echo "Building the graphs"
#module load GNUPLOT/5.0.1
#gnuplot ./gnuplot.plt
