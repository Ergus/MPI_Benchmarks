#!/bin/bash

stats_file="times.res"
echo "Nodes\t time(msec)" > $stats_file

rm -f *.mat outfile_matmul_nodes_*.{out,err}

module load intel python/3.5.0

for i in {0..4}; do
    nodes=$((2**i))
    rm -f *.mat
    ofile="out_mpi_matmul_nodes_$nodes.out"
    efile="out_mpi_matmul_nodes_$nodes.err"
    echo "-> Number of nodes: $nodes ofile: $ofile efile: $efile"
    bsub -K -n $nodes -o $ofile -e $efile < submit-matmul.sh
    
    ./prove.py matrix_A.mat matrix_B.mat matrix_C.mat
    if [[ $? -ne 0 ]]; then
	echo "Test with $i nodes failed"
	exit 1
    fi
done

echo "Filtering results from output files"
awk '/algorithm_time/{match(FILENAME,/([[:digit:]]+)/,arr); printf(%s %.3f, arr[1], /10^6 ) }' *.out | sort -n > data.csv

echo "Building the graphs"
module load GNUPLOT/5.0.1
gnuplot ./gnuplot.plt
