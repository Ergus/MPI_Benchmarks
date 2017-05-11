#!/bin/bash

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 executable matrix_dim iterations local_threads alpha"
    exit 1
fi

exec=$1
dim=$2
iters=$3
lthreads=$4
alpha=$5

outdir="output_${exec}_${dim}_${iters}_$(date +%a%d%m%Y_%H%M)"
mkdir -p $outdir

statsfile="${outdir}/times.res"

cont=0

for i in {0..4}; do
    nodes=$((2**i))
    threads=$((nodes*lthreads))  #  threads= local*node
    for ((iter=0;iter<iters;iter++)); do
	format="${outdir}/${exec}_${nodes}_${lthreads}_${iter}"
	ofile="${format}.out"
	efile="${format}.err"

	echo -e "\nFormat: "$format
	
	command=$(sed -e "s/NODES/$nodes/g" -e "s/EXEC/$exec/g" -e "s/DIM/$dim/g" \
		      -e "s/THREADS/$threads/g" -e "s/PREFIX/${format/\//\\/}/g" \
		      -e "s/ALPHA/$alpha/g" submit.sh)
	
	echo "-> Number of nodes: $nodes ofile: $ofile efile: $efile"
	echo "$command" | bsub -o $ofile -e $efile
	((cont++))
    done
done

echo "Submited in total $cont jobs"
