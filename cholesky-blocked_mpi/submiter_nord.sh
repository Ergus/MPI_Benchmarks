#!/bin/bash

if [ $# -ge 1 ] && [ -f $1 ]; then
   echo "Adding Jobs for: "${1}
else
	echo "No file: "$1
	echo "Usage $0 cholesky_executable [tasks]"
	exit
fi

[ $# -ge 2 ] && tasks=$2 || tasks=1

nodes=(1 2 4 8 16 32)
blocksize=(64 128 256 512)
dims=(1024 2048 4096 8192 16384)

# nodes=(1 2)
# blocksize=(64)
# dims=(1024)

now=$(date +%F_%T)
resdir="results_${now}"

mkdir ${resdir}

echo "Submitting group"
echo "nodes ${nodes[@]}"
echo "blocks ${blocksize[@]}"
echo "dim ${dims[@]}"

for dim in ${dims[@]}; do
	for bs in ${blocksize[@]}; do
		for node in ${nodes[@]}; do
			for ((it=0; it<$tasks; ++it)); do
				if [ $((node*bs<=dim)) = 1 ]; then
					printf "Submitting combination dim: %s bs: %s nodes: %s\n" \
						   $dim $bs $node
					jobname="mpi_${dim}_${bs}_${node}"

					sed -e "s|__EXEC|$1|" \
						-e "s|__DIM|${dims}|" \
						-e "s|__BS|${bs}|" \
						submit_nord.sh | bsub -W 10 \
										 -J "${jobname}" \
										 -o "${resdir}/${jobname}_${it}_%J.out" \
										 -e "${resdir}/${jobname}_${it}_%J.err" \
										 -n ${node}

				else
					printf "Jump combination dim: %s bs: %s nodes: %s\n" \
						   $dim $bs $node
				fi
			done
		done
	done
done
