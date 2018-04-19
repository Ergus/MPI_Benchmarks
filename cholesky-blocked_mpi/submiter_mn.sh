#!/bin/bash

source parser.sh

add_argument -a x -l exe -h "Executable file"
add_argument -a i -l iter -h "Repetitions default[1]" -d 1
add_argument -a e -l extrae -h "Use extrae? default false" -b
parse_args "$@"

printargs

nodes=(1 2 4 8 16 32)
blocksize=(64 256 512)
dims=(1024 4096 8192 16384)

# nodes=(1 2)
# blocksize=(64)
# dims=(1024)

now=$(date +%F_%H-%M-%S)
name=${ARGS[x]/%.nanos6}
resdir="results/${name}_${now}"

mkdir -p ${resdir}

echo "Directory ${resdir}"
echo "Submitting group"
echo "nodes ${nodes[@]}"
echo "blocks ${blocksize[@]}"
echo "dim ${dims[@]}"

for dim in ${dims[@]}; do
	for node in ${nodes[@]}; do
		for bs in ${blocksize[@]}; do

			if [ $((node*bs<=dim)) = 1 ]; then
				printf "Submitting combination dim: %s bs: %s nodes: %s\n" \
					   $dim $bs $node
				jobname="mpi_${dim}_${bs}_${node}"
				filename="${resdir}/${jobname}"

 				sbatch --ntasks=${node} \
					   --time=01:00:00 \
					   --array=1-${ARGS[i]}  \
 					   --job-name=${jobname} \
 					   --output="${resdir}/%x_%2a_%j.out" \
 					   --error="${resdir}/%x_%2a_%j.err" \
 					   ./submit_mn.sh ${ARGS[x]} ${dim} ${bs} 0
			else
				printf "Jump combination dim: %s bs: %s nodes: %s\n" \
					   $dim $bs $node
			fi
		done
	done
done

