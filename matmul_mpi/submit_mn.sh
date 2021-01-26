#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

source @PROJECT_BINARY_DIR@/argparse.sh
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a R -l repeats -h "Repetitions per program default[1]" -t int

parse_args "$@"
printargs >&2

dims=(1024 4096 8192)
blocksizes=(64 128 256)

REPEATS="${ARGS[R]}"
COMMAND="${ARGS[x]} ${ARGS[d]} ${ARGS[b]} 5"

# Start run here printing run info header
echo "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo "# QOS: ${SLURM_JOB_QOS}"
echo "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"

echo "# Repetitions: ${REPEATS}"
echo "# Command: ${COMMAND}"
env | grep NANOS6 | sed -e 's/^#*/# /'
echo "# ======================================"

for dim in ${dims[@]}; do
	for bs in ${blocksizes[@]}; do
		if [ $((SLURM_JOB_NUM_NODES*bs<=dim)) = 1 ]; then
			for ((it=0; it<${REPEATS}; ++it)) {
				echo "# Starting it: ${it} at: $(date)"
				start=${SECONDS}
				srun ${ARGS[x]} $dim ${ARGS[b]} 5
				end=${SECONDS}
				echo "# Ending:  $(date)"
				echo "# Elapsed: $((end-start))"
				echo "# --------------------------------------"
			}
		else
			echo "# Jump combination nodes: $node, dim: $dim, bs: $bs"
		fi
	done
done
