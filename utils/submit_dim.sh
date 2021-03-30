#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

# Declare command line arguments.
# In this script there are not default values to prevent errors.
source @PROJECT_BINARY_DIR@/argparse.sh

# Arguments to use in this script
add_argument -a R -l repeats -h "Repetitions per program" -t int
add_argument -a W -l weakscaling -h "Do weak scaling (dim*sqrt(nodes))" -t int

# Arguments for the executable
add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations" -t int

# Parse input command line arguments
parse_args "$@"

init=${SECONDS}

DIM=${ARGS[D]}
BS=${ARGS[B]}
REPEATS=${ARGS[R]}
ITS=${ARGS[I]}

# Start run here printing run info header
echo "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo "# QOS: ${SLURM_JOB_QOS}"
echo "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo "# Walltime: $(squeue -h -j $SLURM_JOBID -o "%l")"

# Print command line arguments
printargs "# "

# Print nanos6 environment variables
env | grep NANOS6 | sed -e 's/^#*/# /'
echo "# ======================================"

if [ $((SLURM_JOB_NUM_NODES*BS<=DIM)) != 1 ]; then
	echo "# Jump combination nodes: $node, dim: $rows bs: $BS"
	exit
fi

for exe in @TEST@_*; do

	COMMAND="$exe $DIM $BS $ITS"
	echo -e "# Starting command: ${COMMAND}"
	echo "# ======================================"
	for ((it=0; it<${REPEATS}; ++it)) {
		echo "# Starting it: ${it} at: $(date)"
		start=${SECONDS}
		srun ./${COMMAND}
		end=${SECONDS}
		echo "# Ending: $(date)"
		echo "# Elapsed: $((end-start))"
		echo "# --------------------------------------"
	}
		echo ""
done

# We arrive here only when not wall time was reached.
# else the function job_finish_hook is called.
finalize=${SECONDS}
echo "# Done:  $(date)"
echo "# Elapsed: $((finalize-init))"
echo "# ======================================"
