#!/bin/bash

#SBATCH --qos=bsc_cs
#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

source ../argparse.sh
add_argument -a x -l exe -h "Executable file"
add_argument -a d -l dim -h "Dimension"
add_argument -a b -l bsize -h "Block size"
add_argument -a i -l iter -h "Repetitions default[1]" -d 1
add_argument -a e -l extrae -h "Use extrae? default false" -b
parse_args "$@"
printargs >&2

module purge
module load bsc/1.0
module load gcc/7.2.0
module load liberep/git
module load mkl/2018.1
module load hwloc/1.11.8
module load impi/2018.1

echo -e "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo -e "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo -e "# Nodes_List: ${SLURM_JOB_NODELIST} QOS: ${SLURM_JOB_QOS}"
echo -e "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo -e "# Command: ${ARGS[x]} ${ARGS[d]} ${ARGS[b]} 0"
echo -e "# --------------------------------------\n"

for ((it=0; it<ARGS[i]; ++it)) {
		echo "# Starting it: ${it} at: " $(date)
		start=$SECONDS
		srun ${ARGS[x]} ${ARGS[d]} ${ARGS[b]} 0
		end=$SECONDS
		echo "# Ending: " $(date)
		echo "# Elapsed: "$((end-start))
		echo -e "\n# --------------------------------------\n"
	}
