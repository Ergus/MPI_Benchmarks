#!/usr/bin/env bash
#SBATCH --job-name=spawn
#SBATCH --qos=bsc_cs
#SBATCH --time=8:00:00
#SBATCH --nodes=40
#SBATCH --error=%x_%A.err 
#SBATCH --output=%x_%A.out
#SBATCH --tasks-per-node=1

#----------- Header --------------
echo -e "# Job: ${SLURM_JOB_NAME}\tid: ${SLURM_JOB_ID}"
echo -e "# Total Nodes: ${SLURM_JOB_NUM_NODES}\tTasks per Node: ${SLURM_NTASKS_PER_NODE}\tCPUS in master: ${SLURM_JOB_CPUS_PER_NODE}"
echo -e "# Nodes List: ${SLURM_JOB_NODELIST}\tQOS: ${SLURM_JOB_QOS}"
echo -e "# Account: ${SLURM_JOB_ACCOUNT}\tSubmitter host: ${SLURM_SUBMIT_HOST}\tRunning Host: ${SLURMD_NODENAME}"
echo -e "# Start: $(date)"
T=$(date +%s)
echo -n "# " && printf '%.0s-' {1..50} && echo # prints a line
#------------ Main ----------------

if (($# < 1)); then
    (>&2 echo "Error missing parameters")
    exit
fi

readonly END=$1
for i in {1..8} {9..39..2}; do
    for j in {1..9} {10..40..2}; do
	for ((a=1; a<=END; ++a)); do
	    mpirun -np 1 ./dynamicOOP -s $i -s $j
	done
    done
done

#------------- Foot ----------------
echo -n "# " && printf '%.0s-' {1..50} && echo # prints a line
echo -e "# Job: "${SLURM_JOB_NAME}"\tid: "${SLURM_JOB_ID}"\tend: "$(date)
T=$(($(date +%s)-T))
printf "# Elapsed time: %02d:%02d:%02d:%02d\n" "$((T/86400))" "$((T/3600%24))" "$((T/60%60))" "$((T%60))"
echo "# Time in seconds: ${T}"
echo -n "# " && printf '%.0s-' {1..50} && echo # prints a line
