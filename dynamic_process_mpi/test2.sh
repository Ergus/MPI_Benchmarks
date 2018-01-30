#!/usr/bin/env bash
#SBATCH --job-name=dynamic_delete
#SBATCH --qos=debug
#SBATCH --time=00-00:10:00
#SBATCH --nodes=16
#SBATCH --error=%x_%A.err 
#SBATCH --output=%x_%A.out
#SBATCH --tasks-per-node=1

if (($# < 1)); then
    (>&2 echo "Error missing parameters")
    exit
fi

for i in 1 2 4 8 16; do
    END=$1
    times=$((SLURM_JOB_NUM_NODES/i))     # Estimates nproc/step
    char='%.0s -d '$i                    # generates a substring '%0s -s step'
    input=$(printf "${char}" $(eval echo {1..$times}))
    echo "Running times="$times
    echo $input
    for ((a=1; a<=END; ++a)); do
	mpirun -np 1 ./dynamicOOP -s 16 -i ${input} | tee ${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_JOB_NUM_NODES}_$i.out
    done
done
