#!/usr/bin/env bash
#SBATCH --job-name=mpi_dynamic_delete
#SBATCH --qos=debug
#SBATCH --time=00-00:05:00
#SBATCH --error=%x_%A.err 
#SBATCH --output=%x_%A.out
#SBATCH --nodes=16
#SBATCH --tasks-per-node=1

for i in 1 2 4 8 16; do
    times=$((SLURM_JOB_NUM_NODES/i))     # Estimates nproc/step
    char='%.0s -d '$i                    # generates a substring '%0s -s step'
    input=$(printf "${char}" $(eval echo {1..$times}))
    echo "Running times="$times
    echo $input
    mpirun -np 1 ./dynamicOOP -s 16 ${input} -i | tee ${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_JOB_NUM_NODES}_$i.out
done
