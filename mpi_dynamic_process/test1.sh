#!/usr/bin/env bash
#SBATCH --job-name=mpi_dynamic
#SBATCH --qos=debug
#SBATCH --time=00-00:05:00
#SBATCH --error=%x_%A.err 
#SBATCH --output=%x_%A.out
#SBATCH --nodes=4
#SBATCH --tasks-per-node=1

mpirun -np 1 ./dynamicOOP -s 1 -s 1 -s 1 -i -s 1 -s 1 -s 1 -i
mpirun -np 1 ./dynamicOOP -s 1 -s 1 -s 1 -i -s 1 -s 1 -s 1 -i
mpirun -np 1 ./dynamicOOP -s 1 -s 1 -s 1 -i -s 1 -s 1 -s 1 -i
mpirun -np 1 ./dynamicOOP -s 1 -s 1 -s 1 -i -s 1 -s 1 -s 1 -i
