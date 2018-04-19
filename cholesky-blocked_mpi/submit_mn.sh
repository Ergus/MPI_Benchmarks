#!/bin/bash

#SBATCH --qos=bsc_cs
#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

module purge
module load bsc/1.0
module load gcc/7.2.0
module load liberep/git
module load mkl/2018.1
module load hwloc/1.11.8
module load impi/2018.1

echo "# Starting: " $(date) >&2
start=$SECONDS
srun $@
end=$SECONDS
echo "# Ending: " $(date) >&2
echo "# Elapsed: "$((end-start)) >&2
