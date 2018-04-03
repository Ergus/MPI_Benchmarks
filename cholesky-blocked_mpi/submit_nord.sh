#!/bin/bash

#BSUB -q bsc_cs
#BSUB -x
#BSUB -R "span[ptile=1]"

module purge
module load transfer/1.0
module load bsc/current
module load BASH/4.3.42
module load gcc/6.2.0
module load hwloc/1.11.8
module load impi/2017.1

#export OMP_NUM_THREADS=__NUMTHREADS
export OMP_NUM_THREADS=16

echo "# Starting: " $(date) >&2
start=$SECONDS
mpirun ./__EXEC __DIM __BS 0
end=$SECONDS
echo "# Ending: " $(date) >&2
echo "# Elapsed: "$((end-start)) >&2
