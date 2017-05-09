#!/bin/bash
#BSUB -R "span[ptile=1]"
#BSUB -J mpi_matmul
#BSUB -W 00:30
#BSUB -M 64000
#BSUB -C 10000000

module unload intel
module unload openmpi
module unload mpich
module load GCC/5.1.0 openmpi/1.6.4-mpithread BOOST/1.63.0
export NANOS6=debug
export NANOS6_SCHEDULER=distributed
mpirun ./matmul_mpi 4096 128 print
if [[ $? -ne 0 ]]; then
    >&2 echo "mpirun for $i nodes failed"
    exit 1
fi
