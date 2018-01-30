#!/bin/bash
#BSUB -n NODES
#BSUB -R "span[ptile=1]"
#BSUB -J EXEC_NODES
#BSUB -W 00:30
#BSUB -M 64000
#BSUB -C 10000000

module unload intel
module unload openmpi
module unload mpich
module load GCC/5.1.0 openmpi/1.6.4-mpithread BOOST/1.63.0
export NANOS6=debug
export NANOS6_SCHEDULER=distributed
echo "Starting executable: EXEC dim: DIM threads: THREADS nodes: NODES prefix: PREFIX"
mpirun ./EXEC DIM THREADS PREFIX

module load intel python/3.5.0

./prove.py PREFIX_{A,B,C}.mat
match=$?
rm -r PREFIX_{A,B,C}.mat
if [[ match -ne 0 ]]; then
    >&2 echo 'Matrices do not match'
    exit 1
fi
