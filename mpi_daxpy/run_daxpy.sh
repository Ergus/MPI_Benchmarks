#!/bin/bash

#BSUB -R "span[ptile=1]"
#BSUB -cwd ${CMAKE_BINARY_DIR}

module unload intel openmpi
module load gcc/5.1.0
module load openmpi/1.8.1-multithread

ulimit -c unlimited
mpirun -np NR_PROCS ./mpi_daxpy/daxpy-gather-omp SIZE TS ALPHA PRINT
