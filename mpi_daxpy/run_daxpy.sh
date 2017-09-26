#!/bin/bash

#BSUB -R "span[ptile=1]"
#BSUB -cwd ${CMAKE_BINARY_DIR}

module unload openmpi
module load openmpi/1.6.4-mpithread
module load gcc/5.1.0

ulimit -c unlimited
mpirun -np NR_PROCS ./mpi_daxpy/daxpy-gather-omp SIZE TS ALPHA PRINT
