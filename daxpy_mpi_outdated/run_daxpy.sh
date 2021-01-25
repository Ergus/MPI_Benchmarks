#!/bin/bash

#BSUB -R "span[ptile=1]"
#BSUB -cwd @CMAKE_BINARY_DIR@

module unload intel openmpi
module load gcc/5.1.0 openmpi/2.1.2

if ( WITHEXTRAE ); then
    module load extrae_thread/3.5.1

    export EXTRAE_CONFIG_FILE="EXTRAEFILE"
    export LD_PRELOAD=${EXTRAE_HOME}/lib/libompitrace.so
fi

ulimit -c unlimited
ulimit -s 10240
mpirun -np NR_PROCS ./daxpy_mpi/daxpy-gather-omp SIZE TS ALPHA
