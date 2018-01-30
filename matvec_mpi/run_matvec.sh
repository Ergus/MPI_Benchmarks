#!/bin/bash

#BSUB -R "span[ptile=1]"
#BSUB -cwd @CMAKE_BINARY_DIR@

module unload openmpi intel
module load gcc/5.1.0 openmpi/3.0.0

ulimit -c unlimited
ulimit -s 10240

if ( WITHEXTRAE ); then
    module load extrae_thread/3.5.2

    export EXTRAE_CONFIG_FILE="EXTRAEFILE"
    export LD_PRELOAD=${EXTRAE_HOME}/lib/libompitrace.so
fi

mpirun -np NODES ./matvec_mpi/matvec-allgather SIZE THREADS ITS
