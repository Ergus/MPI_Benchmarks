#!/bin/bash

#BSUB -R "span[ptile=1]"
#BSUB -cwd @CMAKE_BINARY_DIR@

module unload openmpi intel

if ( $EXTRAE ); then
    module load gcc/6.2.0 openmpi/1.8.1-multithread
    module load extrae

    export EXTRAE_HOME="/apps/BSCTOOLS/extrae/3.5.1/openmpi_1_8_1_multithread"
    export EXTRAE_CONFIG_FILE="EXTRAEFILE"
    export LD_PRELOAD=${EXTRAE_HOME}/lib/libompitrace.so
else
    module load gcc/5.1.0 openmpi/1.8.1-multithread
fi

ulimit -c unlimited
mpirun -np NR_PROCS ./mpi_matvec/matvec-allgather ROWS TS ITS
