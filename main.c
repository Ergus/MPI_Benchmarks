#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define myassert(x){                                            \
        if(!(x)){                                               \
            fprintf(stderr,"Error: %s returned %d\n",#x, x);    \
            printme();                                          \
            MPI_Abort(MPI_COMM_WORLD, -1);                      \
            }                                                   \
        }

int main(){

    int rank, size;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    const int arraysize=5*size;
    int *array=(int*) malloc(arraysize*sizeof(int));
    printf("Process %d: arraysize: %d\n",rank,arraysize);
    
    for(int i=0, cont=rank*5;i<5;++i, ++cont){
        array[cont]=cont;
        }

    MPI_Allgather(MPI_IN_PLACE, 5, MPI_INT,
                  array, 5, MPI_INT, MPI_COMM_WORLD);

    
    
//    const int arraysize=(rank==0 ? 5*size : 5);
//    
//    int *array=(int*) malloc(arraysize*sizeof(int));
//
//    printf("Process %d: arraysize: %d\n",rank,arraysize);
//    
//
//    for(int i=0, cont=rank*5;i<5;++i, ++cont){
//        array[i]=cont;
//        }
 //    
//    int displs[size], recvcounts[size];
//    for(int i=0;i<size;i++){
//        displs[i]=i*5;
//        recvcounts[i]=5;
//        }    
//
//    MPI_Gatherv(rarray, 5, MPI_INT,
//                array, recvcounts, displs, MPI_INT,
//                0, MPI_COMM_WORLD);

//    if(rank==0)
//        MPI_Gather(MPI_IN_PLACE, 5, MPI_INT,
//                   array, 5, MPI_INT,
//                   0, MPI_COMM_WORLD);
//    else
//        MPI_Gather(array, 5, MPI_INT,
//                   NULL, 5, MPI_INT,
//                   0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
        for(int i=0;i<size*5;++i){
            printf("array[%d] = %d\n",i,array[i]);
            }
        }
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==1){
        for(int i=0;i<size*5;++i){
            printf("array[%d] = %d\n",i,array[i]);
            }
        }
    fflush(stdout);

    MPI_Finalize();
    }
