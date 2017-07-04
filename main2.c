#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define printme() {                                             \
        fprintf(stderr,"%s in %s:%d (process %s)\n",            \
                __PRETTY_FUNCTION__, __FILE__, __LINE__,        \
                getenv("OMPI_COMM_WORLD_RANK"));                \
        }

#define modcheck(a, b){                                         \
        const int themod=(a) % (b);                             \
        if(themod) {                                            \
            fprintf(stderr,"Error: %s %% %s = %d\n", #a, #b, themod);   \
            printme();                                          \
            MPI_Abort(MPI_COMM_WORLD, -1);                      \
            }                                                   \
        }

#define myassert(x){                                            \
        if(!(x)){                                               \
            fprintf(stderr,"Error: %s returned %d\n",#x, x);    \
            printme();                                          \
            MPI_Abort(MPI_COMM_WORLD, -1);                      \
            }                                                   \
        }


int main(int argc, char** argv){

    printf("argv[1]=%s argv[2]=%s\n",argv[1],argv[2]);
    
    MPI_Init(&argc, &argv);

    int a=8, b=2;
    
    modcheck(a, b);
    myassert(a<b);
    
    MPI_Finalize();
    return 0;
    }
