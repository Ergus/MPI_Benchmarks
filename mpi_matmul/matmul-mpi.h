
#ifndef matmul_h
#define  matmul_h

// C headers
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define imin(x, y) (((x) < (y)) ? (x) : (y))
#define imax(x, y) (((x) > (y)) ? (x) : (y))

#define frand()(4.*(double)rand()/(RAND_MAX)-2.) //uniform rng in [-2,2]

#define printme() {                                             \
        fprintf(stderr,"Func: %s in %s:%d (process %s)\n",      \
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


#ifdef __cplusplus
extern "C" {
#endif
    
    typedef struct{
        int rank, worldsize;
        size_t dim, tthreads;
        size_t ldim, lthreads, first_local_thread;
        bool IprintA, IprintB, IprintC;
        } envinfo;

    extern envinfo _env;

    void Initialize(int *argc,char ***argv,size_t dim,size_t N);

    void init(double* array, size_t rows, size_t cols);

    void matmul(double *A, double *B, double *C, int rows, int cols);

    void __print(double* mat, size_t dim, const char *name, const char* prefix);

    void Finalize();

#ifdef __cplusplus
}
#endif

#define printmatrix(mat, dim, pref) __print(mat, dim, #mat, pref);

#endif
