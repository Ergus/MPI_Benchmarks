
#ifndef matmul_h
#define  matmul_h

#include "benchmarks.h"

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

    #define printmatrix(mat, dim, pref) {__print(mat, dim, #mat, pref); printf("printed %s\n", #mat); }

#ifdef __cplusplus
}
#endif

#endif
