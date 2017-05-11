
#ifndef daxpy_h
#define daxpy_h

#include "benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    typedef struct{
        int rank, worldsize;
        size_t dim, tthreads;
        size_t ldim, lthreads, first_local_thread;
        bool IprintX; // Y always printed in rank 0, but X depends
        } envinfo;

    extern envinfo _env;

    void Initialize(int *argc,char ***argv,size_t dim,size_t N);

    void init(double* array, size_t ldim);

    void daxpy(double *lY,double a,double *X,size_t ldim);  // remember this is X+=a*Y

    void __print(double* mat, size_t dim, const char *name, const char* prefix);

    void Finalize();

    #define printmatrix(mat, dim, pref) __print(mat, dim, #mat, pref); 

#ifdef __cplusplus
}
#endif

#endif
