
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

    void Initialize(int *argc, char ***argv, size_t dim, size_t TS);

    void init(double* __restrict__ array,
			  size_t rows, size_t cols);

    void matvec(double * __restrict__ A,
				double * __restrict__ B,
				double * __restrict__ C,
				int rowsA, int colsA);

    void __print(double* mat, size_t rows, size_t cols,
				 const char *name, const char* prefix);

    void Finalize();

#define printmatrix(mat, rows, cols, pref) __print(mat, rows, cols,#mat, pref)

#ifdef __cplusplus
}
#endif

#endif
