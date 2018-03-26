
#ifndef matmul_h
#define  matmul_h

#include "benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct{
		double *arr;
		size_t len, seed;
	}input1;

	typedef struct{
		double *A, *B, *C;
		size_t rows, cols, i;
	}input2;

	typedef struct{
		int rank, worldsize;
		size_t dim, tthreads;
		size_t ldim, lthreads, first_local_thread;
		bool IprintA, IprintB, IprintC;
	} envinfo;

    extern envinfo _env;

    void Initialize(int *argc, char ***argv,
                    const size_t dim, const size_t lthreads);

	void *init_thread(void *in);

	void init(double * const array,
	          const size_t rows, const size_t cols);

	void *matvec_thread(void *in);

    void matvec(const double * const __restrict__ A,
				const double * const __restrict__ B,
				double * __restrict__ C,
				const size_t rowsA, const size_t colsA);

    void __print(const double * const mat,
                 const size_t rows, const size_t cols,
				 const char *name, const char* prefix);

    void Finalize();

#define printmatrix(mat, rows, cols, pref) __print(mat, rows, cols,#mat, pref)

#ifdef __cplusplus
}
#endif

#endif
