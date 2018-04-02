#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_MKL_H
#   include <mkl.h>
#else
#   include <cblas.h>
#   include <lapacke.h>
#endif

void check_factorization(const size_t n, double A1[n][n], double A2[n][n],
	                         const size_t lda, const double eps);

void fill_block_ij(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i, const size_t j);

void fill_block_ii(const size_t nblocks, const size_t bsize,
                        double matrix[nblocks][nblocks][bsize][bsize],
                        const size_t i);

void initialize_matrix_blocked(const size_t nblocks, const size_t bsize,
                               double matrix[nblocks][nblocks][bsize][bsize]);

void blocked2flat(size_t nblocks, size_t bsize, size_t ld,
                       double blocked[nblocks][nblocks][bsize][bsize],
                       double flat[nblocks*bsize][nblocks*bsize]);

void write_matrix(const char filename[64], size_t ld, double matrix[ld][ld]);

void write_matrix_blocked(const char filename[64],
	                          const size_t nblocks, const size_t bsize,
	                          double matrix[nblocks][nblocks][bsize][bsize]);

// Blas wrappers

enum blas_cmach_type {
	blas_base      = 151,
	blas_t         = 152,
	blas_rnd       = 153,
	blas_ieee      = 154,
	blas_emin      = 155,
	blas_emax      = 156,
	blas_eps       = 157,
	blas_prec      = 158,
	blas_underflow = 159,
	blas_overflow  = 160,
	blas_sfmin     = 161
};

static inline double BLAS_dpow_di(double x, int n)
{
	double rv = 1.0;

	if (n < 0) {
		n = -n;
		x = 1.0 / x;
	}

	for (; n; n >>= 1, x *= x) {
		if (n & 1) rv *= x;
	}

	return rv;
}

static inline double BLAS_dfpinfo(enum blas_cmach_type cmach)
{
	const double b = 2.0;
	const int    t = 53, m = -1021; // l = 1024, 
	const double eps = BLAS_dpow_di(b, -t);
	const double r   = BLAS_dpow_di(b, m-1);

	//double o = ((1.0 - eps) * BLAS_dpow_di(b, l-1)) * b;

	switch (cmach) {
	case blas_eps:   return eps;
	case blas_sfmin: return r;
	default:
		fprintf(stderr, "%s %d %d %d\n", __func__, -1, cmach, 0);
		abort();
		break;
	}

	return 0.0;
}

void oss_potrf(const size_t bsize, double A[bsize][bsize]);
void oss_trsm(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize]);
void oss_syrk(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize]);
void oss_gemm(const size_t bsize, double A[bsize][bsize],
	              double B[bsize][bsize], double C[bsize][bsize]);

#endif // UTIL_H
