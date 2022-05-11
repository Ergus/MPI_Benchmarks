/*
 * Copyright (C) 2021  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CHOLESKY_OMP_MPI_H
#define CHOLESKY_OMP_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "benchmarks_mpi.h"

static inline void omp_potrf(int ts, double A[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_POTRF);
	assert(A != NULL);

	static int INFO;
	static const char L = 'L';

	dpotrf_(&L, &ts, (double *)A, &ts, &INFO);
	inst_event(BLAS_EVENT, BLAS_NONE);
}

static inline void omp_trsm(int ts, double A[ts][ts], double B[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_TRSM);
	assert(A != NULL);
	assert(B != NULL);

    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    double DONE = 1.0;
    dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE,
	       (double *)A, &ts,
	       (double *)B, &ts);

	inst_event(BLAS_EVENT, BLAS_NONE);
}

static inline void omp_gemm(
	int ts,
	double A[ts][ts],
	double B[ts][ts],
	double C[ts][ts]
) {
	inst_event(BLAS_EVENT, BLAS_GEMM);
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);

    const char TR = 'T', NT = 'N';
    double DONE = 1.0, DMONE = -1.0;
    dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE,
	       (double *)A, &ts,
	       (double *)B, &ts, &DONE,
	       (double *)C, &ts);

	inst_event(BLAS_EVENT, BLAS_NONE);
}

static inline void omp_syrk(int ts, double A[ts][ts], double B[ts][ts])
{
	inst_event(BLAS_EVENT, BLAS_SYRK);
	assert(A != NULL);
	assert(B != NULL);

    static char LO = 'L', NT = 'N';
    static double DONE = 1.0, DMONE = -1.0;
    dsyrk_(&LO, &NT, &ts, &ts, &DMONE,
	       (double *)A, &ts, &DONE,
	       (double *)B, &ts);
	inst_event(BLAS_EVENT, BLAS_NONE);
}

#ifdef __cplusplus
}
#endif


#endif /* CHOLESKY_OMP_MPI_H */
