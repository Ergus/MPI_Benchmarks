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

#ifndef MATVEC_OMP_MPI_H
#define MATVEC_OMP_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "benchmarks_mpi.h"

void init_AB(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;

	#pragma omp parallel
	{
		struct drand48_data drand_buf;
		srand48_r(env->first_local_thread + omp_get_thread_num(),
		          &drand_buf);

		#pragma omp for
		for (size_t i = 0; i < env->ldim; ++i) {

			double cum = 0.0, sum = 0.0, x;

			for (size_t j = 0; j < env->dim; ++j) {
				drand48_r(&drand_buf, &x);
				A[i * env->dim + j] = x;
				cum += fabs(x);
				sum += x;
			}

			const double valii = A[i * env->dim + first_row + i];

			if (signbit(valii)) {
				A[i * env->dim + first_row + i] = valii - cum;
				B[first_row + i] = sum - cum;
			} else {
				A[i * env->dim + first_row + i] = valii + cum;
				B[first_row + i] = sum + cum;
			}
		}
	}
}


void init_x(double *x, const size_t dim, const double val)
{
	#pragma omp parallel for
	for (size_t i = 0; i < dim; ++i) { // loop nodes
		x[i] = val;
	}
}


void jacobi_modify(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;

	#pragma omp parallel
	{
		struct drand48_data drand_buf;
		srand48_r(env->first_local_thread + omp_get_thread_num(),
		          &drand_buf);

		#pragma omp for
		for (size_t i = 0; i < env->ldim; ++i) {

			double cum = 0.0, sum = 0.0, x;

			const double Aii = A[i * env->dim + first_row + i];

			for (size_t j = 0; j < env->dim; ++j) {
				if (first_row + i == j) {
					A[i * env->dim + j] = 0.0;
				} else {
					A[i * env->dim + j] = - (A[i * env->dim + j] / Aii);
				}
			}
			B[first_row + i] /= Aii;
		}
	}
}


#ifdef __cplusplus
}
#endif

#endif // MATVEC_OMP_MPI_H
