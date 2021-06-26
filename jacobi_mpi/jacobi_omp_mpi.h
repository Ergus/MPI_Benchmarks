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

	void jacobi_base(
		const double * __restrict__ A,
		double Bi,
		const double * __restrict__ xin,
		double * __restrict__ xouti, size_t dim
	);

	void jacobi(const double *A, const double *B,
	            const double *xin, double *xout, size_t ts, size_t dim
	) {
		for (size_t i = 0; i < ts; ++i) {
			inst_event(9910002, dim);

			jacobi_base(&A[i * dim], B[i], xin, &xout[i], dim);

			inst_event(9910002, 0);
		}
	}

#ifdef __cplusplus
}
#endif

#endif // MATVEC_OMP_MPI_H
