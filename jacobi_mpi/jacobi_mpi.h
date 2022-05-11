/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
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

#ifndef JACOBI_MPI_H
#define JACOBI_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "benchmarks_mpi.h"

	void init_AB_i(size_t dim, double Ai[dim], double Bi[1], size_t global_i)
	{
		struct drand48_data drand_buf;
		srand48_r(global_i, &drand_buf);     // Init rng per row.
		double cum = 0.0, sum = 0.0, x;

		for (size_t j = 0; j < dim; ++j) {
			drand48_r(&drand_buf, &x);
			Ai[j] = x;
			cum += fabs(x);
			sum += x;
		}

		const double valii = Ai[global_i];
		if (signbit(valii)) {
			Ai[global_i] = valii - cum;
			Bi[0] = sum - cum;
		} else {
			Ai[global_i] = valii + cum;
			Bi[0] = sum + cum;
		}

	}

	void jacobi_modify_i(size_t dim, double Ai[dim], double Bi[1], size_t global_i)
	{
		const double iAii = 1 / fabs(Ai[global_i]);

		for (size_t j = 0; j < dim; ++j) {
			Ai[j] = (global_i == j) ? 0.0 : -1.0 * Ai[j] * iAii;
		}
		Bi[0] *= iAii;
	}

	int jacobi_Allgather_p2p(const void* b_send, int nsend, MPI_Datatype type_send,
	                         void* b_recv, int nrecv, MPI_Datatype type_recv,
	                         MPI_Comm comm
	) {
		dbprintf("# Call %s\n", __func__);
		int worldsize = -1;
		MPI_Comm_size(comm, &worldsize);
		if (worldsize == 1)
			return MPI_SUCCESS;

		myassert(b_send != MPI_IN_PLACE); // Simplest supported.

		int rank = -1;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		int size_type_send = 0, size_type_recv = 0;
		MPI_Type_size(type_send, &size_type_send);
		MPI_Type_size(type_recv, &size_type_recv);

		myassert(nsend * size_type_send == nrecv * size_type_recv); // This is required.

		const int nrequests = 2 * (worldsize - 1);
		MPI_Request *reqs =
			(MPI_Request *) malloc(nrequests * sizeof(MPI_Request));

		int it = 0;

		for (int i = 0; i < worldsize; ++i) {
			if (i == rank) {
				continue;
			}
			MPI_Isend(b_send, nsend, type_send, i, rank, comm, &reqs[it++]);
		}

		for (int i = 0; i < worldsize; ++i) {
			void *b_recv_i = &((char*)b_recv)[i * nrecv * size_type_recv];

			if (i == rank) {
				memcpy(b_recv_i, b_send, nsend * size_type_send);
			} else {
				MPI_Irecv(b_recv_i, nrecv, type_recv, i, i, comm, &reqs[it++]);
			}
		}

		myassert(it == nrequests);

		MPI_Waitall(nrequests, reqs, MPI_STATUSES_IGNORE);

		free(reqs);

		return MPI_SUCCESS;
	}

	void jacobi(size_t dim, size_t ts, double A[ts][dim], const double B[ts],
	            const double xin[dim], double xout[ts]
	) {
		#if BLAS == 0
		inst_event(USER_EVENT, USER_JACOBI);
		for (size_t i = 0; i < ts; ++i) {
			xout[i] = B[i];

			for (size_t j = 0; j < dim; ++j) {
				xout[i] += A[i][j] * xin[j];
			}
		}
		inst_event(USER_EVENT, USER_NONE);

		#elif BLAS == 1

		myassert(dim < (size_t) INT_MAX);

		const char TR = 'T';
		const int M = (int) dim;
		const int N = (int) ts;
		const double alpha = 1.0;
		const double beta = 1.0;
		const int inc = 1;

		inst_event(BLAS_EVENT, BLAS_COPY);
		dcopy_(&N, B, &inc, xout, &inc);

		inst_event(BLAS_EVENT, BLAS_GEMV);
		dgemv_(&TR, &M, &N, &alpha, (double *)A, &M, xin, &inc, &beta, xout, &inc);

		inst_event(BLAS_EVENT, BLAS_NONE);
		#else // BLAS
		#error No valid BLAS value
		#endif // BLAS
	}


#ifdef __cplusplus
}
#endif

#endif // JACOBI_MPI_H
