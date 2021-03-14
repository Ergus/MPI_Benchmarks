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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>

#include "ArgParserC/argparser.h"

#include "benchmarks_mpi.h"

#if ISMATVEC
#define PREFIX "matvec"
#else
#define PREFIX "matmul"
#endif

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

void jacobi_omp(const double *A, const double *B,
                double *xin, double *xout, const envinfo *env
) {
	const size_t first_row = env->ldim * env->rank;

	#pragma omp parallel for
	for (size_t i = 0; i < env->ldim; ++i) {
		size_t fi = first_row + i;
		xout[i] = B[fi];

		for (size_t j = 0; j < env->dim; ++j) {
			xout[i] += A[i * env->dim + j] * xin[j];
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

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const int ROWS = create_cl_int("Rows");
	const int TS = create_cl_int("Tasksize");
	const int ITS = create_optional_cl_int("Iterations", 1);
	const int PRINT = create_optional_cl_int("Print", 0);

	envinfo env;

	Initialize(&env, &argc, &argv, ROWS, TS);

	printf("# Initializing data in process: %d\n", env.rank);;

	timer ttimer = create_timer("Total_time");

	// Allocate memory
	const size_t rowsA = (env.printerA == env.rank ? env.dim : env.ldim);
	double *A = (double *) malloc(rowsA * env.dim * sizeof(double));
	double *const lA = (env.printerA == env.rank
	                    ? &A[env.rank * env.ldim * env.dim] : A);

	double *B = (double *) malloc(env.dim * sizeof(double));
	double *x1 = (double *) malloc(env.dim * sizeof(double));
	double *x2 = (double *) malloc(env.ldim * sizeof(double));

	// Initialize arrays local portions
	init_AB(lA, B, &env);
	jacobi_modify(lA, B, &env);
	init_x(x1, env.dim, 0);
	init_x(x2, env.ldim, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");

	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, env.ldim, MPI_DOUBLE,
	              B, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {
		jacobi_omp(lA, B, x1, x2, &env);

		MPI_Allgather(x2, env.ldim, MPI_DOUBLE,
		              x1, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);
	}

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm in process: %d\n", env.rank);
	stop_timer(&ttimer);

	if (env.rank == 0) {
		create_reportable_int("worldsize", env.worldsize);
		create_reportable_int("cpu_count", env.cpu_count);
		create_reportable_int("maxthreads", env.maxthreads);
		report_args();
	}

	if (PRINT) {
		// Gather C to ITS printer
		printf("# Call C Gather in process: %d\n", env.rank);
		if (env.printerC == env.rank) {
			printmatrix(x1, env.dim, 1, PREFIX);
		}

		// Gather A to its printer
		printf("# Call A Gather in process: %d\n", env.rank);
		if (env.printerA == env.rank) {
			MPI_Gather(MPI_IN_PLACE, env.ldim * env.dim, MPI_DOUBLE,
			           A, env.ldim * env.dim, MPI_DOUBLE,
			           env.printerA, MPI_COMM_WORLD);

			printf("# Print A in process: %d\n",env.rank);
			printmatrix(A, env.dim, env.dim, PREFIX);
		} else {
			MPI_Gather(A, env.ldim * env.dim, MPI_DOUBLE,
			           NULL, env.ldim * env.dim, MPI_DOUBLE,
			           env.printerA, MPI_COMM_WORLD);
		}

		if (env.printerB == env.rank) {
			printf("# Print B in process: %d\n", env.rank);
			printmatrix(B, env.dim, 1, PREFIX);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (env.rank) {
			printf("# Done printing results...\n");
		}
	}

	free(x2);
	free(x1);
	free(B);
	free(A);
	free_args();

	Finalize();
	return 0;
}
