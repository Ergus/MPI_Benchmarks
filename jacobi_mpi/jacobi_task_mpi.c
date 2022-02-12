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

#include "jacobi_mpi.h"


void init_AB_task(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += env->ts) { // loop tasks

			#pragma omp task depend(out:A[i * dim]) depend(out:B[first_row + i])
			for (size_t j = i; j < i + env->ts; ++j) {
				const size_t grow = first_row + j;

				struct drand48_data drand_buf;
				srand48_r(grow, &drand_buf);
				double cum = 0.0, sum = 0.0, x;

				for (size_t k = 0; k < dim; ++k) {
					drand48_r(&drand_buf, &x);
					A[j * dim + k] = x;
					cum += fabs(x);
					sum += x;
				}
				// Diagonal element condition.
				const double valii = A[j * dim + grow];
				if (signbit(valii)) {
					A[j * dim + grow] = valii - cum;
					B[grow] = sum - cum;
				} else {
					A[j * dim + grow] = valii + cum;
					B[grow] = sum + cum;
				}
			}
		}
	}
}


void init_x_task(
	double *x,
	const size_t dim,
	const double val,
	const envinfo *env
) {
	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < dim; i += env->ts) { // loop nodes

			#pragma omp task depend(out:x[i])
			{
				for (size_t j = i; j < i + env->ts; ++j) { // loop nodes
					x[j] = val;
				}
			}
		}
	}
}


void jacobi_modify_task(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += env->ts) {

			#pragma omp task depend(inout:A[i * dim]) \
				depend(inout:B[first_row + i])
			{
				for (size_t j = i; j < i + env->ts; ++j) {
					const size_t grow = first_row + j;
					const double iAjj = 1 / fabs(A[j * dim + grow]);

					for (size_t k = 0; k < dim; ++k) {
						A[j * dim + k]
							= (grow == k) ? 0.0 : -1.0 * A[j * dim + k] * iAjj;
					}
					B[grow] *= iAjj;
				}
			}
		}
	}
}

// A * xin + B = xout
void jacobi_task_mpi(const double *A, const double *B,
                     const double *xin, double *xout,
                     const envinfo *env, size_t it
) {
	if (it == 0) {
		printf("# jacobi with tasks (node: %d)\n", env->rank);
	}

	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += env->ts) {

			#pragma omp task depend(in:A[i * dim])			   \
				depend(in:xin[0])							   \
				depend(in:B[first_row + i])					   \
				depend(out:xout[i])
			{
				jacobi(&A[i * dim], &B[first_row + i], xin, &xout[i], env->ts, env->dim);
			}
		}
	}
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const size_t ROWS = create_cl_size_t("Rows");
	const size_t TS = create_cl_size_t("Tasksize");
	const size_t ITS = create_optional_cl_size_t("Iterations", 1);
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
	init_AB_task(lA, B, &env);
	jacobi_modify_task(lA, B, &env);
	init_x_task(x1, env.dim, 0, &env);
	init_x_task(x2, env.ldim, 0, &env);

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");

	if (env.worldsize > 0) {
		// Gather B to all
		MPI_Allgather(MPI_IN_PLACE, env.ldim, MPI_DOUBLE,
		              B, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);
	}

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {
		jacobi_task_mpi(lA, B, x1, x2, &env, i);

		if (env.worldsize > 0) {
			// Gather x to all
			MPI_Allgather(x2, env.ldim, MPI_DOUBLE,
			              x1, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);
		}
	}

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm in process: %d\n", env.rank);
	stop_timer(&ttimer);

	if (env.rank == 0) {
		create_reportable_int("worldsize", env.worldsize);
		create_reportable_int("cpu_count", env.cpu_count);
		create_reportable_int("omp_max_threads", env.maxthreads);
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
