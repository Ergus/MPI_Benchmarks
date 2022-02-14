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

void init_AB_taskfor(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel for
	for (size_t i = 0; i < env->ldim; ++i) {
		const size_t grow = first_row + i;

		struct drand48_data drand_buf;
		srand48_r(grow, &drand_buf);
		double cum = 0.0, sum = 0.0, x;

		for (size_t j = 0; j < dim; ++j) {
			drand48_r(&drand_buf, &x);
			A[i * dim + j] = x;
			cum += fabs(x);
			sum += x;
		}

		const double valii = A[i * dim + grow];
		if (signbit(valii)) {
			A[i * dim + grow] = valii - cum;
			B[i] = sum - cum;
		} else {
			A[i * dim + grow] = valii + cum;
			B[i] = sum + cum;
		}
	}
}


void init_x_taskfor(double *x, const size_t dim, const double val)
{
	#pragma omp parallel for
	for (size_t i = 0; i < dim; ++i) { // loop nodes
		x[i] = val;
	}
}


void jacobi_modify_taskfor(double *A, double *B, const envinfo *env)
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel for
	for (size_t i = 0; i < env->ldim; ++i) {
		const size_t grow = first_row + i;
		const double iAii = 1 / fabs(A[i * dim + grow]);

		for (size_t k = 0; k < dim; ++k) {
			A[i * dim + k]
				= (grow == k) ? 0.0 : -1.0 * A[i * dim + k] * iAii;
		}
		B[i] *= iAii;
	}
}


// A * xin + B = xout
void jacobi_parallelfor_mpi(const double *A, const double *B,
                            const double *xin, double *xout,
                            const envinfo *env, size_t it
) {
	if (it == 0) {
		printf("# jacobi with parallel_for (node: %d)\n", env->rank);
	}

	const size_t first_row = env->ldim * env->rank;

	#pragma omp parallel for
	for (size_t i = 0; i < env->ldim; ++i) {
		jacobi(&A[i * env->dim], &B[first_row + i], xin, &xout[i], 1, env->dim);
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
	double *A = (double *) malloc(env.ldim * env.dim * sizeof(double));
	double *B = (double *) malloc(env.ldim * sizeof(double));
	double *lx = (double *) malloc(env.ldim * sizeof(double));
	double *gx = (double *) malloc(env.dim * sizeof(double)); // Init in a gather

	// Initialize arrays local portions
	init_AB_taskfor(A, B, &env);
	jacobi_modify_taskfor(A, B, &env);
	init_x_taskfor(lx, env.ldim, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {
		jacobi_Allgather(lx, env.ldim, MPI_DOUBLE,
		                 gx, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);

		jacobi_parallelfor_mpi(A, B, lx, gx, &env, i);
	}

	MPI_Barrier(MPI_COMM_WORLD);
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
		printmatrix_mpi(A, env.ldim, env.dim, PREFIX, &env);
		printmatrix_mpi(B, env.ldim, 1, PREFIX, &env);

		double *x = lx;
		printmatrix_mpi(x, env.ldim, 1, PREFIX, &env);

		MPI_Barrier(MPI_COMM_WORLD);
		if (env.rank == 0) {
			printf("# Done printing results...\n");
		}
	}

	free(lx);
	free(gx);
	free(B);
	free(A);
	free_args();

	Finalize();
	return 0;
}
