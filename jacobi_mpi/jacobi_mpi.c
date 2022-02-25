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

void init_x(double *x, const size_t dim, const double val)
{
	#pragma omp parallel for
	for (size_t i = 0; i < dim; ++i) { // loop nodes
		x[i] = val;
	}
}

#if TASKTYPE == 0 || TASKTYPE == 1  // Serial and Parallel for.
void init_AB(const envinfo *env, double A[env->ldim][env->dim], double B[env->ldim])
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#if TASKTYPE == 1
	#pragma omp parallel for
	#endif // TASKTYPE == 1
	for (size_t i = 0; i < env->ldim; ++i) {
		init_AB_i(dim, A[i], &B[i], first_row + i);
	}
}

void jacobi_modify(const envinfo *env, double A[env->ldim][env->dim], double B[env->ldim])
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#if TASKTYPE==1
	#pragma omp parallel for
	#endif // TASKTYPE == 1
	for (size_t i = 0; i < env->ldim; ++i) {
		jacobi_modify_i(dim, A[i], &B[i], first_row + i);
	}
}

// A * xin + B = xout
void jacobi_mpi(const envinfo *env,
                double A[env->ldim][env->dim], const double B[env->ldim],
                const double xin[env->dim], double xout[env->ldim],
                size_t it
) {
	if (it == 0) {
		printf("# jacobi with parallelfor (node: %d)\n", env->rank);
	}

	const size_t first_row = env->ldim * env->rank;

	#if TASKTYPE == 1
	#pragma omp parallel for
	#endif // TASKTYPE == 1
	for (size_t i = 0; i < env->ldim; ++i) {
		jacobi(env->dim, 1, &A[i], &B[i], xin, &xout[i]);
	}
}

#elif TASKTYPE == 2  // With tasks

void init_AB(const envinfo *env, double A[env->ldim][env->dim], double B[env->ldim])
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += env->ts) { // loop tasks
			#pragma omp task depend(out:A[i]) depend(out:B[i])
			{
				for (size_t j = i; j < i + env->ts; ++j) {
					init_AB_i(dim, A[j], &B[j], first_row + j);
				}
			}
		}
	}
}

void jacobi_modify(const envinfo *env, double A[env->ldim][env->dim], double B[env->ldim])
{
	const size_t first_row = env->ldim * env->rank;
	const size_t dim = env->dim;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += env->ts) {
			#pragma omp task depend(inout:A[i]) depend(inout:B[i])
			{
				for (size_t j = i; j < i + env->ts; ++j) {
					jacobi_modify_i(dim, A[j], &B[j], first_row + j);
				}
			}
		}
	}
}

// A * xin + B = xout
void jacobi_mpi(const envinfo *env,
                double A[env->ldim][env->dim], const double B[env->ldim],
                const double xin[env->dim], double xout[env->ldim],
                size_t it
) {
	if (it == 0) {
		printf("# jacobi with tasks (node: %d)\n", env->rank);
	}

	const size_t dim = env->dim;
	const size_t ts = env->ts;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < env->ldim; i += ts) {
			#pragma omp task depend(in:A[i * dim])			   \
				depend(in:xin[0])							   \
				depend(in:B[i])								   \
				depend(out:xout[i])
			{
				jacobi(dim, ts, &A[i], &B[i], xin, &xout[i]);
			}
		}
	}
}
#else // TASKTYPE
#error "Invalid TASKTYPE type."
#endif // TASKTYPE


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
	double (*A)[ROWS] = malloc(env.ldim * env.dim * sizeof(double));
	double *B = (double *) malloc(env.ldim * sizeof(double));
	double *lx = (double *) malloc(env.ldim * sizeof(double));
	double *gx = (double *) malloc(env.dim * sizeof(double));  // Init in a gather

	// Initialize arrays local portions
	init_AB(&env, A, B);
	jacobi_modify(&env, A, B);
	init_x(lx, env.ldim, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {

#if TASKTYPE == 0 // No gather required.
#elif P2P == 0
		dbprintf("# Call MPI_Allgather\n");
		MPI_Allgather(lx, env.ldim, MPI_DOUBLE,
		              gx, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);
#elif P2P == 1
		jacobi_Allgather_p2p(lx, env.ldim, MPI_DOUBLE,
		                     gx, env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);
#else
		#error Invalid p2p value
#endif

		jacobi_mpi(&env, A, B, gx, lx, i);
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
		printmatrix_mpi((double *) A, env.ldim, env.dim, PREFIX, &env);
		printmatrix_mpi((double *) B, env.ldim, 1, PREFIX, &env);

		double *x = lx;
		printmatrix_mpi((double *) x, env.ldim, 1, PREFIX, &env);

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
