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

#include "benchmarks_mpi.h"

void matrix_init(double * const __restrict__ array,
                 const size_t rows, const size_t cols, int seed
) {
	const size_t fullsize = rows * cols;

	#pragma omp parallel
	{
		struct drand48_data drand_buf;
		srand48_r(seed + omp_get_thread_num(), &drand_buf);
		double x;

		#pragma omp for
		for(size_t i = 0; i < fullsize; ++i) {
			drand48_r(&drand_buf, &x);
			array[i] = x;
		}
	}
}


#if ISMATVEC // ==================================================

void matmul_base(const double *A, const double *B, double * const C,
                 size_t ts, size_t dim, size_t colsBC
) {
	#if BLAS == 0
	inst_event(USER_EVENT, USER_MATVEC);
	for (size_t i = 0; i < ts; ++i) {
		C[i] = 0.0;

		for (size_t j = 0; j < dim; ++j) {
			C[i] += A[i * dim + j] * B[j];
		}
	}
	inst_event(USER_EVENT, USER_NONE);

	#elif BLAS == 1

	inst_event(BLAS_EVENT, BLAS_GEMV);
	myassert(dim < (size_t) INT_MAX);
	const char TR = 'T';
	const int M = (int) dim;
	const int N = (int) ts;
	const double alpha = 1.0;
	const double beta = 0.0;
	const int incx = 1;

	dgemv_(&TR, &M, &N, &alpha, A, &M, B, &incx, &beta, C, &incx);

	inst_event(BLAS_EVENT, BLAS_NONE);
	#else // BLAS
	#error No valid BLAS value
	#endif // BLAS
}

#else // ISMATVEC ================================================

void matmul_base(const double *A, const double *B, double * const C,
                 size_t ts, size_t dim, size_t colsBC
) {
	#if BLAS == 0
	inst_event(USER_EVENT, USER_MATMUL);
	for (size_t i = 0; i < ts; ++i) {
		for (size_t k = 0; k < colsBC; ++k)
			C[i * colsBC + k] = 0.0;

		for (size_t j = 0; j < dim; ++j) {
			const double temp = A[i * dim + j];

			for (size_t k = 0; k < colsBC; ++k) {
				C[i * colsBC + k] += (temp * B[j * colsBC + k]);
			}
		}
	}
	inst_event(USER_EVENT, USER_NONE);

	#elif BLAS == 1

	inst_event(BLAS_EVENT, BLAS_GEMM);
	myassert(dim < (size_t) INT_MAX);
	const char TA = 'N';
	const char TB = 'N';
	const int M = (int) dim;
	const int N = (int) ts;
	const int K = (int) colsBC;
	const double ALPHA = 1.0;
	const int LDA = M;
	const int LDB = K;
	const double BETA = 0.0;
	const int LDC = M;

	dgemm_(&TA, &TB, &M, &N, &K, &ALPHA,
	       B, &LDB,
	       A, &LDA, &BETA,
	       C, &LDC);

	inst_event(BLAS_EVENT, BLAS_NONE);
	#else // BLAS
	#error No valid BLAS value
	#endif // BLAS
}

#endif // ISMATVEC ================================================


#if TASKTYPE == 0 // parallel for

void matmul_mpi(const double *A, const double *B, double * const C,
                const envinfo *env, size_t colsBC, size_t it
) {
	if (it == 0) {
		printf("# %s with parallel_for (node: %d)\n",
		       (ISMATVEC ? "matvec" : "matmul"), env->rank);
	}

	const size_t dim = env->dim;
	const size_t ldim = env->ldim;
	const size_t ts = env->ts;

	#pragma omp parallel for
	for (size_t i = 0; i < ldim; i += ts) {
		matmul_base(&A[i * dim], B, &C[i * colsBC], ts, dim, colsBC);
	}
}

#elif TASKTYPE == 1 // task

void matmul_mpi(const double *A, const double *B, double * const C,
                const envinfo *env, size_t colsBC, size_t it
) {
	if (it == 0) {
		printf("# %s with tasks (node: %d)\n",
		       (ISMATVEC ? "matvec" : "matmul"), env->rank);
	}

	const size_t dim = env->dim;
	const size_t ldim = env->ldim;
	const size_t ts = env->ts;

	#pragma omp parallel
	#pragma omp single
	{
		for (size_t i = 0; i < ldim; i += ts) {
			#pragma omp task
			matmul_base(&A[i * dim], B, &C[i * colsBC], ts, dim, colsBC);
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
	double *A = (double *) malloc(env.ldim * env.dim * sizeof(double));

	const size_t colsBC = (ISMATVEC ? 1 : env.dim);

	const int nelsBC = env.ldim * colsBC;
	double *B = (double *) malloc(env.dim * colsBC * sizeof(double));
	double *const lB = &B[env.rank * nelsBC];	// B is fully distributed

	double *C = (double *) malloc(env.ldim * colsBC * sizeof(double));

	// Initialize arrays local portions
	matrix_init(A, env.ldim, env.dim, env.first_local_thread);
	matrix_init(lB, env.ldim, colsBC, env.first_local_thread);
	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");

	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, nelsBC, MPI_DOUBLE,
	              B, nelsBC, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {
		matmul_mpi(A, B, C, &env, colsBC, i);
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
		printmatrix_mpi(lB, env.ldim, colsBC, PREFIX, &env);
		printmatrix_mpi(C, env.ldim, colsBC, PREFIX, &env);

		if (env.rank == 0) {
			printf("# Done printing results...\n");
		}
	}

	free(C);
	free(B);
	free(A);
	free_args();

	Finalize();
	return 0;
}
