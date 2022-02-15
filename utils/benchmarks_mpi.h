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

#ifndef BENCHMARKS_MPI_H
#define BENCHMARKS_MPI_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>  // basename
#include <assert.h>
#include <unistd.h>
#include <errno.h>

#include <mpi.h>
#include <omp.h>

#include "cmacros/macros.h"
#include "ArgParserC/argparser.h"


// Declare some blas routines.
#include <mkl.h>
#include <limits.h>

void dcopy_(const int *n, const double *dx, const int *incx, double *dy, const int *incy);

void dgemv_ (const char *trans, const int *m, const int *n,
             const double *alpha, const double *A, const int *lda,
             const double *x, const int *incx,
             const double *beta, double *y, const int *incy);

void dgemm_(const char *transa, const char *transb,
            const int *l, const int *n, const int *m,
            const double *alpha, const void *a, const int *lda,
            const void *b, const int *ldb,
            const double *beta, void *c, const int *ldc);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double *a, int *lda, double *b, int *ldb);

void dsyrk_(char *uplo, char *trans, int *n, int *k,
            double *alpha, double *a, int *lda,
            double *beta, double *c, int *ldc);

#if __WITH_EXTRAE // #####################

#include <extrae.h>
#include "extrae_user_events.h"

typedef extrae_type_t inst_type_t;
typedef extrae_value_t inst_value_t;

#define inst_define_event_type(type,name,nvalues,values,descriptions) \
	Extrae_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val) Extrae_event(evt, val)

#define USER_EVENT 9910002
#define USER_EVENT_VALUES						\
	EVENT(USER_NONE)							\
	EVENT(USER_MATVEC)							\
	EVENT(USER_MATMUL)							\
	EVENT(USER_JACOBI)

	enum user_values_t {
		#define EVENT(evt) evt,
		BLAS_EVENT_VALUES
		#undef EVENT
		USER_NEVENTS
	};

#define BLAS_EVENT 9910003
#define BLAS_EVENT_VALUES						\
	EVENT(BLAS_NONE)							\
	EVENT(BLAS_POTRF)							\
	EVENT(BLAS_TRSM)							\
	EVENT(BLAS_GEMM)							\
	EVENT(BLAS_GEMV)							\
	EVENT(BLAS_COPY)							\
	EVENT(BLAS_SYRK)

	enum blas_values_t {
		#define EVENT(evt) evt,
		BLAS_EVENT_VALUES
		#undef EVENT
		BLAS_NEVENTS
	};

	void inst_register_events()
	{
		extrae_type_t user_event = USER_EVENT;
		extrae_type_t blas_event = BLAS_EVENT;

		unsigned user_nvalues = USER_NEVENTS;
		unsigned blas_nvalues = BLAS_NEVENTS;

		#define EVENT(evt) (extrae_value_t) evt,
		static extrae_value_t blas_values[BLAS_NEVENTS] = {
			BLAS_EVENT_VALUES
		};

		static extrae_value_t user_values[USER_NEVENTS] = {
			USER_EVENT_VALUES
		};
		#undef EVENT

		#define EVENT(evt) #evt,
		static char *blas_names[BLAS_NEVENTS] = {
			BLAS_EVENT_VALUES
		};

		static char *blas_names[BLAS_NEVENTS] = {
			USER_EVENT_VALUES
		};
		#undef EVENT

		inst_define_event_type(&user_event, "user_event", &user_nvalues, user_values, user_names);
		inst_define_event_type(&blas_event, "blas_event", &blas_nvalues, blas_values, blas_names);
	}

#else // __WITH_EXTRAE // #####################

	typedef size_t inst_type_t;
	typedef size_t inst_value_t;

#define inst_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val)

#define BLAS_EVENT 0
#define inst_register_events()

#endif // __WITH_EXTRAE // #####################

typedef struct {
	int rank, worldsize;
	size_t maxthreads, cpu_count, ts;
	size_t dim, ldim, first_local_thread;
	int printerA, printerB, printerC;
} envinfo;

void Initialize(envinfo * _env, int *argc, char ***argv, size_t dim, size_t ts)
{
	int provided = 0;
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
	assert (provided == MPI_THREAD_MULTIPLE);

	MPI_Comm_rank(MPI_COMM_WORLD, &(_env->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(_env->worldsize));

	omp_set_dynamic(0);	// Disable dynamic teams
	omp_set_schedule(omp_sched_dynamic, ts);
	_env->maxthreads = omp_get_max_threads();

	// test that the cpuset >= number of local threads
	_env->cpu_count = count_sched_cpus();
	myassert(_env->cpu_count >= 0);

	// It seems like OMP gets the maxthreads from the total number of CPU
	// and ignored the count_sched_cpus. So for this benchmarks I prefer to
	// set them manually and emit a warning. I also tried to manage this
	// using the external variable OMP_NUM_THREADS, but something is missing
	// and the variable does not take effect sometimes.
	if (_env->maxthreads > _env->cpu_count) {
		dbprintf("# Setting omp_set_num_threads because maxthreads=%zu and cpu_count=%zu\n",
		         _env->maxthreads, _env->cpu_count);
		omp_set_num_threads((int) _env->cpu_count);
		// To assert that it worked.
		_env->maxthreads = omp_get_max_threads();
	}

	_env->ts = ts;
	myassert(_env->ts > 0);
	modcheck(dim, _env->ts)

		myassert(dim >= (size_t)_env->worldsize);	// more rows than task size
	modcheck(dim, _env->worldsize);	// we need to split exactly

	_env->dim = dim;
	_env->ldim = dim / _env->worldsize;	// rows for local process
	_env->first_local_thread = _env->rank * _env->maxthreads; // first local thread id

	// Helper to print.
	_env->printerA = 0;	// prints A
	_env->printerB = imin(1, _env->worldsize - 1);	// prints B
	_env->printerC = imin(2, _env->worldsize - 1);	// prints C

	inst_register_events();
}

void Finalize() {
	MPI_Finalize();
}


void __printslice(const double * const mat,
                  const size_t rows, const size_t cols,
                  const char prefix[64], const char name[64],
                  const envinfo * const env
) {
	for (int i = 0; i < env->worldsize; ++i) {
		if (i == env->rank) {
			FILE *fp = NULL;

			if (env->rank == 0) {
				printf("# Printing %s\n", name);
				fp = get_file(prefix, name, "w+");
				print_matrix_header(fp, name, env->dim, cols);
			} else {
				fp = get_file(prefix, name, "a");
			}
			print_matrix_data(fp, mat, rows, cols);

			if (fp != stdout)
				fclose(fp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

#define printmatrix_mpi(mat, rows, cols, prefix, env)	\
	__printslice(mat, rows, cols, prefix, #mat, env)

#ifdef __cplusplus
}
#endif

#endif // BENCHMARKS_MPI_H
