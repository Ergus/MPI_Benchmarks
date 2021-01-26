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

#ifndef BENCHMARKS_MPI
#define BENCHMARKS_MPI

#ifdef __cplusplus
extern "C" {
#endif

#include "cmacros/macros.h"
#include "mpi.h"

	typedef struct {
		int rank, worldsize;
		size_t maxthreads, cpu_count;
		size_t dim, ldim, first_local_thread;
		int printerA, printerB, printerC;
	} envinfo;

	void Initialize(envinfo * _env, int *argc, char ***argv, size_t dim, size_t TS)
	{
		MPI_Init(argc, argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &(_env->rank));
		MPI_Comm_size(MPI_COMM_WORLD, &(_env->worldsize));

		omp_set_dynamic(0);	// Disable dynamic teams
		omp_set_schedule(omp_sched_static, TS);
		_env->maxthreads = omp_get_max_threads();

		// test that the cpuset >= number of local threads
		_env->cpu_count = count_sched_cpus();
		myassert(_env->cpu_count >= 0);

		myassert(dim >= (size_t)_env->worldsize);	// more rows than task size
		modcheck(dim, _env->worldsize);	// we need to split exactly

		_env->dim = dim;
		_env->ldim = dim / _env->worldsize;	// rows for local process
		_env->first_local_thread = _env->rank * _env->maxthreads; // first local thread id

		// Helper to print.
		_env->printerA = 0;	// prints A
		_env->printerB = imin(1, _env->worldsize - 1);	// prints B
		_env->printerC = imin(2, _env->worldsize - 1);	// prints C
	}

	void Finalize() {
		MPI_Finalize();
	}

#ifdef __cplusplus
}
#endif

#endif				// BENCHMARKS_MPI