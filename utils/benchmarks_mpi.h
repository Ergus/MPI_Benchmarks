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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <libgen.h>  // basename
#include <assert.h>
#include "mpi.h"

#include "cmacros/macros.h"
#include "ArgParserC/argparser.h"

// Define extrae if macro is set.
#if __WITH_EXTRAE
#include "extrae_user_events.h"

typedef extrae_type_t inst_type_t;
typedef extrae_value_t inst_value_t;
#define inst_define_event_type(type,name,nvalues,values,descriptions) \
	Extrae_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val) Extrae_event(evt, val)
#else // __WITH_EXTRAE
typedef size_t inst_type_t;
typedef size_t inst_value_t;
#define inst_define_event_type(type,name,nvalues,values,descriptions)
#define inst_event(evt, val)
#endif // __WITH_EXTRAE

	typedef struct {
		int rank, worldsize;
		size_t maxthreads, cpu_count, TS;
		size_t dim, ldim, first_local_thread;
		int printerA, printerB, printerC;
	} envinfo;

	void Initialize(envinfo * _env, int *argc, char ***argv, size_t dim, size_t TS)
	{
		int provided = 0;
		MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
		assert (provided == MPI_THREAD_MULTIPLE);

		MPI_Comm_rank(MPI_COMM_WORLD, &(_env->rank));
		MPI_Comm_size(MPI_COMM_WORLD, &(_env->worldsize));

		omp_set_dynamic(0);	// Disable dynamic teams
		omp_set_schedule(omp_sched_static, TS);
		_env->maxthreads = omp_get_max_threads();

		// test that the cpuset >= number of local threads
		_env->cpu_count = count_sched_cpus();
		myassert(_env->cpu_count >= 0);

		_env->TS = TS;
		myassert(_env->TS > 0);

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
