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

#ifndef matmul_h
#define  matmul_h

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sched.h>
#include "benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct {
		int rank, worldsize;
		size_t dim, tthreads, cpu_count;
		size_t ldim, lthreads, first_local_thread;
		bool IprintA, IprintB, IprintC;
	} envinfo;

    extern envinfo _env;

    void Initialize(int *argc, char ***argv,
	                const size_t dim, const size_t lthreads);
    void Finalize();

	void init(double * const __restrict__ array,
	          const size_t rows, const size_t cols);

	void __print(const double * const mat,
	             const size_t rows, const size_t cols,
	             const char name[64]);

#define printmatrix(mat, rows, cols) __print(mat, rows, cols, #mat)

#ifdef __cplusplus
}
#endif

#endif
