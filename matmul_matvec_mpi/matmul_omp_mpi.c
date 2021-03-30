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
#include <libgen.h>

#include "ArgParserC/argparser.h"

#include "benchmarks_mpi.h"

void matmul_omp(const double *A, const double *B, double * const C,
                size_t lrowsA, size_t dim, size_t colsBC)
{
	#pragma omp parallel for
	for (size_t i = 0; i < lrowsA; ++i) {
		for (size_t k = 0; k < colsBC; ++k)
			C[i * colsBC + k] = 0.0;

		for (size_t j = 0; j < dim; ++j) {
			const double temp = A[i * dim + j];

			for (size_t k = 0; k < colsBC; ++k) {
				C[i * colsBC + k] += (temp * B[j * colsBC + k]);
			}
		}
	}
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
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
	double *A = (double *)malloc(rowsA * env.dim * sizeof(double));
	double *const lA = (env.printerA == env.rank
	                    ? &A[env.rank * env.ldim * env.dim] : A);

	const size_t colsBC = (ISMATVEC ? 1 : env.dim);
	const int nelsBC = env.ldim * colsBC;
	double *B = (double *)malloc(env.dim * colsBC * sizeof(double));
	double *const lB = &B[env.rank * nelsBC];	// B is fully distributed

	const size_t rowsC = (env.printerC == env.rank ? env.dim : env.ldim);
	double *C = (double *)malloc(rowsC * colsBC * sizeof(double));
	double *const lC = (env.printerC == env.rank
	                    ? &C[env.rank * nelsBC] : C);

	// Initialise arrays local portions
	matrix_init(lA, env.ldim, env.dim, env.first_local_thread);
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
		matmul_omp(lA, B, lC, env.ldim, env.dim, colsBC);
	}

	MPI_Barrier(MPI_COMM_WORLD);

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
			MPI_Gather(MPI_IN_PLACE, nelsBC, MPI_DOUBLE,
			           C, nelsBC, MPI_DOUBLE,
			           env.printerC, MPI_COMM_WORLD);

			printf("# Print C in process %d\n", env.rank);
			printmatrix(C, env.dim, colsBC, PREFIX);
		} else {
			MPI_Gather(C, nelsBC, MPI_DOUBLE,
			           NULL, nelsBC, MPI_DOUBLE,
			           env.printerC, MPI_COMM_WORLD);
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
			printmatrix(B, env.dim, colsBC, PREFIX);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (env.rank) {
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
