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

#include "benchmarks_mpi.h"

#include "ArgParserC/argparser.h"

#include <stdio.h>

#if ISMATVEC
#define PREFIX "matvec"
#else
#define PREFIX "matmul"
#endif

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

	const int ROWS = create_cl_int("Rows");
	const int LTHREADS = create_cl_int("lthreads");
	const int ITS = create_optional_cl_int("iterations", 1);
	const int PRINT = create_optional_cl_int("print", 0);

	envinfo _env;

	Initialize(&_env, &argc, &argv, ROWS, LTHREADS);

	printf("# Initialize in process: %d\n", _env.rank);;

	timer *ttimer = create_timer("Total time");

	// Allocate memory
	const size_t rowsA = (_env.printerA == _env.rank ? _env.dim : _env.ldim);
	double *A = (double *)malloc(rowsA * _env.dim * sizeof(double));
	double *const lA = (_env.printerA == _env.rank ?
	                    &A[_env.rank * _env.ldim * _env.dim] : A);

	const size_t colsBC = (ISMATVEC ? 1 : _env.dim);
	const int nelsBC = _env.ldim * colsBC;
	double *B = (double *)malloc(_env.dim * colsBC * sizeof(double));
	double *const lB = &B[_env.rank * nelsBC];	// B is fully distributed

	const size_t rowsC = (_env.printerC == _env.rank ? _env.dim : _env.ldim);
	double *C = (double *)malloc(rowsC * colsBC * sizeof(double));
	double *const lC = (_env.printerC == _env.rank
	                    ? &C[_env.rank * nelsBC] : C);

	// Initialise arrays local portions
	matrix_init(lA, _env.ldim, _env.dim, _env.first_local_thread);
	matrix_init(lB, _env.ldim, colsBC, _env.first_local_thread);

	//==========================================================================

	timer *atimer = create_timer("Algorithm time");

	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, nelsBC, MPI_DOUBLE,
	              B, nelsBC, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
	for (size_t i = 0; i < ITS; ++i) {
		matmul_omp(lA, B, lC, _env.ldim, _env.dim, colsBC);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	stop_timer(atimer);
	stop_timer(ttimer);
	// =========================================================================

	printf("# Finished algorithm...\n");

	if (_env.rank == 0) {
		const double performance =
			ITS * _env.dim * _env.dim * colsBC * 2000.0 / getNS_timer(atimer);

		create_reportable_double("performance", performance);
		report_args();
	}


	if (PRINT) {
		// Gather C to ITS printer
		printf("# Call C Gather in process: %d\n", _env.rank);
		if (_env.printerC == _env.rank) {
			MPI_Gather(MPI_IN_PLACE, nelsBC, MPI_DOUBLE,
			           C, nelsBC, MPI_DOUBLE,
			           _env.printerC, MPI_COMM_WORLD);

			printf("# Print C in process %d\n", _env.rank);
			printmatrix(C, _env.dim, colsBC, PREFIX);
		} else {
			MPI_Gather(C, nelsBC, MPI_DOUBLE,
			           NULL, nelsBC, MPI_DOUBLE,
			           _env.printerC, MPI_COMM_WORLD);
		}

		// Gather A to its printer
		printf("# Call A Gather in process: %d\n", _env.rank);
		if (_env.printerA == _env.rank) {
			MPI_Gather(MPI_IN_PLACE, _env.ldim * _env.dim, MPI_DOUBLE,
			           A, _env.ldim * _env.dim, MPI_DOUBLE,
			           _env.printerA, MPI_COMM_WORLD);

			printf("# Print A in process: %d\n",_env.rank);
			printmatrix(A, _env.dim, _env.dim, PREFIX);
		} else {
			MPI_Gather(A, _env.ldim * _env.dim, MPI_DOUBLE,
			           NULL, _env.ldim * _env.dim, MPI_DOUBLE,
			           _env.printerA, MPI_COMM_WORLD);
		}

		if (_env.printerB == _env.rank) {
			printf("# Print B in process: %d\n", _env.rank);
			printmatrix(B, _env.dim, colsBC, PREFIX);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (_env.rank) {
			printf("# Done printing results...\n");
		}
	}

	free(C);
	free(B);
	free(A);

	Finalize();
	return 0;
}
