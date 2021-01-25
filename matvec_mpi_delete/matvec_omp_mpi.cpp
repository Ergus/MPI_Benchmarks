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

#include "matvec.h"

#include "ArgParserC/argparser.h"

#include <iostream>

void matvec_omp(const double *A, const double *b, double *x,
                size_t rowsA, size_t colsA)
{
	#pragma omp parallel for
	for (size_t i = 0; i < rowsA; ++i) {
		x[i] = 0.0;
		for (size_t j = 0; j < colsA; ++j)
			x[i] += A[i * colsA + j] * b[j];
	}
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const int ROWS = create_cl_int("Rows");
	const int LTHREADS = create_cl_int("lthreads");
	const int ITS = create_optional_cl_int("iterations", 1);
	const int PRINT = create_optional_cl_int("print", 0);

	Initialize(&argc, &argv, ROWS, LTHREADS);

	std::cout << "Initialized in process " << _env.rank << std::endl;

	const size_t dimA = _env.IprintA ? _env.dim : _env.ldim;
	const size_t dimB = _env.dim;	// Everybody has full B
	const size_t dimC = _env.IprintC ? _env.dim : _env.ldim;

	timer *ttimer = create_timer("Total time");

	double *A = (double *)malloc(dimA * _env.dim * sizeof(double));
	double *B = (double *)malloc(dimB * sizeof(double));
	double *C = (double *)calloc(dimC, sizeof(double));

	double *const lA = (_env.IprintA ? &A[_env.rank * nels] : A); // Just in case.
	double *const lB = &B[_env.rank * _env.ldim];	              // B is fully distributed
	double *const lC = (_env.IprintC ? &C[_env.rank * _env.ldim] : C);

	// Initialize arrays local portions
	init(lA, _env.ldim, _env.dim);
	init(lB, _env.ldim, 1);

	//==========================================================================

	timer *atimer = create_timer("Algorithm time");

	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, _env.ldim, MPI_DOUBLE,
	              B, _env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
	for (size_t i = 0; i < ITS; ++i)
		matvec_omp(lA, B, lC, _env.ldim, _env.dim);

	MPI_Barrier(MPI_COMM_WORLD);

	stop_timer(atimer);
	// =========================================================================

	std::cout << "Finished algorithm..." << std::endl;

	if (_env.rank == 0) {
		const double performance = ITS * _env.dim * _env.dim * 2000.0 / getNS_timer(atimer);

		create_reportable_double("performance", performance);
		report_args();
	}

	if (PRINT) {
		// Gather C to ITS printer
		std::cout << "Call C Gather in process %d\n" << _env.rank << std::endl;
		if (_env.IprintC) {

			MPI_Gather(MPI_IN_PLACE, _env.ldim, MPI_DOUBLE,
			           C, _env.ldim, MPI_DOUBLE, _env.rank, MPI_COMM_WORLD);

			printmatrix(C, _env.dim, 1, "matvec");
		} else {
			MPI_Gather(C, _env.ldim, MPI_DOUBLE,
			           NULL, _env.ldim, MPI_DOUBLE,
			           imin(2, _env.worldsize - 1), MPI_COMM_WORLD);
		}

		// Gather A to ITS printer
		std::cout << "Call A Gather in process %d\n" << _env.rank << std::endl;
		if (_env.IprintA) {
			MPI_Gather(MPI_IN_PLACE, _env.ldim * _env.dim, MPI_DOUBLE,
			           A, _env.ldim * _env.dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			printmatrix(A, _env.dim, _env.dim, "matvec");
		} else {
			MPI_Gather(lA, _env.ldim * _env.dim, MPI_DOUBLE,
			           NULL, _env.ldim * _env.dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		// Print B not need to gather.
		if (_env.IprintB) {
			std::cout << "Print B in process %d: " << _env.rank << std::endl;
			printmatrix(B, 1, _env.dim, "matvec");
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (_env.rank)
			std::cout << "Done printing results..." << std::endl;
	}

	free(C);
	free(B);
	free(A);

	Finalize();
	return 0;
}
