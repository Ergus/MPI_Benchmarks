#include "matvec-mpi.h"

#include <string>

// Common local headers
#include "CommandLineParameter.hpp"
#include "Report.hpp"
#include "Reportable.hpp"
#include "Timer.hpp"

using namespace std;

int main(int argc, char** argv)
{
	CommandLine::initialize(argc, argv, "mpi_matmul", "cluster");
	CommandLineParameter<size_t> dim("dim", "Array size");
	CommandLineParameter<size_t> TS("TS", "Threads number");
	CommandLineParameter<size_t> its("its", "Number of iterations");
	OptionalCommandLineParameter<bool> print("print", 0, "Print Matrix");
	OptionalCommandLineParameter<string> pref("prefix", "matvec", "Prefix");
	CommandLine::validate();

	Initialize(&argc, &argv, dim, TS);
	printf("Initialized in process %d\n", _env.rank);

	//! _env.ldim is the fraction dim/N
	const size_t ldimA=_env.ldim*dim;

	const size_t dimA = (_env.IprintA ? dim : _env.ldim);
	const size_t dimC = (_env.IprintC ? dim : _env.ldim);

	printf("Allocating Memory in process %d\n", _env.rank);
	// Allocate memory
	double* A = (double*) malloc(dimA * dim * sizeof(double));
	double* B = (double*) malloc(dim * sizeof(double));
	double* C = (double*) calloc(dimC, sizeof(double));

	double* const lA = A;						  // C printer is always 0
	double* const lB = &B[_env.rank*_env.ldim];	  // B is fully distributed
	double* const lC = (_env.IprintC ? &C[_env.rank*_env.ldim] : C);

	Timer timer("algorithm_time", "Execution time");
	// Initialise arrays local portions
	init(lA, _env.ldim, dim);
	init(lB, _env.ldim, 1);

	printf("Call Allgather in process %d\n", _env.rank);

	if (_env.rank == 0)
		timer.start();
	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, _env.ldim, MPI_DOUBLE,
				  B, _env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
	dprintf("Multiplication %d\n", _env.rank);
	for(size_t i = 0; i < its; ++i)
		matvec(lA, B, lC, _env.ldim, dim);

	MPI_Barrier(MPI_COMM_WORLD);
	timer.stop();

	// Gather C to root
	dprintf("Call C Gather in process %d\n", _env.rank);

	if (_env.IprintC) {
		MPI_Gather(MPI_IN_PLACE, _env.ldim, MPI_DOUBLE,
				   C, _env.ldim, MPI_DOUBLE,
				   _env.rank, MPI_COMM_WORLD);
	} else {
		MPI_Gather(C, _env.ldim, MPI_DOUBLE,
				   NULL, _env.ldim, MPI_DOUBLE,
				   imin(2,_env.worldsize-1), MPI_COMM_WORLD);
	}

	// Gather A to its printer
	dprintf("Call A Gather in process %d\n", _env.rank);

	if (_env.IprintA) {
		MPI_Gather(MPI_IN_PLACE, _env.ldim*dim, MPI_DOUBLE,
				   A, _env.ldim*dim, MPI_DOUBLE,
				   0, MPI_COMM_WORLD);
	} else {
		MPI_Gather(lA, _env.ldim*dim, MPI_DOUBLE,
				   NULL, _env.ldim*dim, MPI_DOUBLE,
				   0, MPI_COMM_WORLD);
	}

	if (_env.rank==0)
		Report::emit();

	if (print){
	  printf("Printing Matrices in process %d\n", _env.rank);
	  const string prefix = pref;
	  printmatrix(A, _env.dim, _env.dim, prefix.c_str());
	  printmatrix(B, 1, _env.dim,prefix.c_str());
	  printmatrix(C, _env.dim, 1, prefix.c_str());
	}

	printf("Freeing Memory in process %d\n", _env.rank);

	free(C);
	free(B);
	free(A);

	Finalize();
	return 0;
}
