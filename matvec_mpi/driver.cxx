#include "matvec-mpi.h"

#include <iostream>
#include <string>

// Common local headers
#include "CommandLineParameter.hpp"
#include "Report.hpp"
#include "Reportable.hpp"
#include "Timer.hpp"

using namespace std;

int main(int argc, char** argv)
{
	CommandLine::initialize(argc, argv, "matmul_mpi", "cluster");
	CommandLineParameter<size_t> dim("dim", "Array size");
	CommandLineParameter<size_t> nth("nth", "Number of threads per node");
	CommandLineParameter<size_t> its("its", "Number of iterations");
	OptionalCommandLineParameter<bool> print("print", 0, "Print Matrices");
	OptionalCommandLineParameter<string> pref("prefix", "MTV", "Prefix");
	CommandLine::validate();

	cout << "Initialized in process" << _env.rank << endl;
	Initialize(&argc, &argv, dim, nth);

	const size_t ldimA=_env.ldim*dim;        //_env.ldim is the fraction dim/N
	const size_t dimA = (_env.IprintA ? dim : _env.ldim);
	const size_t dimC = (_env.IprintC ? dim : _env.ldim);

	Timer timer("algorithm_time", "Execution time");
	Timer total_time("total_time", "Total execution time");

	double* A = (double*) malloc(dimA * dim * sizeof(double));
	double* B = (double*) malloc(dim * sizeof(double));
	double* C = (double*) calloc(dimC, sizeof(double));

	double* const lA = A;						  // C printer is always 0
	double* const lB = &B[_env.rank*_env.ldim];	  // B is fully distributed
	double* const lC = (_env.IprintC ? &C[_env.rank*_env.ldim] : C);

	// Initialize arrays local portions
	init(lA, _env.ldim, dim);
	init(lB, _env.ldim, 1);

	//==========================================================================	
	if (_env.rank == 0){
		cout << "Starting algorithm" << endl;
		timer.start(); // =============================
	}
	// Gather B to all
	MPI_Allgather(MPI_IN_PLACE, _env.ldim, MPI_DOUBLE,
	              B, _env.ldim, MPI_DOUBLE, MPI_COMM_WORLD);

	// Multiplication
 	for(size_t i = 0; i < its; ++i)
		matvec(lA, B, lC, _env.ldim, _env.dim);

	timer.stop();
	// =========================================================================
	total_time.stop();

	MPI_Barrier(MPI_COMM_WORLD);

	cout << "Finished algorithm..." << endl;

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
	if (_env.IprintA) {
		MPI_Gather(MPI_IN_PLACE, _env.ldim*dim, MPI_DOUBLE,
		           A, _env.ldim*dim, MPI_DOUBLE,
		           0, MPI_COMM_WORLD);
	} else {
		MPI_Gather(lA, _env.ldim*dim, MPI_DOUBLE,
		           NULL, _env.ldim*dim, MPI_DOUBLE,
		           0, MPI_COMM_WORLD);
	}

	if (print) {
		const string prefix = pref;
		printmatrix(A, _env.dim, _env.dim, prefix.c_str());
		printmatrix(B, 1, _env.dim,prefix.c_str());
		printmatrix(C, _env.dim, 1, prefix.c_str());
		cout << "Done printing results..." << endl;
	}

	if (_env.rank==0){
		double performance = its * dim * dim * 2.0;
		performance = performance * 1000.0 / (double) timer; // In nanoseconds
		ReportEntry<double> passes("RESULT", "performance", "Megaflops per second", performance, "MFlops/s");

		Report::emit();
	}

	free(C);
	free(B);
	free(A);

	Finalize();
	return 0;
}
