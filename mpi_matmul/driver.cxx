#include "matmul-mpi.h"

#include <string>

// Common local headers
#include "CommandLineParameter.hpp"
#include "Report.hpp"
#include "Reportable.hpp"
#include "Timer.hpp"

using namespace std;
int main(int argc, char** argv){


    CommandLine::initialize(argc, argv, "mpi_matmul", "cluster");
    CommandLineParameter<long> dim("dim", "Array size", "number of elements");
    CommandLineParameter<long> N("N", "Threads number", "number of paralell threads");
    CommandLineParameter<string> clprefix("prefix", "prefix", "prefix for files");
    CommandLine::validate();

    const string prefix=clprefix;

    printf("Initialization in process %d\n", _env.rank);
    Initialize(&argc, &argv, dim, N);
    
    const size_t nels = _env.ldim*dim;  //local elements to multiply
    
    const size_t dimA = (_env.IprintA ? dim : _env.ldim);
    const size_t dimB = dim;            // Everybody has full B
    const size_t dimC = (_env.IprintC ? dim : _env.ldim);

    printf("In process: %d dimA=%lu dimB=%lu dimC=%lu\n",_env.rank, dimA, dimB, dimC);
    printf("Allocating Memory in process %d\n", _env.rank);
    // Allocate memory
    double* A = (double*) malloc( dimA*dim*sizeof(double));
    double* B = (double*) malloc( dimB*dim*sizeof(double));
    double* C = (double*) calloc( dimC*dim, sizeof(double));
    
    double* lA=(_env.IprintA ? &A[_env.rank*nels] : A );
    double* lB=&B[_env.rank*nels];  // B is fully distributed
    double* lC= C;                  // C printer is always 0

    // Initialise arrays local portions
    init(lA,_env.ldim,dim);
    init(lB,_env.ldim,dim);

    Timer timer("algorithm_time", "Execution time");

    printf("Call Allgather in process %d\n", _env.rank);
    if (_env.rank==0) timer.start();    
    // Gather B to all
    MPI_Allgather(MPI_IN_PLACE, nels, MPI_DOUBLE,
                  B, nels, MPI_DOUBLE, MPI_COMM_WORLD);

    // Multiplication
    printf("Multiplication %d\n", _env.rank);
    matmul(lA,B,C,_env.ldim,dim);

    // Gather C to root
    printf("Call C Gather in process %d\n", _env.rank);
    if(_env.IprintC){
        MPI_Gather(MPI_IN_PLACE, nels, MPI_DOUBLE,
                   C, nels, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(C, nels, MPI_DOUBLE,
                   NULL, nels, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }
    timer.stop();

    // Gather A to its printer
    printf("Call A Gather in process %d\n", _env.rank);
    if(_env.IprintA){
        MPI_Gather(MPI_IN_PLACE, nels, MPI_DOUBLE,
                   A, nels, MPI_DOUBLE,
                   _env.rank, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(A, nels, MPI_DOUBLE,
                   NULL, nels, MPI_DOUBLE,                       
                       imin(2,_env.worldsize-1), MPI_COMM_WORLD);
        }

    if(_env.rank==0) Report::emit();    

    printf("Printing Matrices in process %d\n", _env.rank);
    printmatrix(A, _env.dim, prefix.c_str());
    printmatrix(B, _env.dim, prefix.c_str());
    printmatrix(C, _env.dim, prefix.c_str());

    printf("Freeing Memory in process %d\n", _env.rank);

    free(C);
    free(B);
    free(A);
    
    Finalize();
    return 0;
    }
