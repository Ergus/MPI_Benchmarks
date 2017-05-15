#include "daxpy-mpi.h"

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
    CommandLineParameter<double> a("a", "a value for y=a*x+y", "double scalar value");
    CommandLineParameter<string> clprefix("prefix", "prefix", "prefix for files");
    CommandLine::validate();

    const string prefix=clprefix;   // This is needed don't know why

    printf("Initialization in process %d\n", _env.rank);
    Initialize(&argc, &argv, dim, N);

    const size_t ldim=_env.ldim;                                // ldim is very used
    const size_t dimX = (_env.IprintX ? dim : ldim);            // X printed in 0 or 1
    const size_t dimY = (_env.rank==0 ? dim : ldim);            // Y always printed in rank=0

    printf("In process: %d dimX=%lu dimY=%lu \n",_env.rank, dimX, dimY);
    printf("Allocating Memory in process %d\n", _env.rank);
    // Allocate memory
    double* X = (double*) malloc( dimX*sizeof(double));
    double* Y = (double*) malloc( dimY*sizeof(double));
    
    double* lX=(_env.IprintX ? &X[_env.rank*ldim] : X );

    // Initialise arrays local portions
    init(lX,ldim);
    init( Y,ldim);

    printf("First Gather for X to print before daxpy execution process %d\n",_env.rank);
    if(_env.rank==0){
        MPI_Gather(MPI_IN_PLACE, ldim, MPI_DOUBLE,
                   Y, ldim, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(Y, ldim, MPI_DOUBLE,
                   NULL, ldim, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }
    string tmppref=prefix+"_in";
    printmatrix(Y, _env.dim, tmppref.c_str());

    Timer timer("algorithm_time", "Execution time");
    
    if (_env.rank==0) timer.start();

    daxpy(Y,a,lX,ldim);  // remember this is Y+=a*X
    timer.stop();

    if(_env.rank==0){
        MPI_Gather(MPI_IN_PLACE, ldim, MPI_DOUBLE,
                   Y, ldim, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(Y, ldim, MPI_DOUBLE,
                   NULL, ldim, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }

    // Gather A to its printer
    printf("Gather for Y to print in process %d\n",_env.rank);    
    if(_env.IprintX){
        MPI_Gather(MPI_IN_PLACE, ldim, MPI_DOUBLE,
                   X, ldim, MPI_DOUBLE,
                   _env.rank, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(X, ldim, MPI_DOUBLE,
                   NULL, ldim, MPI_DOUBLE,
                   imin(1,_env.worldsize-1), MPI_COMM_WORLD);
        }
    
    if(_env.rank==0) Report::emit();    

    printf("Printing Matrices in process %d\n", _env.rank);
    tmppref=prefix+"_out";
    printmatrix(Y, _env.dim, tmppref.c_str());    
    printmatrix(X, _env.dim, prefix.c_str());

    printf("Freeing Memory in process %d\n", _env.rank);
    free(Y);
    free(X);
    
    Finalize();
    return 0;
    }
