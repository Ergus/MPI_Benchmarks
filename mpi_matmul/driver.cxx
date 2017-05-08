
#include "matmul-mpi.h"

int main(int argc, char** argv){

    size_t dim, N;
    sscanf(argv[1], "%zu", &dim);
    sscanf(argv[2], "%zu", &N);
    
    Initialize(&argc, &argv, dim, N);
    
    const size_t nels = _env.ldim*dim;  //local elements to multiply
    
    const size_t dimA = (_env.IprintA ? dim : _env.ldim);
    const size_t dimB = dim;            // Everybody has full B
    const size_t dimC = (_env.IprintC ? dim : _env.ldim);

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
    
    // Gather B to all
    MPI_Allgather(MPI_IN_PLACE, nels, MPI_DOUBLE,
                  B, nels, MPI_DOUBLE, MPI_COMM_WORLD);

    // Multiplication
    matmul(lA,B,C,_env.ldim,dim);

    // Gather A to its printer
    if(_env.IprintA){
        MPI_Gather(MPI_IN_PLACE, nels, MPI_DOUBLE,
                   A, nels, MPI_DOUBLE,
                   _env.rank, MPI_COMM_WORLD);
        }
    else{
        MPI_Gather(A, nels, MPI_DOUBLE,
                   NULL, nels, MPI_DOUBLE,
                   min(2,_env.worldsize-1), MPI_COMM_WORLD);
        }

    // Gather C to root
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

    printmatrix(A, _env.dim, "matrix");
    printmatrix(B, _env.dim, "matrix");
    printmatrix(C, _env.dim, "matrix");

    free(C);
    free(B);
    free(A);
    
    Finalize();
    return 0;
    }
