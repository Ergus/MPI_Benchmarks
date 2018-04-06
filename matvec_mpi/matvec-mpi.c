#include "matvec-mpi.h"

envinfo _env;

void Initialize(int *argc, char ***argv,
                const size_t dim, const size_t lthreads)
{
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(_env.rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(_env.worldsize));

	myassert(dim >= _env.worldsize);                 // more rows than task size
	myassert(lthreads > 0);                          // more rows than task size
	modcheck(dim, _env.worldsize);                   // we need to split exactly

	_env.lthreads = lthreads;		                 // threads per process
	_env.tthreads = lthreads * _env.worldsize;		 // Total threads
	_env.first_local_thread = _env.rank * lthreads;  // first local thread id
	omp_set_dynamic(0);                              // Disable dynamic teams
	omp_set_num_threads(lthreads);                   // parallel regions size

	_env.dim = dim;					                 // dimension
	_env.ldim = dim / _env.worldsize;                // rows for local process

	_env.IprintA = (_env.rank == 0);						// prints A
	_env.IprintB = (imin(1,_env.worldsize-1) == _env.rank); // prints B
	_env.IprintC = (imin(2,_env.worldsize-1) == _env.rank); // prints C

	// test that the cpuset > number of threads
	cpu_set_t mask;
	int ret = sched_getaffinity(0, sizeof(mask), &mask);
	myassert(ret == 0);

	_env.cpu_count = CPU_COUNT(&mask);
	myassert(_env.cpu_count >= lthreads);
}

void Finalize()
{
	MPI_Finalize();
}

void init(double * const __restrict__ array,
          const size_t rows, const size_t cols)
{
	size_t i;
	const size_t fullsize = rows * cols;
	#pragma omp parallel private (i)
	{
		srand(_env.first_local_thread + omp_get_thread_num());

		#pragma omp for
		for(i = 0; i < fullsize; ++i) {
			array[i] = frand();
		}
	}
}

void matvec(const double * const __restrict__ A,
            const double * const __restrict__ B,
            double * const __restrict__ C,
            const size_t rowsA, const size_t colsA)
{
	size_t i, j;

	#pragma omp parallel for private(i,j)
	for (i = 0; i < rowsA; ++i) {
		double sum=0.0;
		for (j = 0; j < colsA; ++j) {
			sum+=A[i*colsA+j]*B[j];
		}
		C[i]=sum;
	}
}

void __print(const double * const mat,
             const size_t rows, const size_t cols,
             const char *name, const char *prefix)
{
	int i,j;
	if (_env.rank<3) {
		if ( (_env.IprintA && (strcmp(name,"A") == 0)) ||
		     (_env.IprintB && (strcmp(name,"B") == 0)) ||
		     (_env.IprintC && (strcmp(name,"C") == 0))
		     ) {

			printf("Printing %s in process %d\n",name,_env.rank);

			char filename[128];
			sprintf(filename,"%s_%s.mat", prefix, name);
			FILE *fp = fopen(filename, "w+");
			myassert(fp);

			fprintf(fp, "# name: %s\n", name);
			fprintf(fp, "# type: matrix\n");
			fprintf(fp, "# rows: %lu\n", rows);
			fprintf(fp, "# columns: %lu\n", cols);

			for (i = 0; i < rows; ++i) {
				for(j = 0; j < cols; ++j) {
					fprintf(fp, "%3.8lf ", mat[i*cols+j]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}
	}
}

