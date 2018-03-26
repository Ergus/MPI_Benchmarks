#include "matvec-mpi-pthreads.h"
#include <pthread.h>
#include <assert.h>
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

	_env.dim = dim;					                 // dimension
	_env.ldim = dim / _env.worldsize;                // rows for local process

	_env.IprintA = (_env.rank == 0);						// prints A
	_env.IprintB = (imin(1,_env.worldsize-1) == _env.rank); // prints B
	_env.IprintC = (imin(2,_env.worldsize-1) == _env.rank); // prints C
}

void Finalize()
{
	MPI_Finalize();
}

void *init_thread(void *in)
{
	input1 *lin = (input1 *) in;

	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(lin->seed, &status);
	double rnd;

	for (size_t i = 0; i < lin->len; ++i) {
		drand48_r(&status, &rnd);
		lin->arr[i] = rnd;
	}
}

void init(double * const array,
          const size_t rows, const size_t cols)
{
	size_t i;
	const size_t lthreads = _env.lthreads;
	const size_t fullsize = rows * cols;
	const size_t thread_size = fullsize / lthreads;

	input1 *inputs = malloc(lthreads * sizeof(input1));
	pthread_t *threads = malloc(lthreads * sizeof(pthread_t));

	for (i = 0; i < lthreads; ++i) {
		size_t seed = _env.first_local_thread+i;

		inputs[i] = (input1){ &array[i*thread_size], thread_size, seed};
		pthread_create(&threads[i] , NULL, *init_thread, (void*) &inputs[i]);
	}

	for (i = 0; i < lthreads; ++i)
		pthread_join(threads[i],NULL);

	free(inputs);
	free(threads);
}


void *matvec_thread(void *in)
{
	input2 *lin = (input2 *) in;
	const size_t rows = lin->rows, cols = lin->cols;
	size_t i, j;

	for (i = 0; i < rows; ++i) {
		double sum=0.0;
		for (j = 0; j < cols; ++j) {
			sum += lin->A[i*cols+j] * lin->B[j];
		}
		lin->C[i] = sum;
	}
}

void matvec(const double * const __restrict__ A,
            const double * const __restrict__ B,
            double* const __restrict__ C,
            const size_t rowsA, const size_t colsA)
{
	size_t i, j;

	const size_t lthreads = _env.lthreads;
	const size_t threads_rows = rowsA / lthreads;

	assert(rowsA%lthreads == 0);

	input2 *inputs = malloc(lthreads * sizeof(input1));
	pthread_t *threads = malloc(lthreads * sizeof(pthread_t));

	for(i = 0; i < lthreads; ++i) {
		inputs[i] = (input2){ &A[i*threads_rows*colsA],
		                      B,
		                      &C[i*threads_rows],
		                      threads_rows,
		                      colsA, i};
		pthread_create(&threads[i] , NULL, *matvec_thread, (void*) &inputs[i]);
	}

	for (i = 0; i < lthreads; ++i)
		pthread_join(threads[i],NULL);

	free(inputs);
	free(threads);
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

