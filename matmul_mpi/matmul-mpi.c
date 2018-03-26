
#include "matmul-mpi.h"

envinfo _env;

void Initialize(int *argc, char ***argv, size_t dim, size_t N)
{
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(_env.rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(_env.worldsize));

	myassert(dim>=N);              // more columns than threads
	myassert(N>=_env.worldsize);   // more threads than processes
	modcheck(dim, _env.worldsize); // columns % processes == 0
	modcheck(N, _env.worldsize);   // threads % processes == 0

	_env.dim=dim;
	_env.tthreads=N;

	const size_t lthreads = N/_env.worldsize;

	_env.lthreads =  lthreads;           // threads per process
	_env.first_local_thread = _env.rank * lthreads; // first local thread id

	_env.ldim = dim/_env.worldsize;     // rows for local process
	modcheck(_env.ldim, _env.lthreads); // ldim % lthreads == 0

	_env.IprintA = (imin(2,_env.worldsize-1)==_env.rank); // prints A
	_env.IprintB = (imin(1,_env.worldsize-1)==_env.rank); // prints B
	_env.IprintC = (0==_env.rank);                       // prints C

	omp_set_num_threads(lthreads);
}


void Finalize()
{
	MPI_Finalize();
}

void init(double* array, size_t rows, size_t cols)
{
	const size_t fullsize=rows*cols;
	const size_t start = _env.first_local_thread;

	#pragma omp parallel
	{
		const size_t id = omp_get_thread_num();
		size_t i;
		srand(start+id);

		#pragma omp for private(i)
		for(i=0; i<fullsize; ++i){
			array[i] = frand();
		}
	}
}

void matmul(double *A, double *B, double *C, int rows, int cols)
{
	size_t i,j,k;

	#pragma omp parallel for private(i,j,k)
	for(i=0; i<rows; ++i){
		for(j=0; j<cols; ++j){
			const double temp=A[i*cols+j];
			for(k=0; k<cols; ++k){
				C[i*cols+k]+= (temp*B[j*cols+k]);
			}
		}
	}
}

void __print(double* mat, size_t dim, const char *name, const char* prefix)
{
	int i,j;

	if(_env.rank<3){
		if ( ( _env.IprintA && (strcmp(name,"A")==0) ) ||
		     ( _env.IprintB && (strcmp(name,"B")==0) ) ||
		     ( _env.IprintC && (strcmp(name,"C")==0) )
			){

			printf("Printing %s in process %d\n",name,_env.rank);

			char filename[128];
			sprintf(filename,"%s_%s.mat", prefix, name);
			FILE* fp = fopen(filename, "w+");
			myassert(fp);

			fprintf(fp, "# name: %s\n", name);
			fprintf(fp, "# type: matrix\n");
			fprintf(fp, "# rows: %lu\n", dim);
			fprintf(fp, "# columns: %lu\n", dim);

			for(i=0; i<dim; ++i) {
				for(j=0; j<dim; ++j) {
					fprintf(fp, "%3.8lf ", mat[i*dim+j]);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}
	}
}
