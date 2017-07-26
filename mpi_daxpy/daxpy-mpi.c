
#include "daxpy-mpi.h"

envinfo _env;

void Initialize(int *argc, char ***argv, size_t dim, size_t N){

    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(_env.rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(_env.worldsize));

    myassert(dim>=N);              // more columns than threads
    myassert(N>=_env.worldsize);   // more threads than processes
    modcheck(dim, _env.worldsize); // columns % processes == 0
    modcheck(N, _env.worldsize);   // threads % processes == 0
    
    _env.dim=dim;

    #ifdef _OPENMP
    _env.tthreads=N;
    printf("Con OpenMP\n");   
    #else
    _env.tthreads=_env.worldsize;
    printf("Sin OpenMP\n");
    #endif
    
    const size_t lthreads = _env.tthreads/_env.worldsize;

    omp_set_num_threads(lthreads);                  // this is an empty macro if no OpenMP
    
    _env.lthreads =  lthreads;                      // threads per process
    _env.first_local_thread = _env.rank * lthreads; // first local thread id
        
    _env.ldim = dim/_env.worldsize;                 // elements for local process
    modcheck(_env.ldim, _env.lthreads);             // ldim % lthreads == 0

    _env.IprintX = (imin(1,_env.worldsize-1)==_env.rank); // prints Y
        
    }


void Finalize(){
    MPI_Finalize();
    }

void init(double* array, size_t ldim){

    const size_t start = _env.first_local_thread;
    size_t i;
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        const size_t id = omp_get_thread_num();
        srand(start+id);
        #ifdef _OPENMP
        #pragma omp for private(i)
        #endif
        for(i=0; i<ldim; ++i){
            array[i] = frand();
            }        
        }
    }

void daxpy(double *lY,double a,double *X,size_t ldim){  // remember this is Y+=a*X
  size_t i;
  #ifdef _OPENMP
  #pragma omp parallel for private(i)
  #endif
  for(i=0;i<ldim;++i){
    lY[i]+= a*X[i];
    }
  }

void __print(double* mat, size_t dim, const char *name, const char* prefix){

  size_t i;
  if(_env.rank<3 && dim<=512){
    if ( ( (_env.rank==0) && (strcmp(name,"Y")==0) ) ||
         (  _env.IprintX  && (strcmp(name,"X")==0) )
         ){

      printf("Printing %s in process %d\n",name,_env.rank);
            
      char filename[128];
      sprintf(filename,"%s_%s.mat", prefix, name);
      FILE* fp = fopen(filename, "w+");            
      myassert(fp);

      fprintf(fp, "# name: %s\n", name);
      fprintf(fp, "# type: matrix\n");
      fprintf(fp, "# rows: 1\n");
      fprintf(fp, "# columns: %lu\n", dim);
        
      for(i=0; i<dim; ++i) {
        fprintf(fp, "%3.8lf ", mat[i]);
        }
      printf("\n");
      fclose(fp);
      }
    }
  }

