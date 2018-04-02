
#include "k-mean-mpi.h"

envinfo _env;

void Initialize(int *argc, char ***argv,
                size_t _npoints, size_t _ncenters, size_t _tthreads){

    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(_env.rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(_env.worldsize));

    myassert(_tthreads>=_env.worldsize);   // more threads than processes
    myassert(_ncenters>=_tthreads);        // more centres than threads
    myassert(_npoints>=_tthreads);         // more points than threads
    
    modcheck(_tthreads, _env.worldsize);   // threads % processes == 0
    modcheck(_ncenters, _tthreads);        // centres % threads == 0
    modcheck(_npoints, _tthreads);         // points % threads == 0    


    _env.npoints=_npoints;
    _env.ncenters=_ncenters;
    _env.tthreads=_tthreads;
    
    const size_t lthreads = _env.tthreads/_env.worldsize;

    omp_set_num_threads(lthreads);                  // this is an empty macro if no OpenMP
    
    _env.lthreads =  lthreads;                      // threads per process
    _env.first_thread = _env.rank * lthreads;       // first local thread id
        
    _env.lcenters = _ncenters/_env.worldsize;       // centres for local process
    modcheck(_env.lcenters, _env.lthreads);         // lcentres % lthreads == 0
    _env.first_center=_env.lcenters*_env.rank;      // first center for this rank

    _env.lpoints = _npoints/_env.worldsize;         // points for local process
    modcheck(_env.lpoints, _env.lthreads);          // lpoints % lthreads == 0
    _env.first_point=_env.lpoints*_env.rank;       // first center for this rank
  }

void Finalize(){
    MPI_Finalize();
    }

/// Allocates and point vector
/** \param [in] dim array dimension **/
pos* alloc_points(size_t dim){
  pos *out=(pos*) malloc(sizeof(pos)); // Allocates memory
  check(out);                          // checks for proper allocation

  out->coords = (point*) malloc(sizeof(point)*dim);    
  check(out->coords);                  // checks for proper allocation

  // Sets struct
  *(size_t *)&out->size=dim;
  
  return out;
  }

void init_points(pos* inout, size_t lfirst, size_t ldim, double r){
  for(size_t i=lfirst; i<lfirst+ldim; ++i) {
    inout->coords[i].x=frand(r);
    inout->coords[i].y=frand(r);
    inout->coords[i].near = 0;
    }
  }

void free_points(pos* in){
  free(in->coords);
  free(in);
  }

void __print(double* mat, size_t dim, const char *name, const char* prefix){

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
        
            for(int i=0; i<dim; ++i) {
                fprintf(fp, "%3.8lf ", mat[i]);
                }
            printf("\n");
            fclose(fp);
            }
        }
    }

