#include "Node.h"

Node_t::Node_t(int &argc, char** &argv, MPI_Comm _parent):
  nargc(argc),nargv(argv),parent(_parent){
  
  MPI_Comm_dup(MPI_COMM_WORLD, &intra);  // Creates a duplicated comm
  MPI_Comm_size(intra, &wsize);          // gets size in local world
  MPI_Comm_rank(intra, &wrank);          // gets rank in local world
  }

Node_t::~Node_t(){
  MPI_Comm_free(&intra); 
  fprintf(stderr,"Process %d: Exit (world %d)\n",wrank,wsize);
  fflush(stdout);
  MPI_Finalize();
  }

int Node_t::spawn_merge(int n){
  
  MPI_Comm newintra = MPI_COMM_NULL;
  MPI_Comm newinter = MPI_COMM_NULL;
  int *errcode=(int*) malloc(n*sizeof(int));
  int success=0;

  // Spawn now
  fprintf(stderr, "Process %d Spawning in world %d\n", wrank, wsize);
  MPI_Comm_spawn(nargv[0], &nargv[1], n, MPI_INFO_NULL, 0, intra, &newinter, errcode);
  for(int i=0; i<n; ++i){
    if(errcode[i]!=MPI_SUCCESS){
      fprintf(stderr,"Error: In %d spawning %d in world %d\n",wrank,i, wsize);
      continue;
      }
    ++success;
    }
  free(errcode);
  MPI_Comm_free(&intra);
  
  MPI_Intercomm_merge(newinter, false, &newintra);

  intra = newintra;
  MPI_Comm_size(newintra, &wsize);  
  MPI_Comm_free(&newinter);
  fprintf(stderr,"Ending spawn in %d\n", wrank);
  return success;
  }
