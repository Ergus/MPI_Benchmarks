#include "Node.h"

Node_t::Node_t(int &argc, char** &argv, MPI_Comm _parent):
  nargc(argc),nargv(argv),parent(_parent){
  
  MPI_Comm_dup(MPI_COMM_WORLD, &intra);  // Creates a duplicated intracomm
  MPI_Comm_size(intra, &wsize);          // gets size in local world
  MPI_Comm_rank(intra, &wrank);          // gets rank in local world
  }

Node_t::~Node_t(){
  MPI_Comm_free(&intra);                  // Free the intra comm 
  fprintf(stderr,"Process %d: Exit (world %d)\n",wrank,wsize);
  fflush(stdout);

  }

//! Base spawn and merge
int Node_t::spawn_merge(size_t n){
  
  MPI_Comm newintra = MPI_COMM_NULL;           // Variable for intracomm
  MPI_Comm newinter = MPI_COMM_NULL;           // Temporal intercomm
  int *errcode=(int*) malloc(n*sizeof(int));   // Array for individual error report
  int success=0;                               // global error report

  // Spawn now by one.
  fprintf(stderr, "Process %d Spawning in world %d\n", wrank, wsize);
  MPI_Comm_spawn(nargv[0], &nargv[1], n, MPI_INFO_NULL, 0, intra, &newinter, errcode);
  for(int i=0; i<n; ++i){
    if(errcode[i]!=MPI_SUCCESS){
      fprintf(stderr,"Error: In %d spawning %d in world %d\n",wrank,i, wsize);
      continue;
      }
    ++success;
    }
  free(errcode);                                // Free error codes array
  MPI_Comm_free(&intra);                        // Free old intracomm before.
  
  MPI_Intercomm_merge(newinter, false, &newintra); // Create new intra

  intra = newintra;                              // Reassign the intra to the new one
  MPI_Comm_size(newintra, &wsize);               // update wsize
  MPI_Comm_free(&newinter);                      // Free the created intercomm
  fprintf(stderr,"Ending spawn in %d\n", wrank); // Delete this or add in verbose only
  return success;
  }
