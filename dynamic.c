#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <vector>


int spawn_merge(const char* exec, char** args, MPI_Comm actual,MPI_Comm *inter, MPI_Comm *intra){
  int errcode=-1;
  MPI_Comm_spawn( exec, args, 1, MPI_INFO_NULL, 0, actual, inter, &errcode);
  
  MPI_Intercomm_merge(*inter, 0, intra);
  return errcode;
  }

int main( int argc, char *argv[] ){

  int nsp=0;               // Total expawned processes
  int spawned_order;
  int *errcodes=NULL;

  int size;
  
  MPI_Comm parentcomm;     // Parent communicator
  MPI_Comm *intercomm;     // Inter communicator for every remote  
  MPI_Comm intracomm;      // Intra communicator to be expawned

  MPI_Init( &argc, &argv );

  if(argc==1){
    printf("Usage: %s NUM_SPAWNS\n",argv[0]);
    
    MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  MPI_Comm_get_parent( &parentcomm ); // Parent
  
  if (parentcomm == MPI_COMM_NULL){
    nsp = atoi(argv[1]);
    printf("Parent will spawn %d processes\n", nsp);
    
    spawned_order=0;
    errcodes=(int*)malloc(nsp*sizeof(int));
    intercomm=(MPI_Comm*)malloc(nsp*sizeof(MPI_Comm));

    char argument[16];
    char *myargv[2];
    myargv[0]=argument;
    myargv[1]=NULL;

    for(int i=0;i<nsp;++i){      
      printf("Parent, spawning loop -> %d\n", i);
      sprintf(argument, "%d", i+1);
      
      printf("Parent, spawning loop <- %d\n", i);
      
      MPI_Intercomm_merge(intercomm[i], 0, &intracomm);
      }
    
    }
  else{
    nsp=-1;    
    spawned_order = atoi(argv[1]);    
    intercomm=(MPI_Comm*)malloc(sizeof(MPI_Comm));    
    printf("I'm the child %d\n",spawned_order);
    *intercomm = parentcomm;

    MPI_Intercomm_merge(*intercomm, 1, &intracomm);

    MPI_Comm_size(intracomm, &size);    
    printf("Init size %d \n",size);
        
    }

//  // Print off a hello world message
//  int old_rank, old_size, new_rank, new_size;
//  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &old_size);
//  MPI_Comm_rank(newintracomm, &new_rank);
//  MPI_Comm_size(newintracomm, &new_size);
//  printf("Hello from rank %d of %d (was %s %d of %d)\n",
//         new_rank, new_size,
//         was_spawned ? "child" : "parent",
//         old_rank, old_size);
//  MPI_Barrier(newintracomm);
//
  fflush(stdout);
  MPI_Finalize();

  free(errcodes);
  free(intercomm);
  
  return 0;
  }
