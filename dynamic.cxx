#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <thread>

typedef struct Node_t{
    MPI_Comm intra, parent;           // communicators
    int wsize, wrank;                 // environment 
    int nargc;                             
    char** nargv;
  } Node_t;

enum msg_tag{
  TAG_EXIT=0,
  TAG_SPAWN=1,
  };

typedef struct msg_t{
    int value;
  } msg_t;

Node_t Node;

// This is called only in the master
int spawn_merge_parent(const char* exec, char** args, int n){
  
  MPI_Comm newintra = MPI_COMM_NULL;
  MPI_Comm newinter = MPI_COMM_NULL;
  const int wsize=Node.wsize;
  const int wrank=Node.wrank;
  int *errcode=(int*) malloc(n*sizeof(int));
  int success=0;

  // first send message if there are remotes
  if((wrank==0)&&(wsize>1)){
    msg_t msg={n};
    
    for(int i=1;i<wsize;++i){
      fprintf(stderr,"%d--(%d)-->%d (world %d)\n",wrank,msg.value,i,wsize);
      MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, TAG_SPAWN, Node.intra);
      }
    }

  // Spawn now
  fprintf(stderr, "Process %d Spawning in world %d\n", wrank, wsize);
  MPI_Comm_spawn(exec, args, n, MPI_INFO_NULL, 0, Node.intra, &newinter, errcode);
  for(int i=0; i<n; ++i){
    if(errcode[i]!=MPI_SUCCESS){
      fprintf(stderr,"Error: In %d spawning %d in world %d\n",wrank,i, wsize);
      continue;
      }
    ++success;
    }
  free(errcode);
  
  MPI_Intercomm_merge(newinter, false, &newintra);
  
  MPI_Comm_free(&Node.intra);

  Node.intra = newintra;
  MPI_Comm_size(newintra, &Node.wsize);  
  MPI_Comm_free(&newinter);
  fprintf(stderr,"Ending spawn in %d\n", Node.wrank);
  return success;
  }

void stopall(){
  msg_t msg={0};
  
  for(int i=Node.wsize-1; i>0 ;i--){
    fprintf(stderr,"%d--(%d)-->%d (world %d)\n",
            Node.wrank,msg.value,i,Node.wsize);
    MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, 5, Node.intra);
    }
  }
  

void listener(){
	int ret, flag, count;
	MPI_Status status;
  msg_t msg;

  fprintf(stderr, "Process %d listening\n", Node.wrank);
  
  while(true){
    // Prove the message first
    if ((ret = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, Node.intra, &status))) {
      fprintf(stderr, "Error probing for messages\n");
      MPI_Abort(Node.intra, ret);
      }
    
    if (!flag){
      fprintf(stderr, "Received flag wass NULL\n");
      return;
      }

    // Determine size in bytes
    if ((ret = MPI_Get_count(&status, MPI_BYTE, &count))) {
      fprintf(stderr, "Error while determining message's size\n");
      MPI_Abort(Node.intra, ret);
      }
    assert(count!=0);

    // Now receive the message
    if ((ret = MPI_Recv(&msg, count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, Node.intra, MPI_STATUS_IGNORE))) {
      fprintf(stderr, "Error while receiving message\n");
      MPI_Abort(Node.intra, ret);
      }

    fprintf(stderr, "Process %d<--(%d)--%d\n",
            Node.wrank, msg.value, status.MPI_SOURCE);
    
    int delta=msg.value;
    
    if(delta>0){
      fprintf(stderr,"Process %d: Spawing (world %d)\n",
              Node.wrank, Node.wsize);
      spawn_merge_parent(Node.nargv[0], &Node.nargv[1], delta);
      }
    else if(delta<0){
      fprintf(stderr,"Process %d: No action negative not supported yet\n",
              Node.wrank);
      }
    else if(delta==0){
      fprintf(stderr,"Process %d: Exit listening (world %d)\n",
              Node.wrank, Node.wsize);
      break;
      }
    }
  }


// Start the main
int main( int argc, char **argv ){

  int provided;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  Node.nargc=argc;
  Node.nargv=argv;
  
  if(provided != MPI_THREAD_MULTIPLE){
    fprintf(stderr,"Error in MPI initialisation\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    }
  
  MPI_Comm_get_parent( &Node.parent ); // Parent

  // Test that this is the first process started with mpirun not spawned
  if (Node.parent == MPI_COMM_NULL){

    // Check that we only start one master, otherwise it'll be a problem
    MPI_Comm_dup(MPI_COMM_WORLD, &Node.intra);
    MPI_Comm_size(Node.intra, &Node.wsize);
    MPI_Comm_rank(Node.intra, &Node.wrank);

    // Test that only one mpi process was started (general)
    // Test that at least one argument was given (specific)
    if((Node.wsize!=1) || (argc==1)){ 
      printf("Usage: mpirun -np 1 %s NUM_SPAWNS\n",argv[0]);    
      MPI_Abort(MPI_COMM_WORLD,1);
      }
 
    int nsp = atoi(argv[1]);      // Read CL arguments.

    // Start spawning using the function
    printf("Parent will spawn %d processes one by one\n", nsp);
    for(int i=0;i<nsp;++i){
      spawn_merge_parent(Node.nargv[0], &Node.nargv[1], 1);
      }

    // Stop all remotes
    stopall();
    
    }
  else{
    MPI_Intercomm_merge(Node.parent, true,  &Node.intra);
    
    MPI_Comm_size(Node.intra, &Node.wsize);
    MPI_Comm_rank(Node.intra, &Node.wrank);

    std::thread t1(listener);
    t1.join();
    }

  MPI_Comm_free(&Node.intra); 
  fprintf(stderr,"Process %d: Exit (world %d)\n",Node.wrank,Node.wsize);
  fflush(stdout);
  MPI_Finalize();

  return 0;
  }
