#include "Node.h"


Node_master::Node_master(int &argc, char** &argv, MPI_Comm _parent):
  Node_t(argc,argv,_parent){
      
  // Test that only one mpi process was started (general)
  // Test that at least one argument was given (specific)
  if(wsize!=1){ 
    printf("Usage: mpirun -np 1 %s NUM_SPAWNS\n",nargv[0]);    
    MPI_Abort(intra,1);
    }      
  }

void Node_master::run(){

  if(nargc==1){ 
    printf("Usage: mpirun -np 1 %s NUM_SPAWNS\n",nargv[0]);
    printf("missing NUM_SPAWNS\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  int nsp = atoi(nargv[1]);      // Read CL arguments.
  
  printf("Parent will spawn %d processes one by one\n", nsp);
  for(int i=0;i<nsp;++i){
    spawn_merge(1);
    }
  }

Node_master::~Node_master(){
  stopall();        // Stop all remotes
  }

int Node_master::send_to_remotes(msg_t msg){

  // Only send to remotes if there are any
  if(wsize>1){    
    const int remotes=wsize-1;
    // Send unblocking spawn messages to all remotes.
    MPI_Request *requests=(MPI_Request*)malloc(remotes*sizeof(MPI_Request));
    MPI_Status *statuses=(MPI_Status*)malloc(remotes*sizeof(MPI_Status));
    
    for(int i=1;i<wsize;++i){
      fprintf(stderr,"%d--(%d)-->%d (world %d)\n",wrank,msg.type,i,wsize);
      MPI_Isend(&msg, sizeof(msg_t), MPI_BYTE, i, msg.type, intra, &requests[i-1]);
      }
    
    // Check that all sends were fine.
    if(MPI_ERR_IN_STATUS==MPI_Waitall(remotes, requests, statuses)){
      fprintf(stderr,"Error sending unblocking message to all from %d\n",wrank);
      for(int i=0;i<remotes;++i){
        if(statuses[i].MPI_ERROR!=MPI_SUCCESS){
          fprintf(stderr,"Error: %d Sending %d--(%d)-->%d (world %d) error: %d\n",
                  statuses[i].MPI_ERROR,
                  statuses[i].MPI_SOURCE,
                  statuses[i].MPI_TAG,
                  i+1,
                  wsize);
          }
        }
      return 1;
      }
    // Release the arrays
    free(requests);
    free(statuses);
    }
  return 0;
  }

// Parent specific code to spawn
int Node_master::spawn_merge(size_t n){

  // first send message if there are remotes
  if(wsize>1){
    msg_t msg_spawn={TAG_SPAWN,n};
    send_to_remotes(msg_spawn);
    }

  
  return Node_t::spawn_merge(n);
  }
    
 
void Node_master::stopall(){
  msg_t msg_top={TAG_EXIT,0};

  send_to_remotes(msg_top);

  }
