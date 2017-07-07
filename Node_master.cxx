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

int Node_master::spawn_merge(int n){

  // first send message if there are remotes
  if(wsize>1){
    msg_t msg={n};
    
    for(int i=1;i<wsize;++i){
      fprintf(stderr,"%d--(%d)-->%d (world %d)\n",wrank,msg.value,i,wsize);
      MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, TAG_SPAWN, intra);
      }
    }

  return Node_t::spawn_merge(n);
  }
    
 
void Node_master::stopall(){
  msg_t msg={0};
  
  for(int i=wsize-1; i>0 ;i--){
    fprintf(stderr,"%d--(%d)-->%d (world %d)\n", wrank, msg.value, i, wsize);
    MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, 5, intra);
    }
  }
