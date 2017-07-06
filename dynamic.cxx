#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <thread>

typedef std::vector<MPI_Comm> vCom;

typedef struct Node_t{
    vCom inter, intra;                   // vector with all inter and intra communicators
    int nargc;
    char** nargv;
    MPI_Comm parentcomm;                 // Parent communicator, in init time always
    int wsize, wrank;    
  } Node_t;

typedef struct msg_t{
    int value;
  } msg_t;

Node_t Node;

int spawn_merge_parent(const char* exec, char** args, int n){
  
  MPI_Comm newintra = MPI_COMM_NULL;
  MPI_Comm newinter = MPI_COMM_NULL;
  const int wsize=Node.wsize;
  const int wrank=Node.wrank;
  int errcode=-1;

  if((wrank==0) && (wsize>1)){
    msg_t msg={n};
    
    MPI_Request request[wsize-1];
    
    for(int i=1;i<wsize;++i){
      fprintf(stderr,"Process %d send %d to %d in world %d\n", wrank, msg.value,i,wsize);
      //MPI_Isend(&msg, sizeof(msg_t), MPI_BYTE, i, 5, Node.intra.back(), &request[i-1]);
      //MPI_Request_free(&request[i-1]);
      MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, 1, Node.intra.back());

      }

    //MPI_Waitall(wsize-1, request, MPI_STATUSES_IGNORE);
    }

  fprintf(stderr, "Process %d Spawning in world %d\n", wrank, wsize);
  MPI_Comm_spawn( exec, args, 1, MPI_INFO_NULL, 0, Node.intra.back(), &newinter, &errcode);
  MPI_Intercomm_merge(newinter, false, &newintra);

  Node.inter.push_back(newinter);
  Node.intra.push_back(newintra);
  MPI_Comm_size(newintra, &Node.wsize);

  fprintf(stderr,"Ending spawn in %d\n", Node.wrank);
  return errcode;
  }

void stopall(){
  msg_t msg={0};
  
  for(int i=Node.wsize-1; i>0 ;i--){
    fprintf(stderr,"Process: %d Sending stop to %d\n",Node.wrank,i);
    MPI_Send(&msg, sizeof(msg_t), MPI_BYTE, i, 5, Node.intra.back());
    }
  }
  

void listener(){
  int delta=0;
	int ret, flag, count;
	MPI_Status status;
  msg_t msg;

  fprintf(stderr, "Process %d listening on world %d\n", Node.wrank, Node.wsize);
  
  while(true){
    // Prove the message first
    if ((ret = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, Node.intra.back(), &status))) {
      fprintf(stderr, "Error probing for messages\n");
      MPI_Abort(Node.intra.back(), ret);
      }
    
    if (!flag){
      fprintf(stderr, "Received flag wass NULL\n");
      return;
      }

    // Determine size in bytes
    if ((ret = MPI_Get_count(&status, MPI_BYTE, &count))) {
      fprintf(stderr, "Error while determining message's size\n");
      MPI_Abort(Node.intra.back(), ret);
      }
    assert(count!=0);

    // Now receive the message
    if ((ret = MPI_Recv(&msg, count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, Node.intra.back(), MPI_STATUS_IGNORE))) {
      fprintf(stderr, "Error while receiving message\n");
      MPI_Abort(Node.intra.back(), ret);
      }

    delta=msg.value;

    fprintf(stderr, "Process %d received %d from %d\n", Node.wrank, delta, status.MPI_SOURCE);
    
    if(delta>0){
      fprintf(stderr,"Received %d in %d: Spawing\n",delta, Node.wrank);
      spawn_merge_parent(Node.nargv[0], &Node.nargv[1], delta);
      }
    else if(delta<0){
      fprintf(stderr,"Received %d in %d: Negative values not supported yet\n",delta, Node.wrank);
      }
    else if(delta==0){
      fprintf(stderr,"Received %d in %d: Exiting\n",delta, Node.wrank);
      break;
      }
    }
  }

// Start the main
int main( int argc, char *argv[] ){

  int provided;
  
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  
  MPI_Comm_get_parent( &Node.parentcomm ); // Parent
  
  if (Node.parentcomm == MPI_COMM_NULL){

    int nsp, errcodes;

    // Check that we only start one master, otherwise it'll be a problem
    MPI_Comm_size(MPI_COMM_WORLD, &Node.wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &Node.wrank);
    
    if((Node.wsize!=1) || (argc==1)){ // argc only for this example, remove when porting
      printf("Usage: mpirun -np 1 %s NUM_SPAWNS\n",argv[0]);    
      MPI_Abort(MPI_COMM_WORLD,1);
      }
 
    nsp = atoi(argv[1]);
    printf("Parent will spawn %d processes\n", nsp);

    Node.nargc=argc;
    Node.nargv=argv;
    Node.intra.push_back(MPI_COMM_WORLD);

    for(int i=0;i<nsp;++i){
      spawn_merge_parent(Node.nargv[0], &Node.nargv[1], 1);
      }

    stopall();
    
    }
  else{

    MPI_Comm intracomm = MPI_COMM_NULL;

    Node.nargc=argc;
    Node.nargv=argv;
    
    MPI_Intercomm_merge(Node.parentcomm, true,  &intracomm);
    
    MPI_Comm_size(intracomm, &Node.wsize);
    MPI_Comm_rank(intracomm, &Node.wrank);
    Node.inter.push_back(Node.parentcomm);
    Node.intra.push_back(intracomm);

    std::thread t1(listener);
    t1.join();
    }

  fflush(stdout);
  MPI_Finalize();

  return 0;
  }
