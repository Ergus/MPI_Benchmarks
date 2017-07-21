#include "Node.h"

Node_slave::Node_slave(int &argc, char** &argv, MPI_Comm _parent):
  Node_t(argc,argv,_parent){
      
  MPI_Intercomm_merge(parent, true,  &intra);

  MPI_Comm_size(intra, &wsize);
  MPI_Comm_rank(intra, &wrank);
  }

Node_slave::~Node_slave(){  
  }

void Node_slave::listen(){
  int ret, flag, count;
	MPI_Status status;
  msg_t msg;

  fprintf(stderr, "Process %d listening\n", wrank);
  
  while(listening){
    // Prove the message first
    if ((ret = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, intra, &status))) {
      fprintf(stderr, "Error probing for messages\n");
      MPI_Abort(intra, ret);
      }
    
    if (!flag){
      fprintf(stderr, "Received flag was NULL\n");
      return;
      }

    // Determine size in bytes
    if ((ret = MPI_Get_count(&status, MPI_BYTE, &count))) {
      fprintf(stderr, "Error while determining message's size\n");
      MPI_Abort(intra, ret);
      }
    assert(count!=0);

    // Now receive the message
    if ((ret = MPI_Recv(&msg, count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, intra, MPI_STATUS_IGNORE))) {
      fprintf(stderr, "Error while receiving message\n");
      MPI_Abort(intra, ret);
      }

    fprintf(stderr, "%d<--(%d)--%d\n", wrank, msg.type, status.MPI_SOURCE);
    
    const msg_tag type=msg.type;
    
    if(type==TAG_SPAWN){
      fprintf(stderr,"Process %d: Spawning (world %d)\n", wrank, wsize);
      spawn_merge(msg.number);
      }
    else if(type==TAG_REDUCE){     
      fprintf(stderr,"Process %d: Reducing (world %d)\n", wrank, wsize);
      split_kill(msg.number);
      }
    else if(type==TAG_EXIT){
      fprintf(stderr,"Process %d: Exit listening (world %d)\n", wrank, wsize);
      listening=false;
      }
    }
  }

void Node_slave::run(){
  // Create a thread to listen incoming messages.
  std::thread tlisten(&Node_slave::listen, this);

  //Extra code goes here. Always before the join.

  tlisten.join();
  
  }
