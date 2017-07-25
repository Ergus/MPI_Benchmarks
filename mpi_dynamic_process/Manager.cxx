#include "Node.h"

Manager::Manager(int &argc, char** &argv){
  
  int provided;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if(provided != MPI_THREAD_MULTIPLE){
    dprintf("Error in MPI initialisation\n");
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

  MPI_Comm_rank(MPI_COMM_WORLD, &local_wrank); // local world rank
  MPI_Comm_size(MPI_COMM_WORLD, &local_wsize); // local world size
  MPI_Comm_get_parent(&parent);                // Parent to decide

  if(parent==MPI_COMM_NULL && local_wrank==0)
    Node= new Node_master(argc,argv,parent);
  else
    Node= new Node_slave(argc,argv,parent);  
  }

int Manager::run(){
  Node->run(); 
  }

Manager::~Manager(){
  delete Node;
  MPI_Finalize();
  dprintf("Exiting %d in world %d\n",local_wrank,local_wsize);
  }


