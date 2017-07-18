#include "Node.h"

Manager::Manager(int &argc, char** &argv){
  int provided;

  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if(provided != MPI_THREAD_MULTIPLE){
    fprintf(stderr,"Error in MPI initialisation\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    }

  MPI_Comm_get_parent(&parent);  // Parent to decide

  if(parent==MPI_COMM_NULL)
    Node= new Node_master(argc,argv,parent);
  else
    Node= new Node_slave(argc,argv,parent);
  
  }

Manager::~Manager(){
  delete Node;
  MPI_Finalize();
  }
