#include "Node.h"

Node_t::Node_t(int &argc, char** &argv, MPI_Comm _parent):
  nargc(argc),nargv(argv),parent(_parent),listening(true),
  hostid(gethostid()),total_info(NULL){
  
  MPI_Comm_dup(MPI_COMM_WORLD, &intra);              // Creates a duplicated intracomm
  MPI_Comm_size(intra, &wsize);                      // gets size in local world
  MPI_Comm_rank(intra, &wrank);                      // gets rank in local world

  MPI_Info_create(&info);
  MPI_Info_set(info, "pernode", "true");

  }

Node_t::~Node_t(){
  MPI_Info_free(&info);
  MPI_Comm_free(&intra);                           // Free the intra comm
  dprintf("Process %d: Exit (world %d)\n",wrank,wsize);  
  fflush(stdout);
  }

//! Base spawn and merge
int Node_t::spawn_merge(size_t n){
  
  MPI_Comm newintra = MPI_COMM_NULL;               // Variable for intracomm
  MPI_Comm newinter = MPI_COMM_NULL;               // Temporal intercomm
  int *errcode=(int*) malloc(n*sizeof(int));       // Array for individual error report
  timer.fetch("Malloc");
  // Spawn now n processes  
  dprintf("Process %d Spawning in world %d\n", wrank, wsize);
  const int success = MPI_Comm_spawn(nargv[0], &nargv[1], n, MPI_INFO_NULL,
                               0, intra, &newinter, errcode);
  if(success==MPI_ERR_SPAWN){
    for(int i=0; i<n; ++i){
      if(errcode[i]!=MPI_SUCCESS){
        dprintf("Error: In %d spawning %d in world %d\n",wrank,i, wsize);
        continue;
        }
      }
    }
  timer.fetch("Spawn");
  free(errcode);                                   // Free error codes array
  MPI_Comm_free(&intra);                           // Free old intracomm before.
  timer.fetch("Free");
  MPI_Intercomm_merge(newinter, false, &newintra); // Create new intra
  timer.fetch("Merge");
  intra = newintra;                                // Reassign the intra to the new one
  MPI_Comm_size(newintra, &wsize);                 // update wsize
  MPI_Comm_free(&newinter);                        // Free the created intercomm
  dprintf("Ending spawn in %d\n", wrank);          // Delete this or add in verbose only
  return success;
  }

//! Base split and kill
int Node_t::split_kill(size_t n){

  MPI_Comm newintra = MPI_COMM_NULL;               // Variable for intracomm
    
  dprintf("Process %d reducing %zu processes (world %d)\n", wrank, n, wsize);
  listening = (wrank<(wsize-n));  
  dprintf("Process %d split listening=%d (world %d)\n", wrank, listening, wsize);
  
  MPI_Comm_split(intra, (int)listening, wrank, &newintra);
  timer.fetch("Split");
  
  MPI_Comm_free(&intra);                           // Free old intracomm before.
  intra = newintra;                                // Reassign the intra to the new one

  MPI_Comm_size(newintra, &wsize);                 // update wsize
  dprintf("Process %d reduced (world %d) \n",wrank,wsize);
  timer.fetch("Free");
  return MPI_SUCCESS;  
  }


int Node_t::getinfo(){
  info_t linfo={wrank,hostid};
  
  MPI_Gather(&linfo, sizeof(info_t), MPI_BYTE,
             total_info, sizeof(info_t), MPI_BYTE,
             0, intra);  
  }
