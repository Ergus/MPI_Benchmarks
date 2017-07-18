#ifndef NODE_HH
#define NODE_HH

#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <thread>

// Messages types
enum msg_tag{    
  TAG_REDUCE=-1,            // Delete processes
  TAG_EXIT=0,               // Exit
  TAG_SPAWN=1,              // Spawn add new mpi processes
  };

// Simplest message struct
typedef struct msg_t{     
    msg_tag type;           // Message Tag value
    size_t number;          // Number of processes
  } msg_t;

// Base Node class for master and slave
class Node_t{
  public:
    virtual ~Node_t();
    virtual void run()=0;
    
  protected:
    Node_t(int &argc, char** &argv, MPI_Comm _parent);
    
    MPI_Comm intra, parent;     // persistent communicators
    int wsize, wrank;           // environment MPI vars 
    int nargc;                  // command line arguments number
    char** nargv;               // command line arguments vars
    virtual int spawn_merge(size_t n); //
  };

class Node_master:public Node_t{
  public:
    Node_master(int &argc, char** &argv, MPI_Comm _parent);
    ~Node_master();
    void run();
    int send_to_remotes(msg_t msg);
  private:
    void stopall();
    int spawn_merge(size_t n);
    
  };

class Node_slave:public Node_t{
  public:
    Node_slave(int &argc, char** &argv, MPI_Comm _parent);
    ~Node_slave();
    void run();
  private:
    void listen();
  };

class Manager{
  public:
    Manager(int &argc, char** &argv);
    ~Manager();

    Node_t *Node;
    MPI_Comm parent;
  };

#endif
