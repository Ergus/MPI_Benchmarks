#ifndef NODE_HH
#define NODE_HH

#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <thread>

enum msg_tag{
  TAG_EXIT=0,
  TAG_SPAWN=1,
  };

typedef struct msg_t{
    int value;
  } msg_t;

class Node_t{
  public:
    virtual ~Node_t();
    virtual void run()=0;
    
  protected:
    Node_t(int &argc, char** &argv, MPI_Comm _parent);
        
    MPI_Comm intra, parent;     // communicators
    int wsize, wrank;           // environment 
    int nargc;                             
    char** nargv;
    virtual int spawn_merge(int n);
  };

class Node_master:public Node_t{
  public:
    Node_master(int &argc, char** &argv, MPI_Comm _parent);
    ~Node_master();
    void run();
  private:
    void stopall();
    int spawn_merge(int n);
    
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
