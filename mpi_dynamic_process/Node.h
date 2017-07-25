#ifndef NODE_HH
#define NODE_HH

#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <getopt.h>
#include <vector>
#include <thread>

#include "benchmarks.h"

#include "readline/readline.h"
#include "readline/history.h"

// Messages types
enum msg_tag{
  TAG_EXIT=0,               // Exit
  TAG_REDUCE,               // Delete processes
  TAG_SPAWN,                // Spawn add new mpi processes
  };

class Manager;

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
    
    MPI_Comm intra, parent;            // persistent communicators
    int wsize, wrank;                  // environment MPI vars 
    int nargc;                         // command line arguments number
    char** nargv;                      // command line arguments vars
    virtual int spawn_merge(size_t n); // spawns n new mpi processes
    virtual int split_kill(size_t n);  // split the actual communicator and kills
    bool listening;                    // process will listen initialised as true
    friend class Manager;              // The manager has direct access to node.
  };

class Node_master:public Node_t{
  public:
    Node_master(int &argc, char** &argv, MPI_Comm _parent);
    ~Node_master();
    
    void run() override;
  private:
    int send_to_remotes(msg_t msg);    // like a broadcast with unblocking messages 
    void stopall();                    // kill all remotes before exiting    
    int spawn_merge(size_t n) override;         
    int split_kill(size_t n) override;
    
    void automatic();                  // automatic execution (command line)
    void interactive();                // interactive execution (readline)    
  };

class Node_slave:public Node_t{
  public:
    Node_slave(int &argc, char** &argv, MPI_Comm _parent);
    ~Node_slave();
    void run() override;
  private:
    void listen();                     // listen function, running while listening==true
  };

class Manager{
  public:
    Manager(int &argc, char** &argv);
    ~Manager();

    int run();
    
  private:
    Node_t *Node;
    MPI_Comm parent;
    int local_wsize, local_wrank;
  };

#endif
