#include "Node.h"


Node_master::Node_master(int &argc, char** &argv, MPI_Comm _parent):
  Node_t(argc,argv,_parent){

  assert(wrank==0);  
  }

void Node_master::run(){
  dprintf("Starting master\n");
  if(nargc>1)
    automatic();
  else
    interactive();
  }

Node_master::~Node_master(){
  dprintf("Killing all slaves.\n");
  stopall();        // Stop all remotes
  dprintf("Deleting master\n");
  }

int Node_master::send_to_remotes(msg_t msg){
  // Only send to remotes if there are any
  if(wsize>1){    
    const int remotes=wsize-1;
    // Send unblocking spawn messages to all remotes.
    MPI_Request *requests=(MPI_Request*)malloc(remotes*sizeof(MPI_Request));
    MPI_Status *statuses=(MPI_Status*)malloc(remotes*sizeof(MPI_Status));
    
    for(int i=1;i<wsize;++i){
      dprintf("%d--(%d)-->%d (world %d)\n",wrank,msg.type,i,wsize);
      MPI_Isend(&msg, sizeof(msg_t), MPI_BYTE, i, msg.type, intra, &requests[i-1]);
      }
    
    // Check that all sends were fine.
    if(MPI_ERR_IN_STATUS==MPI_Waitall(remotes, requests, statuses)){
      dprintf("Error sending unblocking message to all from %d\n",wrank);
      for(int i=0;i<remotes;++i){
        if(statuses[i].MPI_ERROR!=MPI_SUCCESS){
          dprintf("Error: %d Sending %d--(%d)-->%d (world %d) error\n",
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
  return MPI_SUCCESS;
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
    
// Parent specific code to reduce communicator
int Node_master::split_kill(size_t n){
  // first send message if there are remotes
  dprintf("Reducing %lu processes (world %d)\n", n, wsize);
  if(wsize>1 && n<wsize-1){
    dprintf("Sending reduce to remotes\n");
    msg_t msg_red={TAG_REDUCE,n};
    send_to_remotes(msg_red);
    }

  return Node_t::split_kill(n);
  }

void Node_master::stopall(){
  msg_t msg_top={TAG_EXIT,0};

  send_to_remotes(msg_top);

  }

void Node_master::automatic(){
  // Now parse the cl options
  int opt, value;  
  while ((opt = getopt(nargc, nargv, "c:d:p")) != -1) {
    switch(opt){
      case 'c':
        value=atoi(optarg);
        printf("Processing %c accepted %d\n",opt, value);
        spawn_merge(value);
        break;
      case 'd':
        value=atoi(optarg);
        printf("Processing %c accepted %d\n",opt, value);
        split_kill(value);
        break;
      case 'p':
        fflush(stderr);
        printf("Press enter to continue\n");
        getchar();
        break;
      case '?':
        dprintf("Option %c not recognised\n",opt);
        MPI_Abort(intra, MPI_ERR_OTHER);
      }
    }  
  }

// Dynamic session

const char *character_names[] = {
  "spawn",
  "delete",
  "exit",
  NULL
  };

char *character_name_generator(const char *text, int state){
  static int list_index, len;
  const char *name;

  if (!state) {
    list_index = 0;
    len = strlen(text);
    }

  while ((name = character_names[list_index++])) {
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
      }
    }

  return NULL;
  }

char **character_name_completion(const char *text, int start, int end){
  rl_attempted_completion_over = 1;
  return rl_completion_matches(text, character_name_generator);
  }

void Node_master::interactive(){

  char *line;
  char command[10];
  int value, scanned;

  rl_attempted_completion_function = character_name_completion;
  
  while((line = readline("> ")) != NULL) {
    printf("[%s]\n", line);
    if (*line){
      add_history(line);
      scanned=sscanf(line,"%s %d",command, &value);
      if(scanned==1){
        if(string_in(command,"exit","quit","q")){
          break;
          }
        }
      else if(scanned==2){
        if(string_in(command,"spawn","s")){
          printf("Spawning %d\n",value);
          spawn_merge(value);
          }        
        else if(string_in(command,"reduce","r")){
          printf("Reducing %d scanned %d\n",value, scanned);
          split_kill(value);
          }
        else{
          printf("Invalid command: %s\n",command);
          }
        }
      } 
    free(line);
    }
  
  }
