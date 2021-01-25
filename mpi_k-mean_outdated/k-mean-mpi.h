#ifndef kmean_h
#define kmean_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "benchmarks.h"

//******* Some macros first ************************
#define frand(r)(2.*r*(double)rand()/(RAND_MAX)-r) //uniform rng in [-r,r]


/// This macro makes the memory allocation check and returns a verbose error
#define check(x){                                                               \
    if(!x){                                                                     \
      fprintf(stderr,"Error: Could not allocate memory or open file: %s ",#x);  \
      fprintf(stderr,"(Fun: %s in %s:%d)\n",                                    \
              __PRETTY_FUNCTION__, __FILE__, __LINE__);                         \
      exit(1);                                                                  \
      }                                                                         \
    }

//********* Helper debug function *******************
void dbgprintf(const char * format, ... ){
  #ifdef DEBUG
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
  #endif        
  }

//******** The auxiliar types ************************
#ifdef __cplusplus
extern "C" {
  #endif

  typedef struct{
    int rank, worldsize;
    size_t ncenters, npoints, tthreads;
    size_t lcenters, lpoints, lthreads;
    size_t first_center, first_point, first_thread;
    bool IprintX; // Y always printed in rank 0, but X depends
    } envinfo;

  extern envinfo _env;
  
  /// point struct
  typedef struct { double x, y; ///< are position axis.
      size_t near;              ///< when a point: represents the nearest center position;
                                ///< when a centre: represents the number of closest points
    } point;

  /// Arrays like struct, very simple
  typedef struct {
      const size_t size; 
      point* coords;      
    } pos;
    
  /// L2 norm squared to estimate the distance.
  inline double L2(point p, point c){
    const double x=p.x-c.x;
    const double y=p.y-c.y;
    return x*x+y*y;
    }

  void Initialize(int *argc, char ***argv,
                  size_t _npoints, size_t _ncenters, size_t _tthreads);
  void Finalize();
  pos* alloc_points(size_t dim);
  void init_points(pos* inout, size_t lfirst, size_t ldim, double radius);
  void free_points(pos* in);
  size_t find_close(pos* points, pos* centers);
  double move_centers(pos* points, pos* centers, pos* tmp);
  void kmean(pos* points, pos* centers, double error);
  void print_eps(pos* points, pos* centers, const char* filename);

#ifdef __cplusplus
}
#endif

#endif
