#ifndef MULTITIMER_HPP
#define MULTITIMER_HPP


#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include <vector>
#include <utility>

#include <string.h>
#include <time.h>

using namespace std;

class Multitimer {
  private:
    string _name;
    bool _running;
    timespec _startTime, _accumulated, _temporal, _endTime;
    vector<pair<string,timespec>> _register;
    
    
    inline void getTime(struct timespec &ts){
      int rc = clock_gettime(CLOCK_MONOTONIC, &ts);
      if (rc != 0) {
        int error = errno;
        std::cerr << "Error reading time: " << strerror(error) << std::endl;
        exit(1);
        }
      }
	
    inline double getNS(struct timespec const &ts) const{
      return ts.tv_sec * 1000000000.0 + ts.tv_nsec;
      }

  public:
    Multitimer(string name="unnamed"):
      _name(name),_running(false){
      }

    inline void start(string name=""){
      if(_running){
        std::cerr<<"Timer already running: stop it first"<<std::endl;
        }
      else{
        if(!name.empty()) _name=name;
        _register.clear();
        _running=true;
        _accumulated.tv_sec=0;
        _accumulated.tv_nsec=0;
        getTime(_startTime);
        }
      }
    
    inline void fetch(string note){
      if(_running){
        pause();
        _register.push_back(std::make_pair(note,_temporal));
        _running=true;
        getTime(_startTime);        
        }
      }

    inline void stop(){
      pause();
      _register.push_back(std::make_pair("Total",_accumulated));
      }    
    
    inline void pause(){
      getTime(_endTime);      
      if(_running){
        _running=false;
        
        if (_endTime.tv_nsec < _startTime.tv_nsec) {
          _temporal.tv_nsec = 1000000000L + _endTime.tv_nsec - _startTime.tv_nsec;
          _temporal.tv_sec = _endTime.tv_sec - 1 - _startTime.tv_sec;        
          }
        else {
          _temporal.tv_nsec = _endTime.tv_nsec - _startTime.tv_nsec;
          _temporal.tv_sec = _endTime.tv_sec - _startTime.tv_sec;      
          }
        
        _accumulated.tv_nsec += _temporal.tv_nsec;
        _accumulated.tv_sec += _temporal.tv_sec;
        }
      else{
        std::cerr<<"Timer is not running: start it first to pause"<<std::endl;
        }
      }

    friend ostream& operator<<(ostream& os, Multitimer in){
      os<<in._name<<" ";
      for(auto a: in._register){        
        os<<a.first<<": "<<in.getNS(a.second)<<" ";
        }
      os<<"ns";
      return os;
      }    
  };

#endif // MULTITIMER_HPP
