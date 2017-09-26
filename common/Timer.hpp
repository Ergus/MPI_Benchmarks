#ifndef TIMER_HPP
#define TIMER_HPP


#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include <cstring>
#include <time.h>


#include "Reportable.hpp"


class Timer : Reportable {
	struct timespec _startTime;
	struct timespec _endTime;
	struct timespec _accumulated;
	
	inline void getTime(struct timespec &ts)
	{
		int rc = clock_gettime(CLOCK_MONOTONIC, &ts);
		if (rc != 0) {
			int error = errno;
			std::cerr << "Error reading time: " << strerror(error) << std::endl;
			exit(1);
		}
	}
	
	inline double getNS(struct timespec const &ts) const
	{
		return ts.tv_sec * 1000000000.0 + ts.tv_nsec;
	}
	
public:
	Timer(std::string const &name, std::string const &description)
		: Reportable("MEASUREMENT", name, description)
	{
		_startTime.tv_nsec = 0; _startTime.tv_sec = 0;
		_endTime.tv_nsec = 0; _endTime.tv_sec = 0;
		_accumulated.tv_nsec = 0; _accumulated.tv_sec = 0;
		start();
	}
	
	inline void start()
	{
		getTime(_startTime);
	}
	
	inline void stop()
	{
		getTime(_endTime);
		if (_endTime.tv_nsec < _startTime.tv_nsec) {
			long nsec = 1000000000L + _endTime.tv_nsec - _startTime.tv_nsec;
			_accumulated.tv_nsec += nsec;
			_accumulated.tv_sec += _endTime.tv_sec - 1 - _startTime.tv_sec;
		} else {
			_accumulated.tv_nsec += _endTime.tv_nsec - _startTime.tv_nsec;
			_accumulated.tv_sec += _endTime.tv_sec - _startTime.tv_sec;
		}
	}
	
	inline void reset()
	{
		_startTime.tv_nsec = 0;
		_startTime.tv_sec = 0;
		_endTime.tv_nsec = 0;
		_endTime.tv_sec = 0;
		_accumulated.tv_nsec = 0;
		_accumulated.tv_sec = 0;
	}
	
	inline bool hasBeenStartedAtLeastOnce()
	{
		return (_startTime.tv_nsec != 0L) || (_startTime.tv_sec != 0);
	}
	
	inline bool hasBeenStoppedAtLeastOnce()
	{
		return (_accumulated.tv_nsec != 0L) || (_accumulated.tv_sec != 0);
	}
	
	inline double lap()
	{
		stop();
		start();
		return getNS(_accumulated);
	}
	
	inline double lapAndReset()
	{
		stop();
		double current = getNS(_accumulated);
		reset();
		start();
		return current;
	}
	
	inline operator double()
	{
		return getNS(_accumulated);
	}
	
	inline operator long int()
	{
		return getNS(_accumulated);
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		o.precision(std::numeric_limits<double>::digits10 + 2);
		return o << "double\t" << _name << "\t" << getNS(_accumulated) << "\tns\t" << _description << std::endl;
	}
	
};


#endif // TIMER_HPP
