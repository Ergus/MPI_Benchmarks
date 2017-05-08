#ifndef REPORT_HPP
#define REPORT_HPP


#include <cassert>
#include <list>
#include <iostream>
#include <ostream>

class Reportable;


class Report {
	std::list<Reportable *> _reportEntries;
	
	static Report _singleton;
	
public:
	Report();
	
	static inline void addEntry(Reportable *entry);
	static void emit(std::ostream &o = std::cout);
};


#include "Reportable.hpp"


inline void Report::addEntry(Reportable *entry)
{
	assert(entry != nullptr);
	_singleton._reportEntries.push_back(entry);
}


#endif // REPORT_HPP

