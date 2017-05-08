#ifndef REPORTABLE_HPP
#define REPORTABLE_HPP

#include <ostream>
#include <string>


class Reportable {
protected:
	std::string _area;
	std::string _name;
	std::string _description;
	
	inline std::ostream &emitArea(std::ostream &o);
	
public:
	inline Reportable(std::string const &area, std::string const &name, std::string const &description);
	
	virtual ~Reportable()
	{
	}
	
	inline std::string const &getName() const
	{
		return _name;
	}
	
	inline std::string const &getDescription() const
	{
		return _description;
	}
	
	
	// Format: <area>	<type>	<name>	<value>	<units>	<optional long description>
	virtual std::ostream &emitReport(std::ostream &o) = 0;
};


class ReportableWithDefaultValue: public Reportable {
public:
	inline ReportableWithDefaultValue(std::string const &area, std::string const &name, std::string const &description)
		: Reportable(area, name, description)
	{
	}
	
	virtual std::string getDefaultValueAsString() const = 0;
};



#include "Report.hpp"


Reportable::Reportable(std::string const &area, std::string const &name, std::string const &description)
	: _area(area), _name(name), _description(description)
{
	Report::addEntry(this);
}


std::ostream& Reportable::emitArea(std::ostream& o)
{
	o << _area << "\t";
}


#endif // REPORTABLE_HPP
