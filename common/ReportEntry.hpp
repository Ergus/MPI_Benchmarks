#ifndef REPORT_ENTRY_HPP
#define REPORT_ENTRY_HPP

#include <limits>
#include <ostream>
#include <string>
#include <typeinfo>

#include "Reportable.hpp"
#include "TypeName.hpp"



template <typename T>
class ReportEntry : public Reportable {
protected:
	T _value;
	std::string _units;
	
public:
	ReportEntry(std::string const &area, std::string const &name, std::string const &description, T value = T(), std::string const &units = std::string())
		: Reportable(area, name, description),
		_value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		return o << TypeName<T>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntry<T> &operator=(T value)
	{
		_value = value;
		return *this;
	}
	
	inline operator T() const
	{
		return _value;
	}
	
};


template <>
class ReportEntry<double> : public Reportable {
protected:
	double _value;
	std::string _units;
	
public:
	ReportEntry(std::string const &area, std::string const &name, std::string const &description, double value = 0.0, std::string const &units = std::string())
		: Reportable(area, name, description),
		_value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		o.precision(std::numeric_limits<double>::digits10 + 2);
		return o << TypeName<double>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntry<double> &operator=(double value)
	{
		_value = value;
		return *this;
	}
	
	inline operator double() const
	{
		return _value;
	}
	
};


template <>
class ReportEntry<float> : public Reportable {
protected:
	float _value;
	std::string _units;
	
public:
	ReportEntry(std::string const &area, std::string const &name, std::string const &description, float value = 0.0, std::string const &units = std::string())
		: Reportable(area, name, description),
		_value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		o.precision(std::numeric_limits<float>::digits10 + 2);
		return o << TypeName<float>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntry<float> &operator=(float value)
	{
		_value = value;
		return *this;
	}
	
	inline operator float() const
	{
		return _value;
	}
	
};


#endif // REPORT_ENTRY_HPP
