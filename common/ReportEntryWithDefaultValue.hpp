#ifndef REPORT_ENTRY_WITH_DEFAULT_VALUE_HPP
#define REPORT_ENTRY_WITH_DEFAULT_VALUE_HPP

#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "Reportable.hpp"
#include "ReportEntry.hpp"
#include "TypeName.hpp"


template <typename T>
class ReportEntryWithDefaultValue : public ReportableWithDefaultValue {
protected:
	bool _hasSymbolicDefaultValue;
	std::string _symbolicDefaultValue;
	T _defaultValue;
	T _value;
	std::string _units;
	
public:
	ReportEntryWithDefaultValue(std::string const &area, std::string const &name, std::string const &description, T defaultValue = T(), T value = T(), std::string const &units = std::string())
		: ReportableWithDefaultValue(area, name, description),
		_hasSymbolicDefaultValue(false), _defaultValue(defaultValue), _value(value), _units(units)
	{
	}
	
	template <typename T2>
	ReportEntryWithDefaultValue(std::string const &area, std::string const &name, std::string const &description, ReportEntryWithDefaultValue<T2> const &defaultValue, T value, std::string const &units = std::string())
		: ReportableWithDefaultValue(area, name, description),
		_hasSymbolicDefaultValue(true), _symbolicDefaultValue(defaultValue.getName()), _defaultValue((T2) defaultValue), _value(value), _units(units)
	{
	}
	
	template <typename T2>
	ReportEntryWithDefaultValue(std::string const &area, std::string const &name, std::string const &description, ReportEntry<T2> const &defaultValue, T value, std::string const &units = std::string())
		: ReportableWithDefaultValue(area, name, description),
		_hasSymbolicDefaultValue(true), _symbolicDefaultValue(defaultValue.getName()), _defaultValue((T2) defaultValue), _value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		return o << TypeName<T>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntryWithDefaultValue<T> &operator=(T value)
	{
		_value = value;
		return *this;
	}
	
	inline operator T() const
	{
		return _value;
	}
	
	std::string getDefaultValueAsString() const
	{
		if (_hasSymbolicDefaultValue) {
			return _symbolicDefaultValue;
		} else {
			std::ostringstream oss;
			oss << _defaultValue;
			return oss.str();
		}
	}
};


template <>
class ReportEntryWithDefaultValue<double> : public ReportableWithDefaultValue {
protected:
	double _defaultValue;
	double _value;
	std::string _units;
	
public:
	ReportEntryWithDefaultValue(std::string const &area, std::string const &name, std::string const &description, double defaultValue = 0.0, double value = 0.0, std::string const &units = std::string())
		: ReportableWithDefaultValue(area, name, description),
		_defaultValue(defaultValue), _value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		o.precision(std::numeric_limits<double>::digits10 + 2);
		return o << TypeName<double>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntryWithDefaultValue<double> &operator=(double value)
	{
		_value = value;
		return *this;
	}
	
	inline operator double() const
	{
		return _value;
	}
	
	std::string getDefaultValueAsString() const
	{
		std::ostringstream oss;
		oss.precision(std::numeric_limits<double>::digits10 + 2);
		oss << _defaultValue;
		return oss.str();
	}
	
};


template <>
class ReportEntryWithDefaultValue<float> : public ReportableWithDefaultValue {
protected:
	float _defaultValue;
	float _value;
	std::string _units;
	
public:
	ReportEntryWithDefaultValue(std::string const &area, std::string const &name, std::string const &description, float defaultValue = 0.0, float value = 0.0, std::string const &units = std::string())
		: ReportableWithDefaultValue(area, name, description),
		_defaultValue(defaultValue), _value(value), _units(units)
	{
	}
	
	std::ostream &emitReport(std::ostream &o)
	{
		emitArea(o);
		o.precision(std::numeric_limits<float>::digits10 + 2);
		return o << TypeName<float>::getName() << "\t" << _name << "\t" << _value << "\t" << _units << "\t" << _description << std::endl;
	}
	
	inline ReportEntryWithDefaultValue<float> &operator=(float value)
	{
		_value = value;
		return *this;
	}
	
	inline operator float() const
	{
		return _value;
	}
	
	std::string getDefaultValueAsString() const
	{
		std::ostringstream oss;
		oss.precision(std::numeric_limits<float>::digits10 + 2);
		oss << _defaultValue;
		return oss.str();
	}
	
};


#endif // REPORT_ENTRY_WITH_DEFAULT_VALUE_HPP
