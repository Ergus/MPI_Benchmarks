#ifndef COMMAND_LINE_PARAMETER_HPP
#define COMMAND_LINE_PARAMETER_HPP


#include "CommandLine.hpp"
#include "ReportEntry.hpp"
#include "ReportEntryWithDefaultValue.hpp"

#include <sstream>


template <typename T>
class CommandLineParameter : public ReportEntry<T> {
private:
	bool _bad;
	
public:
	CommandLineParameter(std::string const &name, std::string const &description, std::string const &units = std::string())
		: ReportEntry<T>("PARAMETER", name, description, T(), units), _bad(false)
	{
		if (CommandLine::_currentIndex < CommandLine::_argc) {
			std::istringstream iss(CommandLine::_argv[CommandLine::_currentIndex]);
			iss >> ReportEntry<T>::_value;
			_bad = iss.bad();
		}
		CommandLine::_currentIndex++;
		
		CommandLine::registerMandatoryParameter(this);
	}
	
	inline bool bad() const
	{
		return _bad;
	}
};


template <typename T>
class OptionalCommandLineParameter : public ReportEntryWithDefaultValue<T> {
private:
	bool _bad;
	
public:
	OptionalCommandLineParameter(std::string const &name, T const &defaultValue, std::string const &description, std::string const &units = std::string())
		: ReportEntryWithDefaultValue<T>("PARAMETER", name, description, defaultValue, defaultValue, units),
		_bad(false)
	{
		if (CommandLine::_currentIndex < CommandLine::_argc) {
			std::istringstream iss(CommandLine::_argv[CommandLine::_currentIndex]);
			iss >> ReportEntryWithDefaultValue<T>::_value;
			_bad = iss.bad();
		}
		CommandLine::_currentIndex++;
		
		CommandLine::registerOptionalParameter(this);
	}
	
	template <typename T2>
	OptionalCommandLineParameter(std::string const &name, ReportEntry<T2> const &defaultValue, std::string const &description, std::string const &units = std::string())
		: ReportEntryWithDefaultValue<T>("PARAMETER", name, description, defaultValue, (T2) defaultValue, units),
		_bad(false)
	{
		if (CommandLine::_currentIndex < CommandLine::_argc) {
			std::istringstream iss(CommandLine::_argv[CommandLine::_currentIndex]);
			iss >> ReportEntryWithDefaultValue<T>::_value;
			_bad = iss.bad();
		}
		CommandLine::_currentIndex++;
		
		CommandLine::registerOptionalParameter(this);
	}
	
	template <typename T2>
	OptionalCommandLineParameter(std::string const &name, ReportEntryWithDefaultValue<T2> const &defaultValue, std::string const &description, std::string const &units = std::string())
		: ReportEntryWithDefaultValue<T>("PARAMETER", name, description, defaultValue, (T2) defaultValue, units),
		_bad(false)
	{
		if (CommandLine::_currentIndex < CommandLine::_argc) {
			std::istringstream iss(CommandLine::_argv[CommandLine::_currentIndex]);
			iss >> ReportEntryWithDefaultValue<T>::_value;
			_bad = iss.bad();
		}
		CommandLine::_currentIndex++;
		
		CommandLine::registerOptionalParameter(this);
	}
	
	inline bool bad() const
	{
		return _bad;
	}
};


#endif // COMMAND_LINE_PARAMETER_HPP
