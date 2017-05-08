#include "CommandLine.hpp"
#include "CommandLineParameter.hpp"
#include "ReportEntry.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>


int CommandLine::_argc;
char **CommandLine::_argv;
int CommandLine::_currentIndex;
std::list<int> CommandLine::_failedIndexes;
int CommandLine::_totalMandatoryParameters;
int CommandLine::_totalOptionalParameters;
std::list<Reportable *> CommandLine::_mandatoryParameters;
std::list<ReportableWithDefaultValue *> CommandLine::_optionalParameters;
ReportEntry<std::string> *CommandLine::_benchmarkName = nullptr;
ReportEntry<std::string> *CommandLine::_executableType = nullptr;
ReportEntry<std::string> *CommandLine::_variant = nullptr;
// ReportEntry<std::string> *CommandLine::_executable = nullptr;
ReportEntry<std::string> *CommandLine::_commandLine = nullptr;

ReportEntry<std::string> *CommandLine::_mcxxVersion = nullptr;
ReportEntry<std::string> *CommandLine::_compilerFlags = nullptr;
ReportEntry<std::string> *CommandLine::_compilerVersion = nullptr;

ReportEntry<std::string> *CommandLine::_nanos6Version = nullptr;
ReportEntry<std::string> *CommandLine::_nanos6Branch = nullptr;
ReportEntry<std::string> *CommandLine::_nanos6CompilerVersion = nullptr;
ReportEntry<std::string> *CommandLine::_nanos6CompilerFlags = nullptr;
ReportEntry<int> *CommandLine::_nanos6CPUs = nullptr;

ReportEntry<std::string> *CommandLine::_nanos5Version = nullptr;
ReportEntry<int> *CommandLine::_nanos5CPUs = nullptr;

ReportEntry<int> *CommandLine::_openmpCPUs = nullptr;

ReportEntry<int> *CommandLine::_nativeCPUs = nullptr;

ReportEntry<int> *CommandLine::_mpiRanks = nullptr;
ReportEntry<int> *CommandLine::_mpiVersionMajor = nullptr;
ReportEntry<int> *CommandLine::_mpiVersionMinor = nullptr;

ReportEntry<std::string> *CommandLine::_sourceVersion = nullptr;
ReportEntry<std::string> *CommandLine::_sourceBranch = nullptr;


void CommandLine::initializeInternal(int argc, char **argv, std::string const &benchmarkName, std::string const &executableType, std::string const &variant)
{
	_argc = argc;
	_argv = argv;
	_currentIndex = 1;
	_totalMandatoryParameters = 0;
	_totalOptionalParameters = 0;
	
	_benchmarkName = new ReportEntry<std::string>("DISCRIMINANT", "benchmark", "Name of the benchmark", benchmarkName);
	_executableType = new ReportEntry<std::string>("DISCRIMINANT", "executable_type", "Type of the executable", executableType);
	_variant = new ReportEntry<std::string>("DISCRIMINANT", "variant", "Algorithm or implementation variation", variant);
// 	_executable = new ReportEntry<std::string>("INFO", "executable", "Name of the executable", argv[0]);
	
	std::ostringstream oss;
	oss << argv[0];
	for (int i=1; i < argc; i++) {
		oss << " " << argv[i];
	}
	
	_commandLine = new ReportEntry<std::string>("INFO", "command_line", "Full command line", oss.str());
}


void CommandLine::validate()
{
	bool failed = false;
	
	if (_totalOptionalParameters == 0) {
		if (_totalMandatoryParameters != _argc-1) {
			std::cerr << "Error: invalid number of parameters. Expected " << _totalMandatoryParameters << " and got " << _argc-1 << " instead." << std::endl;
			failed = true;
		}
	} else {
		if ((_totalMandatoryParameters > _argc-1) || (_totalMandatoryParameters+_totalOptionalParameters < _argc-1)) {
			std::cerr << "Error: invalid number of parameters. Expected from " << _totalMandatoryParameters << " to " << _totalMandatoryParameters+_totalOptionalParameters << " and got " << _argc-1 << " instead." << std::endl;
			failed = true;
		}
	}
	
	for (auto failedIndex : _failedIndexes) {
		std::cerr << "Error: invalid value for parameter " << failedIndex << "." << std::endl;
		failed = true;
	}
	
	if (failed) {
		std::cerr << std::endl;
		std::cerr << "Usage: " << _argv[0];
		for (auto parameter : _mandatoryParameters) {
			assert(parameter != nullptr);
			std::cerr << " [" << parameter->getName() << "]";
		}
		for (auto parameter : _optionalParameters) {
			assert(parameter != nullptr);
			std::cerr << " <" << parameter->getName() << ">";
		}
		std::cerr << std::endl;
		
		if (!_mandatoryParameters.empty()) {
			std::cerr << "\tMandatory parameters:" << std::endl;
		}
		for (auto parameter : _mandatoryParameters) {
			assert(parameter != nullptr);
			std::cerr << "\t\t" << parameter->getName() << ":\t" << parameter->getDescription() << std::endl;
		}
		std::cerr << std::endl;
		
		if (!_optionalParameters.empty()) {
			std::cerr << "\tOptional parameters:" << std::endl;
		}
		for (auto parameter : _optionalParameters) {
			assert(parameter != nullptr);
			std::cerr << "\t\t" << parameter->getName() << ":\t" << parameter->getDescription() << " (Default: " << parameter->getDefaultValueAsString() << ")" << std::endl;
		}
		std::cerr << std::endl;
		
		exit(1);
	}
}
