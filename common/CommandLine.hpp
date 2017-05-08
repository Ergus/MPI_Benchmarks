#ifndef COMMAND_LINE_HPP
#define COMMAND_LINE_HPP

#include "Reportable.hpp"

#include <cstdlib>
#include <iostream>
#include <list>
#include <string>

#ifdef _NANOS_H_
#include <cstdio>
#include <sstream>
#include <string>

#include <omp.h>
#elif defined(_OPENMP)
#include <omp.h>
#endif


template <typename T>
class CommandLineParameter;

template <typename T>
class OptionalCommandLineParameter;

template <typename T>
class ReportEntry;


class CommandLine {
protected:
	static int _argc;
	static char **_argv;
	static int _currentIndex;
	static std::list<int> _failedIndexes;
	static int _totalMandatoryParameters;
	static int _totalOptionalParameters;
	static std::list<Reportable *> _mandatoryParameters;
	static std::list<ReportableWithDefaultValue *> _optionalParameters;
	static ReportEntry<std::string> *_benchmarkName;
	static ReportEntry<std::string> *_executableType;
	static ReportEntry<std::string> *_variant;
// 	static ReportEntry<std::string> *_executable;
	static ReportEntry<std::string> *_commandLine;
	
	static ReportEntry<std::string> *_mcxxVersion;
	static ReportEntry<std::string> *_compilerFlags;
	static ReportEntry<std::string> *_compilerVersion;
	
	static ReportEntry<std::string> *_nanos6Version;
	static ReportEntry<std::string> *_nanos6Branch;
	static ReportEntry<std::string> *_nanos6CompilerVersion;
	static ReportEntry<std::string> *_nanos6CompilerFlags;
	static ReportEntry<int> *_nanos6CPUs;
	
	static ReportEntry<std::string> *_nanos5Version;
	static ReportEntry<int> *_nanos5CPUs;
	
	static ReportEntry<int> *_openmpCPUs;
	
	static ReportEntry<int> *_nativeCPUs;
	
	static ReportEntry<int> *_mpiRanks;
	static ReportEntry<int> *_mpiVersionMajor;
	static ReportEntry<int> *_mpiVersionMinor;
	
	static ReportEntry<std::string> *_sourceVersion;
	static ReportEntry<std::string> *_sourceBranch;
	
	
	static void initializeInternal(int argc, char **argv, std::string const &benchmarkName, std::string const &executableType, std::string const &variant = std::string("default"));
	
public:
	CommandLine()
	{
	}
	
	#ifdef EXECUTABLE_TYPE
        static inline void initialize(int argc, char **argv, std::string const &benchmarkName, std::string const &variant = std::string("default"));
	#endif
	
	template <typename T>
	static inline void registerMandatoryParameter(CommandLineParameter<T> *parameter)
	{
		if (_totalOptionalParameters != 0) {
			std::cerr << "Internal program error: optional parameters must come after all mandatory parameters." << std::endl;
			exit(1);
		}
		
		_mandatoryParameters.push_back(parameter);
		if (parameter->bad()) {
			_failedIndexes.push_back(_totalMandatoryParameters);
		}
		_totalMandatoryParameters++;
	}
	
	template <typename T>
	static inline void registerOptionalParameter(OptionalCommandLineParameter<T> *parameter)
	{
		_optionalParameters.push_back(parameter);
		if (parameter->bad()) {
			_failedIndexes.push_back(_totalMandatoryParameters + _totalOptionalParameters);
		}
		_totalOptionalParameters++;
	}
	
	static void validate();
	
	template <typename T>
	friend class CommandLineParameter;
	
	template <typename T>
	friend class OptionalCommandLineParameter;
	
};


#ifdef EXECUTABLE_TYPE
	#include "ReportEntry.hpp"
	
	#include "VersionInfo.hpp"
	
	#ifdef __NANOS6__
	#include <nanos6/debug.h>
	#endif
	
	#ifdef USE_MPI
	#include <mpi.h>
	#endif
	
	void CommandLine::initialize(int argc, char **argv, std::string const &benchmarkName, std::string const &variant)
	{
		initializeInternal(argc, argv, benchmarkName, EXECUTABLE_TYPE, variant);
		
		_mcxxVersion = new ReportEntry<std::string>("SOFT", "mercurium_version", "Mercurium version", CXX_VERSION);
		_compilerFlags = new ReportEntry<std::string>("SOFT", "compiler_flags", "Compiler flags", CXXFLAGS);
		_compilerVersion = new ReportEntry<std::string>("SOFT", "native_compiler_version", "Native compiler version", NATIVE_CXX_VERSION);
		
		#ifdef __NANOS6__
		_nanos6Version = new ReportEntry<std::string>("SOFT", "nanos6_version", "Nanos6 version", nanos_get_runtime_version());
		_nanos6Branch = new ReportEntry<std::string>("SOFT", "nanos6_branch", "Nanos6 branch", nanos_get_runtime_branch());
		_nanos6CompilerVersion = new ReportEntry<std::string>("SOFT", "nanos6_compiler_version", "Nanos6 compiler version", nanos_get_runtime_compiler_version());
		_nanos6CompilerFlags = new ReportEntry<std::string>("SOFT", "nanos6_compiler_flags", "Nanos6 compiler flags", nanos_get_runtime_compiler_flags());
		_nanos6CPUs = new ReportEntry<int>("PARAMETER", "num_cpus", "The number of CPUs used", nanos_get_num_cpus());
		#endif
		
		#ifdef _NANOS_H_
		{
			std::ostringstream oss;
			oss
				<< "eval $("
					<< "ldd " << argv[0]
					<< " | grep libnanox.so.1"
					<< " | awk '{ print $3; }'"
					<< " | sed 's@/lib/[^/]*/libnanox[.]so.*@/bin/nanox --version@'"
				<< ")"
				<< " | head -1"
				<< " | sed 's/^[^(]*(//; s/[)]$//'";
			FILE *versionPipe = popen(oss.str().c_str(), "r");
			char buff[4096];
			
			size_t size = fread(buff, 1, 4095, versionPipe);
			while ((size > 0) && (buff[size-1] == '\n')) {
				size--;
			}
			if (size > 1) {
				buff[size] = '\0';
				std::string version(buff);
				_nanos5Version = new ReportEntry<std::string>("SOFT", "nanos5_version", "Nanos5 version", buff);
			}
			
			_nanos5CPUs = new ReportEntry<int> ("PARAMETER", "num_cpus", "The number of CPUs used", omp_get_max_threads());
		}
		#elif defined(_OPENMP)
			_openmpCPUs = new ReportEntry<int> ("PARAMETER", "num_cpus", "The number of CPUs used", omp_get_max_threads());
		#elif defined(__NANOS6__)
		#else
			_nativeCPUs = new ReportEntry<int> ("PARAMETER", "num_cpus", "The number of CPUs used", 1);
		#endif
		
		#ifdef USE_MPI
			{
				int nranks = 0;
				MPI_Comm_size(MPI_COMM_WORLD, &nranks);
				_mpiRanks = new ReportEntry<int> ("PARAMETER", "mpi_processes", "The number of MPI processes used", nranks);
				
				int major = 0;
				int minor = 0;
				MPI_Get_version(&major, &minor);
				_mpiVersionMajor = new ReportEntry<int> ("SOFT", "mpi_major", "The major version of MPI", major);
				_mpiVersionMinor = new ReportEntry<int> ("SOFT", "mpi_major", "The minor version of MPI", minor);
			}
		#endif
		
		#ifdef SOURCE_VERSION
		_sourceVersion = new ReportEntry<std::string>("INFO", "benchmark_version", "Benchmark version", SOURCE_VERSION);
		_sourceBranch = new ReportEntry<std::string>("INFO", "benchmark_branch", "Benchmark branch", SOURCE_VERSION);
		#endif
	}
#endif


#endif // COMMAND_LINE_HPP
