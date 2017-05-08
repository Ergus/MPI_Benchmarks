#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/utsname.h>

#include <stdio.h>
#include <sys/time.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Report.hpp"
#include "ReportEntry.hpp"


Report Report::_singleton;


static ReportEntry<std::string> *_utsSysname;
static ReportEntry<std::string> *_utsNodename;
static ReportEntry<std::string> *_utsRelease;
static ReportEntry<std::string> *_utsVersion;
static ReportEntry<std::string> *_utsMachine;
static ReportEntry<std::string> *_utsDomainname;
static ReportEntry<long> *_startTime;

Report::Report()
	: _reportEntries()
{
	struct utsname uts;
	
	int rc = uname(&uts);
	if (rc != 0) {
		perror("uname");
		exit(1);
	}
	
	_utsSysname = new ReportEntry<std::string>("HOST", "os", "Operating system name", uts.sysname);
	_utsNodename = new ReportEntry<std::string>("HOST", "hostname", "Host name within network", uts.nodename);
	_utsRelease = new ReportEntry<std::string>("HOST", "os_release", "Operating system release", uts.release);
	_utsVersion = new ReportEntry<std::string>("HOST", "os_version", "Operating system version", uts.version);
	_utsMachine = new ReportEntry<std::string>("HOST", "host_id", "Host hardware identifier", uts.machine);
	_utsDomainname = new ReportEntry<std::string>("HOST", "domainname", "Host domain name", uts.nodename);
	
	struct timeval tv;
	
	rc = gettimeofday(&tv, nullptr);
	if (rc != 0) {
		perror("gettimeofday");
		exit(1);
	}
	
	_startTime = new ReportEntry<long>("TIME", "start_time", "Starting time", tv.tv_sec, "seconds");
}

void Report::emit(std::ostream &o)
{
	for (auto entry : _singleton._reportEntries) {
		assert(entry != nullptr);
		o << "REPORT\t";
		entry->emitReport(o);
	}
}
