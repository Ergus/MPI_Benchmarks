/*
 * Copyright (C) 2022  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <unistd.h>
#include <getopt.h>
#include <alloca.h>
#include "mpi.h"

#include <limits.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <type_traits>
#include <tuple>
#include <stack>

#include <slurm/slurm.h>

#include <cstring>

void mpi_error_handler(MPI_Comm *comm, int *err, ... )
{
	const int ret = *err;

	if (ret == MPI_SUCCESS) {
		return;
	}

	int eclass, len;
	char errstring[MPI_MAX_ERROR_STRING];

	MPI_Error_class(ret, &eclass);
	MPI_Error_string(ret, errstring, &len);

	fprintf(stderr, "MPI_ERROR %d: %s\n", eclass, errstring);
	fflush(stderr);
	MPI_Abort(*comm, ret);
}

template<typename T>
bool string_in(const std::string &str, const T &other)
{
	return (str.compare(other) == 0);
}

template<typename... T>
bool string_in(const std::string &str, const std::string &list0, T& ...rest)
{
	return (str.compare(list0) == 0) || string_in(str, rest...);
}

#define types									\
	hlp(TAG_BASE)								\
	hlp(TAG_EXIT)								\
	hlp(TAG_SHRINK)								\
	hlp(TAG_SPAWN)								\
	hlp(TAG_INFO)								\
	hlp(TAG_INFO_REPLY)

// Macro
enum msg_tag {
#define hlp(in) in,
	types
#undef hlp
	NTAGS
};

const static std::string msg_names[] = {
#define hlp(in) "\""#in"\"",
	types
#undef hlp
};

template<msg_tag TYPE, typename...T>
class msg_t {
public:
	msg_tag _type = TYPE;

private:
	using tp = std::tuple<T...>;
	tp _data;

	template <size_t I = 0, typename... Ts>
	typename std::enable_if<I == sizeof...(Ts), void>::type
	printTuple(__attribute__((unused)) std::ostream &out) const
	{
		return;
	}

	template <size_t I = 0, typename... Ts>
	typename std::enable_if<(I < sizeof...(Ts)), void>::type
	printTuple(std::ostream &out) const
	{
		out << this->get<I>() << " ";
		this->printTuple<I + 1, Ts...>(out);
	}

public:
	msg_t(const T& ... rest) : _type(TYPE), _data(rest...)
	{
		static_assert(std::tuple_size<tp>() == sizeof...(rest));
	}

	constexpr size_t getNfields() const
	{
		return std::tuple_size<tp>();
	}

	template<int idx>
	typename std::tuple_element<idx, tp>::type get() const
	{
		static_assert(idx < std::tuple_size<tp>());
		return std::get<idx>(_data);
	}

	friend std::ostream& operator<<(std::ostream& out, const msg_t &s)
	{
		std::stringstream oss;

		oss << msg_names[s._type] << ": ";
		s.printTuple<0, T...>(oss);
		oss << " size(" << sizeof(s) << ")";

		out << oss.str();
		return out;
	}

	int send(int target, MPI_Comm comm)
	{
		int wrank;
		MPI_Comm_rank(comm, &wrank);
		std::cout << "Process:" << wrank << " --(" << *this << ")--> " << target << std::endl;
		return MPI_Send(this, sizeof(*this), MPI_BYTE, target, TYPE, comm);
	}

    void send_to_all(MPI_Comm comm)
	{
		int wrank, wsize;
		MPI_Comm_rank(comm, &wrank);   // gets rank in local world
		MPI_Comm_size(comm, &wsize);   // gets size in local world

		if (wsize == 1) {
			return;
		}

		const int remotes = wsize - 1;
		// Send unblocking spawn messages to all remotes.
		MPI_Request *reqs = (MPI_Request *) malloc(remotes * sizeof(MPI_Request));
		MPI_Status *statuses = (MPI_Status *) malloc(remotes * sizeof(MPI_Status));

		int count = 0;
		for (int i = 0; i < wsize; ++i) {
			if (i == wrank) {
				continue;
			}

			std::cout << "Process:" << wrank << " --(" << *this << ")--> " << i << std::endl;
			MPI_Isend(this, sizeof(*this), MPI_BYTE, i, this->_type, comm, &reqs[count++]);
		}

		assert(count == remotes);

		const int ret = MPI_Waitall(count, reqs, statuses);

		if (MPI_ERR_IN_STATUS == ret) {
			std::cerr << "Error in request status: " << std::endl;
			for (size_t i = 0; i < count; ++i) {
				if (statuses[i].MPI_ERROR == MPI_SUCCESS) {
					continue;
				}
				std::cerr << "Error request[" << i << "]: " << statuses[i].MPI_ERROR
				          << "Sending " << statuses[i].MPI_SOURCE
				          << "--(" << statuses[i].MPI_TAG <<")" << std::endl;
			}
		}

		free(statuses);
		free(reqs);
	}
};

struct hostname_t {
	char name[HOST_NAME_MAX];

	hostname_t()
	{
		name[0] = '\0';
	}

	hostname_t(const char in[])
	{
		strncpy(name, in, HOST_NAME_MAX);
	}

	hostname_t(const std::string &in) : hostname_t(in.c_str())
	{
	}

	hostname_t &operator=(const std::string &in)
	{
		strncpy(name, in.c_str(), HOST_NAME_MAX);
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& out, const hostname_t &s)
	{
		std::string tmp(s.name);
		out << tmp;
		return out;
	}

	operator std::string() const
	{
		return name;
	}
};

typedef msg_t<TAG_BASE> base_t;
typedef msg_t<TAG_SPAWN, hostname_t> spawn_t;
typedef msg_t<TAG_SHRINK, size_t> shrink_t;
typedef msg_t<TAG_INFO> info_t;
typedef msg_t<TAG_INFO_REPLY> info_rep_t;

//================== Node_t ==========================
class Node_t {

protected:
	struct commInfo {
		MPI_Comm intraComm;
		MPI_Comm interComm;
	};

    MPI_Comm _intra, _parent;           // persistent communicators
    int _wsize, _wrank;                 // environment MPI vars
    int _argc;                          // command line arguments number
    char** _argv;                       // command line arguments vars
    bool _listening;                    // process will listen by default
    hostname_t _hostname;               // id of running host
	MPI_Errhandler _neweh;
	std::stack<commInfo> _spawnedCommInfoStack;

    Node_t(int &argc, char** &argv, MPI_Comm parent)
		: _argc(argc), _argv(argv), _parent(parent), _listening(true)
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &_intra);            // Creates a duplicated _intracomm

		MPI_Comm_create_errhandler(mpi_error_handler, &_neweh);
		MPI_Comm_set_errhandler(_intra, _neweh);

		commInfo startcomm = {.intraComm = _intra, .interComm = _parent};
		_spawnedCommInfoStack.push(startcomm);
		gethostname(_hostname.name, HOST_NAME_MAX);

		if (_parent != MPI_COMM_NULL) {
			MPI_Intercomm_merge(_parent, true,  &_intra);
		}

		MPI_Comm_size(_intra, &_wsize);                    // gets size in local world
		MPI_Comm_rank(_intra, &_wrank);                    // gets rank in local world

		std::cout << "Process " << _wrank << ": start" << std::endl;
	}

    int spawn_merge(std::string hostname)       // spawns n new mpi processes
	{
		std::cout << "Process " << _wrank << ": Spawning to: " << hostname << std::endl;

		MPI_Comm new_intra = MPI_COMM_NULL;               // Variable for _intracomm
		MPI_Comm newinter = MPI_COMM_NULL;                // Temporal intercomm
		int errcode;

		MPI_Info info;
		MPI_Info_create(&info);
		MPI_Info_set(info, "host", hostname.c_str());

		const int success = MPI_Comm_spawn(_argv[0], &_argv[1], 1, info,
		                                   0, _intra, &newinter, &errcode);

		if (success == MPI_ERR_SPAWN) {
			std::cerr << "Process " << _wrank << ": Error spawning in " << hostname << std::endl;
			MPI_Abort(_intra, MPI_ERR_OTHER);
		}

		commInfo startcomm = {.intraComm = _intra, .interComm = newinter};
		_spawnedCommInfoStack.push(startcomm);

		MPI_Intercomm_merge(newinter, false, &new_intra); // Create new _intra

		_intra = new_intra;                               // Reassign the _intra to the new one
		MPI_Comm_size(new_intra, &_wsize);                // update _wsize

		MPI_Info_free(&info);

		std::cout << "Process " << _wrank << ": Spawned done" << std::endl;

		return MPI_SUCCESS;
	}


	///! split the communicator and kills
    int shrink_pop(size_t n)
	{
		std::cout << "Process " << _wrank <<": Reducing (world: " << _wsize << ")" << "\n";

		size_t removed = 0;
		const int newsize = _wsize - n;
		const bool willcontinue = _wrank < newsize;

		int ret, currsize = _wsize;

		while (currsize > newsize && _spawnedCommInfoStack.size() > 0) {
			struct commInfo &spawnedCommInfo = _spawnedCommInfoStack.top();

			MPI_Comm_free(&_intra);                           // Free old _intracomm before.
			MPI_Comm_disconnect(&spawnedCommInfo.interComm);  // disconnect the inter communicator
			_intra = spawnedCommInfo.intraComm;               // Reassign the _intra to the new one

			_spawnedCommInfoStack.pop();
			--currsize;
		}

		if (!willcontinue) {
			assert(_spawnedCommInfoStack.size() == 0);
			assert(_wrank > 0);
		}

		MPI_Comm_size(_intra, &_wsize);                       // update _wsize
		return willcontinue;
	}

	void print_info()
	{
		std::cout << "Process: " << _wrank << "/" << _wsize << " host: " << _hostname << std::endl;
	}

public:
    virtual ~Node_t()
	{
		MPI_Comm_free(&_intra);                          // Free the _intra comm
		std::cout << "Process " << _wrank << ": Exit" << std::endl;
	}

    virtual void run() = 0;
};

//================== Node_master ==========================

class Node_master: public Node_t {

private:

	const std::vector<hostname_t> _hostList;

	static std::vector<hostname_t> getHostList()
	{
		const std::string nodelist = getenv("SLURM_NODELIST");
		std::vector<hostname_t> nodelist_vector;

		hostlist_t hostlist = slurm_hostlist_create(nodelist.c_str());
		if (hostlist == NULL) {
			printf("slurm_hostlist_create returned NULL error %s\n",
			       slurm_strerror(slurm_get_errno()));
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
		}

		// Initialize the node list.
		char *host;
		while ((host = slurm_hostlist_shift(hostlist))) {
			nodelist_vector.push_back(host);
		}
		slurm_hostlist_destroy(hostlist);

		if (nodelist_vector.empty()) {
			printf("Error nodelist_vector is empty.");
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
		}

		return nodelist_vector;
	}

    int master_info()
	{
		info_t msg;
		info_rep_t reply;
		for (int i = 0; i < _wsize; ++i) {
			if (i == _wrank) {
				print_info();
			} else {
				msg.send(i, _intra);

				MPI_Recv(&reply, sizeof(info_rep_t),
				         MPI_BYTE, i, TAG_INFO_REPLY, _intra, MPI_STATUS_IGNORE);
			}
		}
	}

    void process(char opt, size_t n)
	{
		std::cout<< "Master processing: " << opt << " " << n << std::endl;
		switch (opt) {
		case 's':
			assert(n < _hostList.size());
			spawn_t msg_spawn(_hostList[n]);
			msg_spawn.send_to_all(_intra);
			spawn_merge(_hostList[n]);
			break;
		case 'd':
			assert(n < _wsize);
			shrink_t msg_shrink(n);
			msg_shrink.send_to_all(_intra);
			shrink_pop(n);
			break;
		case 'i':
			printf("This will print actual status information\n");
			master_info();
			break;
		case 'e':
			msg_t<TAG_EXIT> msg_stop;
			msg_stop.send_to_all(_intra);
			_listening = false;
			break;
		case '?':
			printf("Option '%c' not recognized\n", opt);
			MPI_Abort(_intra, MPI_ERR_OTHER);
		}
	}

public:
    Node_master(int &argc, char** &argv, MPI_Comm parent)
		:Node_t(argc, argv, parent), _hostList(Node_master::getHostList())
	{
		assert(_parent == MPI_COMM_NULL);
		assert(_wrank == 0);
		assert(_wsize == 1);
		std::cout << "Hostlist: ";
		for (auto const &host : _hostList) {
			std::cout << host << " ";
		}
		std::cout << std::endl;
	}

    ~Node_master()
	{
		printf("Deleting master\n");
	}

    void run() override
	{
		printf("Process Master ready\n");
		size_t value = 0;

		if (_argc > 1) {
			int opt;
			while ((opt = getopt(_argc, _argv, "s:d:pih")) != -1) {
				value = 0;
				if (strchr("sd", opt) != NULL) {
					sscanf(optarg, "%zu", &value);
				} else if (opt == '?') {
					std::cerr << "Invalid command line or argument: " << opt << std::endl;
					continue;
				}

				process(opt, value);
			}
		} else {
			std::string command;
			while (_listening) {
				std::cin >> command;
				value = 0;

				if (string_in(command, "spawn", "delete")) {
					std::cin >> value;
				} else if (!string_in(command, "info", "exit")) {
					std::cerr << "Invalid command or argument: " << command << std::endl;
					continue;
				}

				process(command[0], value);
			}
		}
	}
};

//================== Node_slave ==========================

class Node_slave: public Node_t {

private:
    void run()
	{
		int count;
		MPI_Status status;

		printf("Process %d listening\n", _wrank);

		while (_listening) {
			// Prove the message first
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, _intra, &status);

			// Determine size in bytes
			MPI_Get_count(&status, MPI_BYTE, &count);
			assert(count!=0);

			void *buff = alloca(count);

			// Now receive the message
			MPI_Recv(buff, count, MPI_BYTE,
			         status.MPI_SOURCE, status.MPI_TAG, _intra, MPI_STATUS_IGNORE);

			base_t *msg = reinterpret_cast<base_t *>(buff);
			const msg_tag type = msg->_type;

			std::cout << "Process:" << _wrank << " <--(" << *msg << ")-- " << status.MPI_SOURCE << std::endl;

			switch (type) {
			case TAG_SPAWN:
				const spawn_t *msg_spawn = reinterpret_cast<spawn_t *>(buff);
				spawn_merge(msg_spawn->get<0>());
				break;
			case TAG_SHRINK:
				const shrink_t *msg_shrink = reinterpret_cast<shrink_t *>(buff);
				_listening = shrink_pop(msg_shrink->get<0>());
				break;
			case TAG_EXIT:
				_listening = false;
				break;
			case TAG_INFO:
				print_info();
				info_rep_t reply;
				reply.send(status.MPI_SOURCE, _intra);
				break;
			default:
				std::cout << "Process: " << _wrank
				          << " received unknown message type: " << type
				          << std::endl;
			}
		}
		printf("Process %d: Exit listening (world %d)\n", _wrank, _wsize);
	}

public:
    Node_slave(int &argc, char** &argv, MPI_Comm parent)
		: Node_t(argc, argv, parent)
	{
	}

    ~Node_slave() {
		printf("Killing slave %d.\n", _wrank);
	}
};

int main(int argc, char** argv)
{
    MPI_Comm parent;
    int local_wsize, local_wrank;

	int provided;
	MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	if (provided != MPI_THREAD_MULTIPLE) {
		printf("Error thread provided=%d expected=%d\n", provided, MPI_THREAD_MULTIPLE);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &local_wrank); // local world rank
	MPI_Comm_size(MPI_COMM_WORLD, &local_wsize); // local world size
	MPI_Comm_get_parent(&parent);                // Parent to decide

	Node_t *_singleton;
	if (parent == MPI_COMM_NULL && local_wrank == 0) {
		_singleton = new Node_master(argc, argv, parent);
	} else {
		_singleton = new Node_slave(argc, argv, parent);
	}

	_singleton->run();

	delete _singleton;
	_singleton = nullptr;
	MPI_Finalize();
	printf("Exiting %d in world %d\n",local_wrank,local_wsize);

	return 0;
}

