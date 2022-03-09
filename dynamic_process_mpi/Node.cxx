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

#include <sstream>
#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <type_traits>
#include <tuple>

#define check_mpi(ret) {										\
		if (ret != MPI_SUCCESS) {								\
			std::cerr << "Error in MPI call %d\n" << __LINE__;	\
			MPI_Abort(intra, ret);								\
		}														\
	}


void check_mpi_status(int ret, MPI_Status *statuses, size_t total) {
	if (MPI_ERR_IN_STATUS == ret) {
		std::cerr << "Error in request: " << std::endl;
		for (size_t i = 0; i < total; ++i) {
			if (statuses[i].MPI_ERROR == MPI_SUCCESS) {
				continue;
			}
			std::cerr << "Error request[" << i << "]: " << statuses[i].MPI_ERROR
			          << "Sending " << statuses[i].MPI_SOURCE
			          << "--(" << statuses[i].MPI_TAG <<")" << std::endl;
		}
	}
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
	hlp(TAG_EXIT)								\
	hlp(TAG_REDUCE)								\
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


template<typename...T>
class msg_t {
public:
	msg_tag _type;

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
	msg_t(msg_tag type, const T& ... rest)
		: _type(type), _data(rest...)
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
		std::string msg_names[] = {
			#define hlp(in) "\""#in"\"",
			types
			#undef hlp
		};

		std::stringstream oss;
 
		oss << "type: [" << msg_names[s._type] << "]: ";
		s.printTuple<0, T...>(oss);
		oss << " size(" << sizeof(s) << ")";

		out << oss.str();
		return out;
	}
};


typedef msg_t<int , long> info_t;

//================== Node_t ==========================
class Node_t {
protected:
    MPI_Comm intra, parent;            // persistent communicators
    int wsize, wrank;                  // environment MPI vars
    int nargc;                         // command line arguments number
    char** nargv;                      // command line arguments vars
    bool listening;                    // process will listen by default
    long hostid;                       // id of running host

    Node_t(int &argc, char** &argv, MPI_Comm _parent)
		: nargc(argc), nargv(argv), parent(_parent),
		  listening(true), hostid(gethostid())
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &intra);            // Creates a duplicated intracomm
		MPI_Comm_size(intra, &wsize);                    // gets size in local world
		MPI_Comm_rank(intra, &wrank);                    // gets rank in local world
		std::cout << "Process " << wrank << ": start" << std::endl;
	}

	template<typename...T>
    int send_to_all_remotes(msg_t<T...> &msg)
	{
		// Only send to remotes if there are any
		if (wsize <= 1) {
			return MPI_SUCCESS;
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

			std::cout << "Sending msg:(" << msg << ") to:" << i << std::endl;
			MPI_Isend(&msg, sizeof(msg), MPI_BYTE, i, msg._type, intra, &reqs[count++]);
		}

		const int ret = MPI_Waitall(count, reqs, statuses);
		check_mpi_status(ret, statuses, count);

		free(statuses);
		free(reqs);
		return MPI_SUCCESS;
	}

    int spawn_merge(size_t n)                    // spawns n new mpi processes
	{
		std::cout << "Process "<< wrank << ": Spawning in world " << wsize << std::endl;
		for (size_t i = 0; i < n; ++i) {

			MPI_Comm newintra = MPI_COMM_NULL;               // Variable for intracomm
			MPI_Comm newinter = MPI_COMM_NULL;               // Temporal intercomm
			int errcode;

			int success = MPI_Comm_spawn(nargv[0], &nargv[1], 1, MPI_INFO_NULL,
			                             0, intra, &newinter, &errcode);

			if (success == MPI_ERR_SPAWN) {
				std::cerr << "Process " << wrank << ": Error spawning " << i << std::endl;
			}

			MPI_Comm_free(&intra);                           // Free old intracomm before.
			MPI_Intercomm_merge(newinter, false, &newintra); // Create new intra

			intra = newintra;                                // Reassign the intra to the new one
			MPI_Comm_size(newintra, &wsize);                 // update wsize
			MPI_Comm_free(&newinter);                        // Free the created intercomm
		}
		std::cout << "Process " << wrank << ": Spawned done" << std::endl;

		return MPI_SUCCESS;
	}


	///! split the communicator and kills
    int split_kill(size_t n)
	{
		MPI_Comm newintra = MPI_COMM_NULL;  // Variable for intracomm
		listening = (wrank < (wsize - n));

		printf("Process %d reducing %zu processes (listening = %d world %d)\n",
		        wrank, n, listening, wsize);

		MPI_Comm_split(intra, (int) listening, wrank, &newintra);

		MPI_Comm_free(&intra);              // Free old intracomm before.
		intra = newintra;                   // Reassign the intra to the new one

		MPI_Comm_size(newintra, &wsize);    // update wsize
		printf("Process %d reduced (world %d) \n",wrank,wsize);
		return MPI_SUCCESS;
	}

public:
    virtual ~Node_t()
	{
		MPI_Comm_free(&intra);                          // Free the intra comm
		std::cout << "Process " << wrank << ": Exit" << std::endl;
	}

    virtual void run() = 0;
};

//================== Node_master ==========================

class Node_master: public Node_t {
private:
    int master_getinfo()
	{
		msg_t<> msg_info(TAG_INFO);

		info_t *infos = (info_t *) malloc(wsize * sizeof(info_t));
		MPI_Request *reqs = (MPI_Request *) malloc(wsize * sizeof(MPI_Request));

		int count = 0;
		for (int i = 0; i < wsize; ++i) {
			if (i == wrank) {
				infos[i] = info_t(TAG_INFO_REPLY, wrank , hostid);
				continue;
			}
			MPI_Request reqsend;
			MPI_Isend(&msg_info, sizeof(msg_info), MPI_BYTE, i, msg_info._type, intra, &reqsend);
			MPI_Irecv(&infos[i], sizeof(info_t), MPI_BYTE, i, TAG_INFO_REPLY, intra, &reqs[count++]);
			MPI_Request_free(&reqsend);
		}

		MPI_Status *status = (MPI_Status *) malloc(wsize * sizeof(MPI_Status));

		const int ret = MPI_Waitall(count, reqs, status);
		check_mpi_status(ret, status, count);

		std::cout << "Info: world = " << wsize << std::endl;
		for (int i = 0; i < wsize; ++i) {
			std::cout << infos[i] << std::endl;
		}

		free(status);
		free(reqs);
		free(infos);
	}

    void process(char opt, int n=0)
	{
		switch (opt) {
		case 's':
			msg_t<int> msg_spawn(TAG_SPAWN, n);
			send_to_all_remotes(msg_spawn);

			spawn_merge(n);
			break;
		case 'd':
			msg_t<int> msg_red(TAG_REDUCE, n);
			send_to_all_remotes(msg_red);

			split_kill(n);
			break;
		case 'p':
			fflush(stderr);
			printf("Press enter to continue\n");
			getchar();
			break;
		case 'i':
			printf("This will print actual status information\n");
			master_getinfo();
			break;
		case 'h':
			printf("\tCommand line: %s -[hip] -[sd] [value]\n", nargv[0]);
			printf("\tInteractive: Use tab for available commands\n");
			break;
		case 'e':
			msg_t<> msg_stop(TAG_EXIT);
			send_to_all_remotes(msg_stop);

			listening = false;
			break;
		case '?':
			printf("Option %c not recognised\n", opt);
			MPI_Abort(intra, MPI_ERR_OTHER);
		}
	}

public:
    Node_master(int &argc, char** &argv, MPI_Comm _parent)
		:Node_t(argc, argv, _parent)
	{
		assert(wrank == 0);
	}

    ~Node_master()
	{
		printf("Deleting master\n");
	}

    void run() override
	{
		printf("Starting master\n");

		if (nargc > 1) {
			int opt, value;
			while ((opt = getopt(nargc, nargv, "s:d:pih")) != -1) {
				value = (opt == 's' || opt == 'd') ? atoi(optarg) : 0;
				process(opt, value);
			}
		} else {
			while (listening) {
				std::string command;
				std::cin >> command;

				if (string_in(command, "spawn", "delete")) {
					int value = 0;
					std::cin >> value;
					process(command[0], value);
				} else if (string_in(command, "info", "exit", "help")) {
					process(command[0], 0);
				} else {
					std::cerr << "Invalid command or argument: " << command << std::endl;
				}
			}
		}
	}
};

//================== Node_slave ==========================

class Node_slave: public Node_t {

private:
    void listen()
	{
		int ret, count;
		MPI_Status status;

		printf("Process %d listening\n", wrank);

		while (listening) {

			// Prove the message first
			ret = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, intra, &status);
			check_mpi(ret);

			// Determine size in bytes
			ret = MPI_Get_count(&status, MPI_BYTE, &count);
			check_mpi(ret);
			assert(count!=0);

			void *buff = alloca(count);

			// Now receive the message
			ret = MPI_Recv(buff, count, MPI_BYTE,
			               status.MPI_SOURCE, status.MPI_TAG, intra,
			               MPI_STATUS_IGNORE);
			check_mpi(ret);

			msg_t<> *msg = reinterpret_cast<msg_t<> *>(buff);
			const msg_tag type = msg->_type;

			std::cout << wrank << "<--("<< *msg << ")--" << status.MPI_SOURCE << "\n";

			switch (type) {
			case TAG_SPAWN:
				std::cout << "Process " << wrank <<": Spawning (world: " << wsize << ")" << "\n";
				const msg_t<int> *msg_int1 = reinterpret_cast<msg_t<int> *>(buff);
				spawn_merge(msg_int1->get<0>());
				break;
			case TAG_REDUCE:
				std::cout << "Process " << wrank <<": Reducing (world: " << wsize << ")" << "\n";
				const msg_t<int> *msg_int2 = reinterpret_cast<msg_t<int> *>(buff);
				split_kill(msg_int2->get<0>());
				break;
			case TAG_EXIT:
				listening = false;
				break;
			case TAG_INFO:
				info_t info(TAG_INFO_REPLY, wrank , hostid);
				MPI_Send(&info, sizeof(info), MPI_BYTE, status.MPI_SOURCE, TAG_INFO_REPLY, intra);
				std::cout << wrank << "--("<< info << ")-->" << status.MPI_SOURCE << "\n";
				break;
			default:
				std::cout << "Process: " << wrank
				          << " Received unknown message type: " << type
				          << std::endl;
			}
		}
		printf("Process %d: Exit listening (world %d)\n", wrank, wsize);
	}

public:
    Node_slave(int &argc, char** &argv, MPI_Comm _parent)
		: Node_t(argc, argv, _parent)
	{
		MPI_Intercomm_merge(parent, true,  &intra);

		MPI_Comm_size(intra, &wsize);
		MPI_Comm_rank(intra, &wrank);
	}

    ~Node_slave() {
		printf("Killing all slaves.\n");
	}

    void run() override
	{
		// Create a thread to listen incoming messages.
		std::thread tlisten(&Node_slave::listen, this);

		//Extra code goes here. Always before the join.
		tlisten.join();
	}
};

//================== Manager ==========================

class Manager {
private:
    Node_t *Node;
    MPI_Comm parent;
    int local_wsize, local_wrank;

public:
    Manager(int &argc, char** &argv) : Node(nullptr)
	{
		int provided;
		int ret = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE || ret!=0) {
			printf("Error in MPI initialisation\n");
			printf("Thread support=%d expected=%d returned=%d\n",
			       provided,MPI_THREAD_MULTIPLE, ret);
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
		}

		MPI_Comm_rank(MPI_COMM_WORLD, &local_wrank); // local world rank
		MPI_Comm_size(MPI_COMM_WORLD, &local_wsize); // local world size
		MPI_Comm_get_parent(&parent);                // Parent to decide

		if (parent == MPI_COMM_NULL && local_wrank == 0) {
			Node = new Node_master(argc, argv, parent);
		} else {
			Node = new Node_slave(argc, argv, parent);
		}

		assert(Node != nullptr);
	}

    ~Manager()
	{
		assert(Node != nullptr);
		delete Node;
		MPI_Finalize();
		printf("Exiting %d in world %d\n",local_wrank,local_wsize);
	}

    int run()
	{
		Node->run();
	}

};


int main(int argc, char** argv)
{
	Manager man(argc,argv);
	man.run();

	return 0;
}

