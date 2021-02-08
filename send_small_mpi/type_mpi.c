/*
 * Copyright (C) 2021  Jimmy Aguilar Mena
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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>

#include "mpi.h"

#include "ArgParserC/argparser.h"

#include "benchmarks_mpi.h"
#include "extrae_user_events.h"

int main(int argc, char *argv[])
{
	init_args(argc, argv);

	// The distribution will be: associated with the number of nodes.
	const int fullsize = create_cl_int("Size");

	envinfo env;

	Initialize(&env, &argc, &argv, fullsize, 1);

	// Add some instrumentation for extrae.
	extrae_type_t type = 50000007;
	int nvalues = 2;
	extrae_value_t values[2] = {0, 1};
	char * description_values[2] = {"End", "MPI_Type_indexed"};

	Extrae_define_event_type (&type, "MPI_Type_indexed", &nvalues, values, description_values);

	// Start the code.
	double *array = (double *) malloc(fullsize * sizeof(double));

	size_t counter = 0; // number of elements I own
	for (size_t i = env.rank; i < fullsize; i += env.worldsize) {
		array[i] = (double) env.rank;
		++counter;
	}

	timer ttimer = create_timer("Total time");

	const size_t nsends = (env.worldsize - 1);
	const size_t nrecv = (env.worldsize - 1);
	int total_requests = nsends + nrecv;

	MPI_Request *reqs = (MPI_Request *) malloc(total_requests * sizeof(MPI_Request));
	MPI_Status *stats = (MPI_Status *) malloc(total_requests * sizeof(MPI_Status));

	int req_idx = 0;

	// Allocated auxiliars.
	int *counts = calloc(env.worldsize, sizeof(int));

	int **lengths = (int **) alloca(env.worldsize * sizeof(int *));
	int **displs = (int **) alloca(env.worldsize * sizeof(int *));

	size_t frac = fullsize / env.worldsize;
	size_t rest = fullsize % env.worldsize;

	for (size_t i = 0; i < env.worldsize; ++i) {
		lengths[i] = malloc( (frac + (i < rest)) * sizeof (int) );
		displs[i] = malloc( (frac + (i < rest)) * sizeof (int) );
	}

	// Initialize auxiliars
	for (size_t i = 0; i < fullsize; ++i) {
		size_t owner = i % env.worldsize;
		size_t idx = counts[owner]++;

		lengths[owner][idx] = 1;
		displs[owner][idx] = i;
	}

	// Send and receive
	for (size_t i = 0; i < env.worldsize; ++i) {
		Extrae_event(type, 1);

		MPI_Datatype sendtype;
		MPI_Type_indexed(3, lengths[env.rank], displs[env.rank], MPI_INT, &sendtype);
		MPI_Type_commit(&sendtype);

		Extrae_event(type, 0);

		if (i == env.rank) {  // Send to all
			for (size_t j = 0; j < env.worldsize; ++j) {
				if (j != env.rank) {
					MPI_Isend(array, 1, sendtype, j, i, MPI_COMM_WORLD, &reqs[req_idx++]);
				}
			}
		} else {  // Recv from i
			MPI_Irecv(array, 1, sendtype, i, i, MPI_COMM_WORLD, &reqs[req_idx++]);
		}

		// The type is deleted when all the pending operations finish.
		MPI_Type_free(&sendtype);
	}

	MPI_Waitall(total_requests, reqs, stats);

	stop_timer(&ttimer);

	for (size_t i = 0; i < env.worldsize; ++i) {
		free(lengths[i]);
		free(displs[i]);
	}

	free(reqs);
	free(stats);
	free(array);

	report_args();
	free_args();

	Finalize();

    return 0;
}
