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

#include "benchmarks_mpi.h"

int main(int argc, char *argv[])
{
	init_args(argc, argv);

	// The distribution will be: associated with the number of nodes.
	const int fullsize = create_cl_int("Size");

	envinfo env;

	Initialize(&env, &argc, &argv, fullsize, 1);

	double *array = (double *) malloc(fullsize * sizeof(double));

	size_t counter = 0; // number of elements I own
	for (size_t i = env.rank; i < fullsize; i += env.worldsize) {
		array[i] = (double) env.rank;
		++counter;
	}

	timer ttimer = create_timer("Total time");

	const size_t nsends = (env.worldsize - 1) * counter;
	const size_t nrecv = fullsize - counter;
	int total_requests = nsends + nrecv;

	MPI_Request *reqs = (MPI_Request *) malloc(total_requests * sizeof(MPI_Request));
	MPI_Status *stats = (MPI_Status *) malloc(total_requests * sizeof(MPI_Status));

	int req_idx = 0;

	// Send and receive
	for (size_t i = 0; i < fullsize; ++i) {

		int owner = i % env.worldsize;

		if (owner == env.rank) { // Send this element to the (w-1) remotes.

			for (int remote = 0; remote < env.worldsize; ++remote) {
				if (remote != env.rank) {
					MPI_Isend(&array[i], 1, MPI_DOUBLE, remote, i, MPI_COMM_WORLD, &reqs[req_idx++]);
				}
			}

		} else { // Wait for the element from it's owner.
			MPI_Irecv(&array[i], 1, MPI_DOUBLE, owner, i, MPI_COMM_WORLD, &reqs[req_idx++]);
		}
	}


	MPI_Waitall(total_requests, reqs, stats);
	stop_timer(&ttimer);

	free(reqs);
	free(stats);
	free(array);

	report_args();
	free_args();

	Finalize();

    return 0;
}
