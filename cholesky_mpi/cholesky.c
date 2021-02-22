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

#include "util.h"
#include "benchmarks_mpi.h"
#include "ArgParserC/argparser.h"


void cholesky(size_t nblocks, size_t bsize, size_t rank, size_t wsize,
              double A[nblocks / wsize][nblocks][bsize][bsize])
{
	modcheck(nblocks, wsize);

	const size_t bsize2 = bsize * bsize;

	const size_t lblocks = nblocks / wsize;	    // local blocks
	const size_t first_block = lblocks * rank;	// first block for this rank

	const size_t doubles_region = lblocks * bsize2;
	const size_t nrequests = (wsize - rank - 1);

	double (*buffer)[bsize][bsize] = NULL;
	double (*buffer0)[lblocks][bsize][bsize] = NULL;
	MPI_Request requests[nrequests];

	// Then receive the portions in ranks > 0
	if (rank > 0) {
		buffer0 = malloc(2 * lblocks * bsize2 * sizeof(double));	//buffer[i:nblocks]
		buffer = malloc((wsize - rank) * lblocks * bsize2 * sizeof(double));	//buffer[i:nblocks]
		assert(buffer0 && buffer);

		MPI_Irecv((void *)buffer0[0], doubles_region, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &requests[0]);

		for (size_t region = rank + 1; region < wsize; ++region) {
			const size_t local_first_block = region * lblocks - first_block;

			MPI_Irecv(
				(void *)buffer[local_first_block], doubles_region, MPI_DOUBLE,
				0, region, MPI_COMM_WORLD,
				&requests[region - rank]
			);
		}
	}

	// Wait data until my turn
	// Process rows from previous processes until iteration gets here
	for (size_t Ri = 0; Ri < first_block; ++Ri) {
		const size_t next_Ri = Ri + 1;
		const size_t buff0index = Ri % 2;

		const int source = next_Ri / lblocks;
		const int tag = next_Ri << 16;

		MPI_Wait(&requests[0], MPI_STATUS_IGNORE);

		if (next_Ri < first_block) {
			const size_t buff0nextindex = next_Ri % 2;
			MPI_Irecv((void *)buffer0[buff0nextindex],
			          doubles_region, MPI_DOUBLE, source, tag | rank,
			          MPI_COMM_WORLD, &requests[0]);
		}

		for (size_t j = 0; j < lblocks; ++j) {	// local triangle
			for (size_t k = 0; k < j; ++k) {
				oss_gemm(bsize, buffer0[buff0index][j],
				         buffer0[buff0index][k], A[k][first_block + j]);
			}

			oss_syrk(bsize, buffer0[buff0index][j], A[j][first_block + j]);
		}

		for (size_t region = rank + 1; region < wsize; ++region) {
			const int rtag = tag | region;
			const size_t localindex = region - rank;
			const size_t local_first_block = region * lblocks - first_block;

			MPI_Wait(&requests[localindex], MPI_STATUS_IGNORE);

			for (size_t j = 0; j < lblocks; ++j) {
				for (size_t k = 0; k < lblocks; ++k)
					oss_gemm(bsize, buffer
					         [local_first_block + j],
					         buffer0[buff0index][k],
					         A[k][region * lblocks + j]);
			}

			if (next_Ri < first_block)
				MPI_Irecv((void *) buffer[local_first_block],
				          doubles_region, MPI_DOUBLE,
				          source, rtag, MPI_COMM_WORLD,
				          &requests[localindex]);
		}
	}

	// Process my rows
	for (size_t Rl = 0; Rl < lblocks; ++Rl) { // Process my blocks

		oss_potrf(bsize, A[Rl][first_block + Rl]); // Diagonal Block Factorization

		// regions first to send them
		for (size_t region = rank + 1; region < wsize; ++region) {
			for (size_t j = 0; j < lblocks; ++j) {
				oss_trsm(bsize, A[Rl][first_block + Rl], A[Rl][region * lblocks + j]);
			}
			for (size_t dest = rank + 1; dest <= region; ++dest) {
				const int rtag = (first_block + Rl) << 16 | region;
				MPI_Request send_request;

				MPI_Isend(A[Rl][region * lblocks], doubles_region, MPI_DOUBLE,
				          dest, rtag, MPI_COMM_WORLD, &send_request);
				MPI_Request_free(&send_request);
			}
		}

		// rest of first row (triangle section)
		for (size_t j = Rl + 1; j < lblocks; ++j) {
			oss_trsm(bsize, A[Rl][first_block + Rl], A[Rl][first_block + j]);
		}

		for (size_t j = Rl + 1; j < lblocks; ++j) {	// local triangle
			for (size_t k = Rl + 1; k < j; ++k) {
				oss_gemm(bsize,
				         A[Rl][first_block + j],
				         A[Rl][first_block + k],
				         A[k][first_block + j]);
			}

			oss_syrk(bsize, A[Rl][first_block + j], A[j][first_block + j]);
		}

		// absolute indices here in j!
		for (size_t j = first_block + lblocks; j < nblocks; ++j) {
			for (size_t k = Rl + 1; k < lblocks; ++k) {
				oss_gemm(bsize, A[Rl][j], A[Rl][first_block + k], A[k][j]);
			}
		}
	}

	if (rank > 0) {
		free(buffer);
		free(buffer0);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
	init_args(argc, argv);

	const size_t DIM = create_cl_int("Dim");	// Matrix size
	const size_t BSIZE = create_cl_int("Bsize");	// Block size
	const bool CHECK = create_optional_cl_int("Check", 0);	// check

	envinfo env;

	Initialize(&env, &argc, &argv, DIM, BSIZE);

	//========= End Command Line ===============
	const size_t dim2 = env.dim * env.dim * sizeof(double);
	const size_t lelements = env.ldim * env.dim;
	const size_t nblocks = env.dim / BSIZE;	// Total number of blocks
	MPI_Request req;

	//=============== Arrays ===================
	double (*matrix)[nblocks][BSIZE][BSIZE] = NULL;
	double (*original)[env.dim] = NULL;
	double (*factorized)[env.dim] = NULL;

	if (env.rank == 0) {

		timer ttimer = create_timer("Total time");
		//======= Allocate matrices ===============
		matrix = malloc(dim2);
		assert(matrix != NULL);

		init_matrix(nblocks, BSIZE, matrix);

		MPI_Iscatter((double *)matrix, lelements, MPI_DOUBLE,
		             MPI_IN_PLACE, lelements, MPI_DOUBLE, 0, MPI_COMM_WORLD, &req);

		if (CHECK) {
			original = malloc(dim2);
			assert(original);

			factorized = malloc(dim2);
			assert(factorized);

			blocked2flat(nblocks, BSIZE, env.dim, matrix, original);
			write_matrix_block("orig.txt", nblocks, BSIZE, matrix);
			write_matrix_flat("orig_flat.txt", env.dim, original);
		}

		MPI_Wait(&req, MPI_STATUS_IGNORE);
		// ===========================================
		printf("# Starting algorithm in process: %d\n", env.rank);
		timer atimer = create_timer("Algorithm time");
		cholesky(nblocks, BSIZE, env.rank, env.worldsize, matrix);
		stop_timer(&atimer);
		// ===========================================

		MPI_Gather(MPI_IN_PLACE, lelements, MPI_DOUBLE,
		           (double *)matrix, lelements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		stop_timer(&ttimer);

		const double performance =
			env.dim * env.dim * env.dim * 3000.0 / getNS_timer(&atimer);

		//======== Check if set =====================
		if (CHECK) {
			blocked2flat(nblocks, BSIZE, env.dim, matrix, factorized);

			printf("# Writing after\n");
			write_matrix_block("fact.txt", nblocks, BSIZE, matrix);
			write_matrix_flat("fact_flat.txt", env.dim, factorized);

			printf("# Checking the correctness of the factorization...\n");
			const double EPS = BLAS_dfpinfo(blas_eps);
			check_factorization(env.dim, original, factorized, env.dim, EPS);

			free(factorized);
			free(original);
		}

	} else {
		matrix = malloc(lelements * sizeof(double));
		assert(matrix);

		MPI_Iscatter(NULL, 0, 0, matrix, lelements, MPI_DOUBLE, 0, MPI_COMM_WORLD, &req);
		MPI_Wait(&req, MPI_STATUS_IGNORE);

		cholesky(nblocks, BSIZE, env.rank, env.worldsize, matrix);

		MPI_Gather(matrix, lelements, MPI_DOUBLE, NULL, 0, 0, 0, MPI_COMM_WORLD);
	}

	free(matrix);
	MPI_Finalize();
	fprintf(stderr, "# Process %d Finalized\n", env.rank);
}
