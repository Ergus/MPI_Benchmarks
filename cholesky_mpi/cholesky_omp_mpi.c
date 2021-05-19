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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.	 If not, see <http://www.gnu.org/licenses/>.
 */

#include "cholesky_omp_mpi.h"

void get_block_rank(const size_t nt, int *block_rank, const envinfo *env)
{
	const size_t np = env->worldsize;

	int row = np, col = np;

	if (np != 1) {
		while (1) {
			row = row / 2;
			if (row * col == np) {
				break;
			}
			col = col / 2;
			if (row * col == np) {
				break;
			}
		}
	}
	dbprintf_if(env->rank == 0, "# row = %d, col = %d\n", row, col);

	size_t tmp_rank = 0, offset = 0;
	for (size_t i = 0; i < nt; i++) {
		for (size_t j = 0; j < nt; j++) {
			block_rank[i * nt + j] = tmp_rank + offset;
			++tmp_rank;
			if (tmp_rank >= col)
				tmp_rank = 0;
		}
		tmp_rank = 0;
		offset = (offset + col >= np) ? 0 : offset + col;
	}
}


//! This is to initialize non-diagonal blocks (assuming 1 block/task)
void fill_block(const size_t ts, double *matrix,
                const size_t i, const size_t j, const envinfo *env
) {
	assert(matrix);

	const size_t seed = (i > j) ? i * ts + j : j * ts + i;
	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(seed, &status);
	double rnd1, rnd2;

	for (size_t k = 0; k < ts; ++k) {
		if (i == j) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);
			matrix[k * ts + k] = rnd1 * rnd2 + env->dim;

			for (size_t l = k + 1; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				const double val = rnd1 * rnd2;
				matrix[k * ts + l] = val;
				matrix[l * ts + k] = val;
			}

		} else {
			for (size_t l = 0; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				if (i > j) {
					matrix[k * ts + l] = rnd1 * rnd2;
				} else {
					matrix[l * ts + k] = rnd1 * rnd2;
				}
			}
		}
	}
}


void cholesky_alloc_init(const size_t nt, double *A[nt][nt], double *Ans[nt][nt],
                         const bool check,
                         const int *block_rank, const envinfo *env
) {
	#pragma omp parallel
	{
		#pragma omp single
		{
			const size_t ts = env->ts;

			for (size_t i = 0; i < nt; ++i) {
				for (size_t j = 0; j < nt; ++j) {

					#pragma omp task depend(out: A[i][j]) shared(Ans, A) firstprivate(i, j)
					{
						if (check) {
							Ans[i][j] = (double *) malloc(ts * ts * sizeof(double));
							fill_block(ts, Ans[i][j], i, j, env);
						}

						if (block_rank[i * nt + j] == env->rank) {
							A[i][j] = (double *) malloc(ts * ts * sizeof(double));
							if (check) {
								memcpy(A[i][j], Ans[i][j], ts * ts * sizeof(double));
							} else {
								fill_block(ts, A[i][j], i, j, env);
							}
						} else {
							A[i][j] = NULL;
						}
					} // pragma
				} // for j
			} // for i
		} // omp single

		#pragma omp taskwait
	} // omp parallel
}

void cholesky_single(const size_t nt, double* A[nt][nt], const envinfo *env)
{
	const size_t ts = env->ts;
	const size_t rank = env->rank;

	for (size_t k = 0; k < nt; k++) {
		#pragma omp task depend(out: A[k][k])
		{
			omp_potrf(ts, A[k][k]);
			dbprintf_if(rank == 0, "# potrf:out:A[%d][%d]\n", k, k);
		}

		for (size_t i = k + 1; i < nt; i++) {
			#pragma omp task depend(in: A[k][k]) depend(out: A[k][i])
			{
				omp_trsm(ts, A[k][k], A[k][i]);
				dbprintf_if(rank == 0,
							"# trsm :in:A[%d][%d]:out:A[%d][%d]\n", k, k, k, i);
			}
		}
		for (size_t i = k + 1; i < nt; i++) {
			for (size_t j = k + 1; j < i; j++) {
				#pragma omp task depend(in: A[k][i], A[k][j]) depend(out: A[j][i])
				{
					omp_gemm(ts, A[k][i], A[k][j], A[j][i]);
					dbprintf_if(rank == 0,
								"# gemm :in:A[%d][%d]:A[%d][%d]:out:A[%d][%d]\n",
								k, i, k, j, j, i);
				}
			}
			#pragma omp task depend(in: A[k][i]) depend(out: A[i][i])
			{
				omp_syrk(ts, A[k][i], A[i][i]);
				dbprintf_if(rank == 0,
				            "# syrk :in:A[%d][%d]:out:A[%d][%d]\n",
				            k, i, i, i);
			}
		}
	}
	#pragma omp taskwait
}


void cholesky_mpi(const size_t nt, double *A[nt][nt], double *B, double *C[nt],
				  const int *block_rank, const envinfo *env
) {
	assert(A != NULL);
	assert(B != NULL);
	assert(C != NULL);

	#pragma omp parallel
	{
		#pragma omp single
		{
			const size_t np = env->worldsize;
			const size_t rank = env->rank;
			const size_t ts = env->ts;

			{ // scope
				for (size_t k = 0; k < nt; ++k) {
					// sentinel task to limit communication task parallelism
					if (block_rank[k * nt + k] == rank) {
						#pragma omp task depend(out: A[k][k]) firstprivate(k)
						omp_potrf(ts, A[k][k]);
					}

					if (np > 1) { // A[k][k]
						const int tag = k * nt + k;

						if (block_rank[tag] == rank) {  // I have this diagonal element.

							#pragma omp task depend(in: A[k][k]) firstprivate(k, tag) untied
							{
								MPI_Request reqs[np];
								int nreqs = 0;

								for (size_t dst = 0; dst < np; ++dst) {
									if (dst == rank) {
										continue;
									}
									for (size_t kk = k + 1; kk < nt; ++kk) {
										if (dst == block_rank[k * nt + kk]) {
											MPI_Isend(A[k][k], ts * ts, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &reqs[nreqs++]);
											break;
										}
									}
								}
								for (size_t i = 0; i < nreqs; ++i) {
									wait(&reqs[i]);
								}
							}
						} else {
							for (size_t i = k + 1; i < nt; ++i) {
								if (block_rank[k * nt + i] == rank) {
									#pragma omp task depend(out: B) firstprivate(k, tag) untied
									{
										MPI_Request recv_req;
										//printf("Process %zu: Waiting0 %d from %d\n", rank, tag, block_rank[tag]);
										MPI_Irecv(B, ts * ts, MPI_DOUBLE,
										          block_rank[tag], tag,
										          MPI_COMM_WORLD, &recv_req);
										wait(&recv_req);
										//printf("Process %zu: Received0 %d from %d\n", rank, tag, block_rank[tag]);
									}
									break;
								}
							}
						} // block_rank[k*nt+k] == rank
					} // np != 1 && A[k][k]

					for (size_t i = k + 1; i < nt; ++i) {
						if (block_rank[k * nt + i] == rank) {
							double *ptr = (block_rank[k * nt + k] == rank ? A[k][k] : B);
							#pragma omp task depend(in: ptr) depend(out: A[k][i]) firstprivate(k, i)
							omp_trsm(ts, ptr, A[k][i]);
						}

						if (np > 1) { // A[k][i]
							const int tag = k * nt + i;

							if (block_rank[k * nt + i] == rank) {
								char send_flags[np];
								memset(send_flags, 0, np * sizeof(char));

								send_flags[block_rank[i * nt + i]] = 1;

								for (size_t ii = k + 1; ii < i; ii++) {
									if (!send_flags[block_rank[ii * nt + i]]) {
										send_flags[block_rank[ii * nt + i]] = 1;
									}
								}
								for (size_t ii = i + 1; ii < nt; ii++) {
									if (!send_flags[block_rank[i * nt + ii]]) {
										send_flags[block_rank[i * nt + ii]] = 1;
									}
								}

								for (size_t dst = 0; dst < np; ++dst) {
									if (send_flags[dst] && dst != rank) {
										#pragma omp task depend(in: A[k][i]) firstprivate(k, i, dst, tag) untied
										{
											//printf("Process %zu: Sending1 %d -> %zu\n", rank, tag, dst);
											MPI_Request send_req;
											MPI_Isend(A[k][i], ts * ts, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, &send_req);
											wait(&send_req);
											//printf("Process %zu: Sent1 %d -> %zu\n", rank, tag, dst);
										}
									}
								}
							} else {
								char recv_flag = (block_rank[i * nt + i] == rank);

								if (recv_flag == 0) {
									for (size_t ii = k + 1; ii < i; ++ii) {
										if (block_rank[ii * nt + i] == rank){
											recv_flag = 1;
											break;
										}
									}
								}

								if (recv_flag == 0) {
									for (size_t ii = i + 1; ii < nt; ++ii) {
										if (block_rank[i * nt + ii] == rank) {
											recv_flag = 1;
											break;
										}
									}
								}

								// Decide if I should send
								if (recv_flag != 0) {
									#pragma omp task depend(out: C[i]) firstprivate(k, i, tag) untied
									{
										MPI_Request recv_req;
										MPI_Irecv(C[i], ts * ts, MPI_DOUBLE, block_rank[tag], tag, MPI_COMM_WORLD, &recv_req);
										wait(&recv_req);
									}
								}
							} // block_rank[k*nt+i] == rank
						} // np > 1 && A[k][i]
					}


					for (size_t i = k + 1; i < nt; i++) {
						double *ptr1 = (block_rank[k * nt + i] == rank ? A[k][i] : C[i]);

						for (size_t j = k + 1; j < i; j++) {
							if (block_rank[j*nt+i] == rank) {
								double *ptr2 = (block_rank[k * nt + j] == rank ? A[k][j] : C[j]);
								#pragma omp task depend(in: ptr1, ptr2) depend(out: A[j][i]) firstprivate(j, i)
								omp_gemm(ts, ptr1, ptr2, A[j][i]);
							}
						}

						if (block_rank[i*nt+i] == rank) {
							#pragma omp task depend(in: ptr1) depend(out: A[i][i]) firstprivate(k, i)
							omp_syrk(ts, ptr1, A[i][i]);
						}
					} // for i
				} // for k
			} // scope
		} // pragma omp single

		#pragma omp taskwait
		MPI_Barrier(MPI_COMM_WORLD);

	} // pragma omp parallel
}


int main(int argc, char *argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const int ROWS = create_cl_int("Rows");
	const int TS = create_cl_int("Tasksize");
	const int CHECK = create_optional_cl_int("Check", 0);

	envinfo env;

	Initialize(&env, &argc, &argv, ROWS, TS);

	printf("# Initializing data in process: %d\n", env.rank);;

	timer ttimer = create_timer("Total_time");

	// Allocate memory
	const size_t nt = env.dim / env.ts;

	int *block_rank = malloc(nt * nt * sizeof(int));
	get_block_rank(nt, block_rank, &env);

	double *A[nt][nt], *B, *C[nt], *Ans[nt][nt];

	cholesky_alloc_init(nt, A, Ans, CHECK, block_rank, &env);

	MPI_Alloc_mem(TS * TS * sizeof(double), MPI_INFO_NULL, &B);
	for (size_t i = 0; i < nt; i++) {
		MPI_Alloc_mem(TS * TS * sizeof(double), MPI_INFO_NULL, &C[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");
	cholesky_mpi(nt, (double* (*)[nt])A, B, C, block_rank, &env);
	MPI_Barrier(MPI_COMM_WORLD);

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm in process: %d\n", env.rank);

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_single(nt, (double* (*)[nt])Ans, &env);
		stop_timer(&stimer);
	}

	stop_timer(&ttimer);

	// Verification
	if (CHECK) {
		bool match = true;
		for (size_t i = 0; i < nt && match; i++) {
			for (size_t j = 0; j < nt && match; j++) {
				if (block_rank[i * nt + j] == env.rank) {
					for (size_t k = 0; k < env.ts * env.ts; ++k) {
						if (A[i][j][k] != Ans[i][j][k]){
							match = false;
							fprintf(stderr,
							        "Check failed in A[%zu + %zu][%zu + %zu]=%lf expected: %lf rank %d\n",
							        i, k / env.ts, j, k % env.ts,
							        A[i][j][k], Ans[i][j][k], env.rank);
							break;
						}
					}
				}
			}
		}
	}

	if (env.rank == 0) {
		create_reportable_int("worldsize", env.worldsize);
		create_reportable_int("cpu_count", env.cpu_count);
		report_args();
	}

	// Release memory
	for (size_t i = 0; i < nt; ++i) {
		for (size_t j = 0; j < nt; j++) {
			if (block_rank[i * nt + j] == env.rank) {
				assert(A[i][j] != NULL);
				free(A[i][j]);
			}
			if (CHECK) {
				assert(Ans[i][j] != NULL);
				free(Ans[i][j]);
			}
		}
		assert(C[i] != NULL);
		free(C[i]);
	}
	assert(B);
	free(B);

	free(block_rank);
	free_args();

	Finalize();
	return 0;
}
