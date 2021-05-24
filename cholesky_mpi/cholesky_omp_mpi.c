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

void print_matrix(const size_t nt, const size_t ts,
                  double A[nt][nt][ts][ts], const envinfo *env)
{
	MPI_Barrier(MPI_COMM_WORLD);

	for (size_t it = 0; it < env->worldsize; ++it) {
		if (it == env->rank) {
			printf("Rank %d\n", env->rank);
			for (size_t i = 0; i < nt; i++) {
				for (size_t k = 0; k < ts; ++k) {
					for (size_t j = 0; j < nt; ++j) {
						for (size_t l = 0; l < ts; ++l) {
							printf("%5.2f ", (float)A[i][j][k][l]);
						}
						printf(" ");
					}
					printf("\n");
				}
				printf(" \n");
			}
			printf("--------\n");
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


void print_block(const size_t ts, double A[ts][ts], const envinfo *env)
{
	printf("Rank %d\n", env->rank);

	for (size_t i = 0; i < ts; ++i) {
		for (size_t j = 0; j < ts; ++j) {
			printf("%5.2f ", (float)A[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
}



void get_block_rank(const size_t nt, int block_rank[nt][nt], const envinfo *env)
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
			block_rank[i][j] = tmp_rank + offset;
			++tmp_rank;
			if (tmp_rank >= col)
				tmp_rank = 0;
		}
		tmp_rank = 0;
		offset = (offset + col >= np) ? 0 : offset + col;
	}
}


//! This is to initialize non-diagonal blocks (assuming 1 block/task)
void fill_block(const size_t ts, double block[ts][ts],
                const size_t i, const size_t j, const envinfo *env
) {
	assert(block);

	const size_t seed = (i > j) ? i * ts + j : j * ts + i;
	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(seed, &status);
	double rnd1, rnd2;

	for (size_t k = 0; k < ts; ++k) {
		if (i == j) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);
			block[k][k] = rnd1 * rnd2 + env->dim;

			for (size_t l = k + 1; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				const double val = rnd1 * rnd2;
				block[k][l] = val;
				block[l][k] = val;
			}

		} else {
			for (size_t l = 0; l < ts; ++l) {
				drand48_r(&status, &rnd1);
				drand48_r(&status, &rnd2);

				if (i > j) {
					block[k][l] = rnd1 * rnd2;
				} else {
					block[l][k] = rnd1 * rnd2;
				}
			}
		}
	}
}


void cholesky_init(const size_t nt, const size_t ts,
                   double A[nt][nt][ts][ts], double Ans[nt][nt][ts][ts],
                   const bool check,
                   int block_rank[nt][nt], const envinfo *env
) {
	#pragma omp parallel
	{
		#pragma omp single
		{
			const size_t ts = env->ts;

			for (size_t i = 0; i < nt; ++i) {
				for (size_t j = 0; j < nt; ++j) {

					#pragma omp task depend(out: A[i][j][0:ts][0:ts])	\
						shared(Ans, A) firstprivate(i, j)
					{
						if (check) {
							fill_block(ts, Ans[i][j], i, j, env);
						}

						if (block_rank[i][j] == env->rank) {
							if (check) {
								memcpy(A[i][j], Ans[i][j], ts * ts * sizeof(double));
							} else {
								fill_block(ts, A[i][j], i, j, env);
							}
						}
					} // pragma
				} // for j
			} // for i
		} // omp single

		#pragma omp taskwait
	} // omp parallel
}

void cholesky_single(const size_t nt, const size_t ts,
                     double A[nt][nt][ts][ts], const envinfo *env)
{
	const size_t rank = env->rank;
	#pragma omp parallel
	{
		#pragma omp single
		{

			for (size_t k = 0; k < nt; ++k) {
				#pragma omp task depend(inout: A[k][k]) firstprivate(k)
				omp_potrf(ts, A[k][k]);

				for (size_t i = k + 1; i < nt; ++i) {
					#pragma omp task depend(in: A[k][k]) \
						depend(out: A[k][i]) firstprivate(k, i)
					omp_trsm(ts, A[k][k], A[k][i]);
				}
				for (size_t i = k + 1; i < nt; ++i) {
					for (size_t j = k + 1; j < i; ++j) {
						#pragma omp task depend(in: A[k][i], A[k][j]) \
							depend(out: A[j][i]) firstprivate(k, i, j)
						omp_gemm(ts, A[k][i], A[k][j], A[j][i]);
					}
					#pragma omp task depend(in: A[k][i]) \
						depend(out: A[i][i]) firstprivate(k, i)
					omp_syrk(ts, A[k][i], A[i][i]);
				} // for i
				#pragma omp taskwait
			} // for k
		} // pragma omp single
	}  // pragma omp parallel
}


void cholesky_mpi(const size_t nt, const size_t ts,
                  double A[nt][nt][ts][ts], double B[ts][ts], double C[nt][ts][ts],
				  int block_rank[nt][nt], const envinfo *env
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

			for (size_t k = 0; k < nt; ++k) {
				// sentinel task to limit communication task parallelism
				if (block_rank[k][k] == rank) {
					#pragma omp task depend(inout: A[k][k]) firstprivate(k)
					omp_potrf(ts, A[k][k]);
				}

				if (np > 1) { // A[k][k]
					const int tag = k * nt + k;

					if (block_rank[k][k] == rank) {  // I have this diagonal element.

						#pragma omp task depend(in: A[k][k]) firstprivate(k, tag) untied
						{
							MPI_Request reqs[np];
							size_t nreqs = 0;

							for (size_t dst = 0; dst < np; ++dst) {
								if (dst == rank) {
									continue;
								}
								for (size_t kk = k + 1; kk < nt; ++kk) {
									if (dst == block_rank[k][kk]) {
										MPI_Isend(A[k][k], ts * ts, MPI_DOUBLE,
										          dst, tag,
										          MPI_COMM_WORLD, &reqs[nreqs]);
										++nreqs;
										break;
									}
								}
							}
							for (size_t i = 0; i < nreqs; ++i) {
								wait(&reqs[i]);
							}
						}
					} else { // block_rank[k][k] != rank
						for (size_t i = k + 1; i < nt; ++i) {
							if (block_rank[k][i] == rank) {
								#pragma omp task depend(out: B[0]) firstprivate(k, tag) untied
								{
									MPI_Request recv_req;
									MPI_Irecv(B, ts * ts, MPI_DOUBLE,
									          block_rank[k][k], tag,
									          MPI_COMM_WORLD, &recv_req);
									wait(&recv_req);
								}
								break;
							}
						}
					} // block_rank[k*nt+k] == rank
				} // np != 1 && A[k][k]

				MPI_Barrier(MPI_COMM_WORLD);

				for (size_t i = k + 1; i < nt; ++i) {
					if (block_rank[k][i] == rank) {
						double (*ptr)[ts] = (block_rank[k][k] == rank ? A[k][k] : B);
						#pragma omp task depend(in: ptr[0]) \
							depend(out: A[k][i]) firstprivate(k, i, ptr)
						omp_trsm(ts, ptr, A[k][i]);
					}

					if (np > 1) { // A[k][i]
						const int tag = k * nt + i;

						if (block_rank[k][i] == rank) {
							char send_flags[np];
							memset(send_flags, 0, np * sizeof(char));

							send_flags[block_rank[i][i]] = 1;

							for (size_t ii = k + 1; ii < i; ++ii) {
								if (!send_flags[block_rank[ii][i]]) {
									send_flags[block_rank[ii][i]] = 1;
								}
							}
							for (size_t ii = i + 1; ii < nt; ++ii) {
								if (!send_flags[block_rank[i][ii]]) {
									send_flags[block_rank[i][ii]] = 1;
								}
							}

							for (size_t dst = 0; dst < np; ++dst) {
								if (send_flags[dst] && dst != rank) {
									#pragma omp task depend(in: A[k][i]) \
										firstprivate(k, i, dst, tag) untied
									{
										MPI_Request send_req;
										MPI_Isend(A[k][i], ts * ts, MPI_DOUBLE,
										          dst, tag,
										          MPI_COMM_WORLD, &send_req);
										wait(&send_req);
									}
								}
							}
						} else {
							char recv_flag = (block_rank[i][i] == rank);

							if (recv_flag == 0) {
								for (size_t ii = k + 1; ii < i; ++ii) {
									if (block_rank[ii][i] == rank){
										recv_flag = 1;
										break;
									}
								}
							}

							if (recv_flag == 0) {
								for (size_t ii = i + 1; ii < nt; ++ii) {
									if (block_rank[i][ii] == rank) {
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
									MPI_Irecv(C[i], ts * ts, MPI_DOUBLE,
									          block_rank[k][i], tag,
									          MPI_COMM_WORLD, &recv_req);
									wait(&recv_req);
								}
							}
						} // block_rank[k*nt+i] == rank
					} // np > 1 && A[k][i]
				}


				for (size_t i = k + 1; i < nt; ++i) {
					double (*ptr1)[ts] = (block_rank[k][i] == rank ? A[k][i] : C[i]);

					for (size_t j = k + 1; j < i; ++j) {
						if (block_rank[j][i] == rank) {
							double (*ptr2)[ts] = (block_rank[k][j] == rank ? A[k][j] : C[j]);
							#pragma omp task depend(in: ptr1[0], ptr2[0]) \
								depend(out: A[j][i]) firstprivate(i, j, k, ptr1, ptr2)
							omp_gemm(ts, ptr1, ptr2, A[j][i]);
						}
					}

					if (block_rank[i][i] == rank) {
						#pragma omp task depend(in: ptr1[0]) \
							depend(out: A[i][i]) firstprivate(k, i, ptr1)
						omp_syrk(ts, ptr1, A[i][i]);
					}
				} // for i

				#pragma omp taskwait

			} // for k
		} // pragma omp single

	} // pragma omp parallel
	MPI_Barrier(MPI_COMM_WORLD);
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

	int (*block_rank)[nt] = malloc(nt * nt * sizeof(int));
	get_block_rank(nt, block_rank, &env);

	double (*A)[nt][TS][TS] = malloc(ROWS * ROWS * sizeof(double)), // [nt][nt][TS][TS]
		(*B)[TS] = malloc(TS * TS * sizeof(double)),                // [TS][TS]
		(*C)[TS][TS] = malloc(ROWS * TS * sizeof(double)),          // [nt][TS][TS]
		(*Ans)[nt][TS][TS] = (CHECK ? malloc(ROWS * ROWS * sizeof(double)) : NULL);


	cholesky_init(nt, TS, A, Ans, CHECK, block_rank, &env);

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");
	cholesky_mpi(nt, TS, A, B, C, block_rank, &env);
	MPI_Barrier(MPI_COMM_WORLD);

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm in process: %d\n", env.rank);

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_single(nt, TS, Ans, &env);
		stop_timer(&stimer);
	}

	stop_timer(&ttimer);

	MPI_Barrier(MPI_COMM_WORLD);

	// Verification
	if (CHECK && env.rank == 1) {
		bool match = true;
		for (size_t i = 0; i < nt && match; i++) {
			for (size_t j = 0; j < nt && match; j++) {
				if (block_rank[i][j] == env.rank) {
					for (size_t k = 0; k < env.ts; ++k) {
						for (size_t l = 0; l < env.ts; ++l) {
							if (A[i][j][k][l] != Ans[i][j][k][l]) {
								match = false;
								fprintf(stdout,
								        "Check failed in A[%zu][%zu][%zu][%zu]=%lf expected: %lf rank %d\n",
								        i, j, k, l,
								        A[i][j][k][l], Ans[i][j][k][l],
								        env.rank);
								break;
							}
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
	if (CHECK) {
		assert(Ans);
		free(Ans);
	}

	assert(C);
	free(C);

	assert(B);
	free(B);

	assert(A);
	free(A);

	free(block_rank);
	free_args();

	Finalize();
	return 0;
}
