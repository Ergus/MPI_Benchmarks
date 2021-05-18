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

void cholesky_alloc_init(const size_t nt, double *A[nt][nt], double *Ans[nt][nt],
						 const bool check, const envinfo *env
) {
	#pragma omp parallel
	{
		#pragma omp single
		{
			for (size_t i = 0; i < nt; ++i) {
				for (size_t j = 0; j < nt; ++j) {
					#pragma omp task depend(out: A[i][j]) shared(Ans, A)
					{
						if (check) {
							MPI_Alloc_mem(ts * ts * sizeof(double), MPI_INFO_NULL, &Ans[i][j]);
							block_init(A[i][j], ts, ts);
							if (i == j) {
								Ans[i][i][i * ts + i] = (double) nt;
							}
						}

						if (block_rank[i * nt + j] == env->rank) {
							MPI_Alloc_mem(ts * ts * sizeof(double), MPI_INFO_NULL, &A[i][j]);
							if (check) {
								memcpy(A[i][j], Ans[i][j], ts * ts * sizeof(double));
							} else {
								block_init(A[i][j], ts, ts);
								if (i == j) {
									A[i][i][i * ts + i] = (double) nt;
								}
							}
						} else {
							A[i][j] = NULL;
						}
					} // pragma
				} // for j
			} // for i
		} // omp single
	} // omp parallel
}

void cholesky_single(const size_t nt, double* A[nt][nt])
{
	const size_t ts = env->ts;
	const size_t rank = env->rank;

	for (size_t k = 0; k < nt; k++) {
		#pragma omp task depend(out: A[k][k])
		{
			omp_potrf(A[k][k], ts, ts);
			dbprintf_if(rank == 0, "# potrf:out:A[%d][%d]\n", k, k);
		}

		for (size_t i = k + 1; i < nt; i++) {
			#pragma omp task depend(in: A[k][k]) depend(out: A[k][i])
			{
				omp_trsm(A[k][k], A[k][i], ts, ts);
				dbprintf_if(rank == 0,
							"# trsm :in:A[%d][%d]:out:A[%d][%d]\n", k, k, k, i);
			}
		}
		for (size_t i = k + 1; i < nt; i++) {
			for (size_t j = k + 1; j < i; j++) {
				#pragma omp task depend(in: A[k][i], A[k][j]) depend(out: A[j][i])
				{
					omp_gemm(A[k][i], A[k][j], A[j][i], ts, ts);
					dbprintf_if(rank == 0,
								"# gemm :in:A[%d][%d]:A[%d][%d]:out:A[%d][%d]\n",
								k, i, k, j, j, i);
				}
			}
			#pragma omp task depend(in: A[k][i]) depend(out: A[i][i])
			{
				omp_syrk(A[k][i], A[i][i], ts, ts);
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
	#pragma omp parallel
	{
		#pragma omp single
		{
			const size_t np = env->worldsize;
			const size_t rank = env->rank;
			{ // scope
				for (size_t k = 0; k < nt; k++) {
					// sentinel task to limit communication task parallelism
					if (block_rank[k*nt+k] == rank) {
						#pragma omp task depend(out: A[k][k]) firstprivate(k)
						omp_potrf(A[k][k], ts, ts);
					}

					if (np > 1) { // A[k][k]
						if (block_rank[k*nt+k] == rank) {  // I have this diagonal element.
							#pragma omp task depend(in: A[k][k]) firstprivate(k) untied
							{
								MPI_Request reqs[np];
								int nreqs = 0;
								for (int dst = 0; dst < np; dst++) {
									for (int kk = k+1; kk < nt; kk++) {
										if (dst == block_rank[k*nt+kk]) {
											MPI_Isend(A[k][k], ts*ts, MPI_DOUBLE, dst, k*nt+k, MPI_COMM_WORLD, &reqs[nreqs++]);
											break;
										}
									}
								}
								for (int i = 0; i < nreqs; ++i) {
									wait(&reqs[i]);
								}
							}
						} else {
							for (int i = k + 1; i < nt; i++) {
								if (block_rank[k*nt+i] == rank) {
									#pragma omp task depend(out: B) firstprivate(k) untied
									{
										MPI_Request recv_req;
										MPI_Irecv(B, ts*ts, MPI_DOUBLE, block_rank[k*nt+k], k*nt+k, MPI_COMM_WORLD, &recv_req);
										wait(&recv_req);
									}
									break;
								}
							}
						} // block_rank[k*nt+k] == rank
					} // np != 1 && A[k][k]

					for (size_t i = k + 1; i < nt; i++) {
						if (block_rank[k*nt+i] == rank) {
							double *ptr = (block_rank[k*nt+k] == rank ? A[k][k] : B);
							#pragma omp task depend(in: ptr) depend(out: A[k][i]) firstprivate(k, i)
							omp_trsm(ptr, A[k][i], ts, ts);
						}

						if (np > 1) { // A[k][i]
							if (block_rank[k*nt+i] == rank) {
								char send_flags[np] = {0};

								send_flags[block_rank[i*nt+i]] = 1;

								for (size_t ii = k + 1; ii < i; ii++) {
									if (!send_flags[block_rank[ii*nt+i]]) {
										send_flags[block_rank[ii*nt+i]] = 1;
									}
								}
								for (size_t ii = i + 1; ii < nt; ii++) {
									if (!send_flags[block_rank[i*nt+ii]]) {
										send_flags[block_rank[i*nt+ii]] = 1;
									}
								}

								for (size_t dst = 0; dst < np; dst++) {
									if (send_flags[dst] && dst != rank) {
										#pragma omp task depend(in: A[k][i]) firstprivate(k, i, dst) untied
										{
											MPI_Request send_req;
											MPI_Isend(A[k][i], ts*ts, MPI_DOUBLE, dst, k*nt+i, MPI_COMM_WORLD, &send_req);
											wait(&send_req);
										}
									}
								}
							} else {
								char recv_flag = (block_rank[i*nt+i] == rank);

								if (recv_flag == 0) {
									for (size_t ii = k + 1; ii < i; ++ii) {
										if (block_rank[ii*nt+i] == rank){
											recv_flag = 1;
											break;
										}
									}
								}

								if (recv_flag == 0) {
									for (size_t ii = i + 1; ii < nt; ++ii) {
										if (block_rank[i*nt+ii] == rank) {
											recv_flag = 1;
											break;
										}
									}
								}

								// Decide if I should send
								if (recv_flag != 0) {
									#pragma omp task depend(out: C[i]) firstprivate(k, i) untied
									{
										MPI_Request recv_req;
										MPI_Irecv(C[i], ts*ts, MPI_DOUBLE, block_rank[k*nt+i], k*nt+i, MPI_COMM_WORLD, &recv_req);
										wait(&recv_req);
									}
								}
							} // block_rank[k*nt+i] == rank
						} // np > 1 && A[k][i]
					}


					for (size_t i = k + 1; i < nt; i++) {
						double *ptr1 = (block_rank[k*nt+i] == rank ? A[k][i] : C[i]);

						for (size_t j = k + 1; j < i; j++) {
							if (block_rank[j*nt+i] == rank) {
								double *ptr2 = (block_rank[k*nt+j] == rank ? A[k][j] : C[j]);
								#pragma omp task depend(in: ptr1, ptr2) depend(out: A[j][i]) firstprivate(j, i)
								omp_gemm(ptr1, ptr2, A[j][i], ts, ts);
							}
						}

						if (block_rank[i*nt+i] == rank) {
							#pragma omp task depend(in: ptr) depend(out: A[i][i]) firstprivate(k, i)
							omp_syrk(ptr, A[i][i], ts, ts);
						}
					}
				}
			} // scope
			#pragma omp taskwait
			MPI_Barrier(MPI_COMM_WORLD);

		} // pragma omp single
	} // pragma omp parallel
}


int main(int argc, char *argv[])
{
	init_args(argc, argv);

	const char *PREFIX = basename(argv[0]);
	const int ROWS = create_cl_int("Rows");
	const int TS = create_cl_int("Tasksize");
	const int CHECK = create_optional_cl_int("Print", 0);

	envinfo env;

	Initialize(&env, &argc, &argv, ROWS, TS);

	printf("# Initializing data in process: %d\n", env.rank);;

	timer ttimer = create_timer("Total_time");

	// Allocate memory
	const size_t nt = ROWS / TS;

	int *block_rank = malloc(nt * nt * sizeof(int));
	get_block_rank(nt, block_rank, env);

	double *A[nt][nt], *B, *C[nt], *Ans[nt][nt];

	cholesky_alloc_init(nt, A, Ans, CHECK, env);

	MPI_Alloc_mem(ts * ts * sizeof(double), MPI_INFO_NULL, &B);
	for (size_t i = 0; i < nt; i++) {
		MPI_Alloc_mem(ts * ts * sizeof(double), MPI_INFO_NULL, &C[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// ===========================================
	printf("# Starting algorithm in process: %d\n", env.rank);

	timer atimer = create_timer("Algorithm_time");
	cholesky_mpi(ts, nt, (double* (*)[nt])A, B, C, block_rank);
	MPI_Barrier(MPI_COMM_WORLD);

	stop_timer(&atimer);
	// ===========================================

	printf("# Finished algorithm in process: %d\n", env.rank);

	if (CHECK) {
		timer stimer = create_timer("Single_time");
		cholesky_single(ts, nt, (double* (*)[nt]) Ans);
		stop_timer(&stimer);
	}

	stop_timer(&ttimer);

	/* Verification */
	if (CHECK) {
		for (size_t i = 0; i < nt; i++) {
			for (size_t j = 0; j < nt; j++) {
				if (block_rank[i * nt + j] == mype) {
					for (size_t k = 0; k < ts*ts; k++) {
						if (Ans[i][j][k] != A[i][j][k])
							check = 2;
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
			if (block_rank[i*nt+j] == mype) {
				free(A[i][j]);
			}
			if (check) {
				free(Ans[i][j]);
			}
		}
		free(C[i]);
	}
	free(B);
	free(block_rank);
	free_args();

	Finalize();
	return 0;
}
